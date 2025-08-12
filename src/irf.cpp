#include <RcppArmadillo.h>
#include <random>
#include <thread>
#include <Rmath.h>
#include <R_ext/Random.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

void get_PHI(arma::cube& PHI, arma::cube& Fmat, const int nhor) {
    const int plag = Fmat.n_slices;
    const int bigK = Fmat.n_cols;
    
    arma::cube PHIx(bigK, bigK, plag+nhor+1, arma::fill::zeros);
    PHIx.slice(plag) = arma::mat(bigK,bigK, arma::fill::eye);
    for(int ihor=plag; ihor < plag+nhor; ihor++) {
        arma::mat acc(bigK,bigK, arma::fill::zeros);
        for(int pp=0; pp<plag; pp++) {
            arma::mat Fmatslice = Fmat.slice(pp);
            arma::mat PHIxslice = PHIx.slice(ihor-pp);
            acc = acc + Fmatslice*PHIxslice;
        }
        PHIx.slice(ihor+1) = acc;
    }
    PHI = PHIx.slices(plag,plag+nhor);
}

void get_nullspace(arma::mat& NU, arma::mat& M) {
    arma::mat Q, R; 
    qr(Q, R, M);
    arma::uword r = rank(M);
    if(r == 0) {
        arma::uvec set(M.n_cols); 
        for(size_t c = 0; c < M.n_cols; c++) {
            set(c) = c;
        }
        NU = Q.cols(set);
    } else {
        arma::uvec set(r); 
        for(size_t c = 0; c < r; c++) {
            set(c) = c;
        }
        NU = Q;
        NU.shed_cols(set);
    }
}

//' @name compute_irf
//' @noRd
// [[Rcpp::export]]
Rcpp::List compute_irf(arma::cube A_large, 
                      arma::cube S_large, 
                      arma::cube Ginv_large, 
                      const int type, 
                      const int nhor, 
                      const int thindraws, 
                      const SEXP shocklist_in, 
                      const bool save_rot, 
                      const bool verbose = true) {
    
    // Input processing
    Rcpp::List shocklist(Rcpp::clone(shocklist_in));
    
    // Get parameters from shocklist
    Rcpp::List shock_idx = shocklist["shock.idx"];
    Rcpp::LogicalVector shock_cidx = shocklist["shock.cidx"];
    const int plag = shocklist["plag"];
    const int bigK = A_large.n_rows;
    const int k = A_large.n_cols;
    const int N = shock_idx.size();
    
    // Prepare output containers
    arma::field<arma::cube> irf_output(thindraws);
    irf_output.fill(arma::ones<arma::cube>(bigK,bigK,nhor));
    
    arma::uword size_of_rot = save_rot ? thindraws : 0;
    arma::field<arma::mat> rot_output(size_of_rot);
    if(save_rot) {
        rot_output.fill(arma::ones<arma::mat>(bigK,bigK));
    }
    
    arma::vec counter_output(thindraws);
    
    // OpenMP Setup
    #ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int use_threads = std::min(num_threads, max_threads);
    omp_set_num_threads(use_threads);
    if (verbose) {
        Rcpp::Rcout << "Using " << use_threads << " threads for parallel computation\n";
    }
    #endif
    
    // Main parallel computation
    #pragma omp parallel for shared(irf_output, rot_output, counter_output, A_large, S_large, Ginv_large, shock_idx, shock_cidx) schedule(dynamic)
    for(int irep = 0; irep < thindraws; irep++) {
        // Thread-local variables
        arma::cube irfa(bigK, bigK, nhor);
        arma::cube PHI(bigK, bigK, nhor);
        arma::mat Amat = A_large.slice(irep);
        arma::mat Smat = S_large.slice(irep);
        arma::mat Ginv = Ginv_large.slice(irep);
        arma::cube Fmat(bigK, bigK, plag, arma::fill::zeros);
        arma::mat P0G(bigK, bigK, arma::fill::zeros);
        arma::mat Q_bar(bigK, bigK, arma::fill::eye);
        arma::mat invGSigma_u(bigK, bigK, arma::fill::zeros);
        
        // Thread-local RNG
        std::mt19937 gen(std::random_device{}() + omp_get_thread_num());
        std::normal_distribution<double> rnorm(0.0, 1.0);
        std::uniform_real_distribution<double> runif(0.0, 1.0);
        
        // Construct Fmat
        for(int pp = 0; pp < plag; pp++) {
            Fmat.slice(pp) = Amat.cols(pp*bigK,(pp+1)*bigK-1);
        }
        
        // Create P0G
        bool chol_success = true;
        for(int cc = 0; cc < N; cc++) {
            arma::uvec idx = shock_idx[cc];
            if(shock_cidx[cc]) {
                arma::mat Sig_chol;
                arma::mat sigma = Smat.submat(idx,idx);
                chol_success = chol(Sig_chol, sigma, "lower");
                if(!chol_success) {
                    break;
                }
                P0G.submat(idx,idx) = Sig_chol;
            } else {
                P0G.submat(idx,idx) = Smat.submat(idx,idx);
            }
        }
        
        // Handle Cholesky failure
        if(!chol_success) {
            #pragma omp critical
            {
                counter_output(irep) = Rcpp::as<int>(shocklist["MaxTries"]);
            }
            continue;
        }
        
        // Create PHI
        get_PHI(PHI, Fmat, nhor);
        
        // Type 3 specific computations
        int icounter = 1;
        if(type == 3) {
            int MaxTries = shocklist["MaxTries"];
            arma::mat S_cube = shocklist["S.cube"];
            arma::mat P_cube = shocklist["P.cube"];
            arma::cube Z_cube = shocklist["Z.cube"];
            arma::uvec shock_order = shocklist["shock.order"];
            arma::uvec shock_horz = shocklist["shock.horz"];
            Rcpp::LogicalVector nozero = shocklist["no.zero.restr"];
            
            int H_restr = shock_horz.n_elem;
            int N_restr = bigK * H_restr;
            
            // Build irf_restr
            arma::mat irf_restr(0, bigK);
            for(int hh = 0; hh < H_restr; hh++) {
                int horz = shock_horz(hh);
                arma::mat irf_hh = PHI.slice(horz) * Ginv * P0G;
                irf_restr = join_cols(irf_restr, irf_hh);
            }
            
            // Sort Zcube
            arma::cube Z_cube_sorted(size(Z_cube), arma::fill::zeros);
            for(int i = 0; i < bigK; i++) {
                int idx = shock_order(i);
                Z_cube_sorted.slice(i) = Z_cube.slice(idx);
            }
            
            // Main iteration loop
            double condall = 0.0;
            arma::vec signCheck(bigK, arma::fill::zeros);
            
            while(condall == 0.0 && icounter < MaxTries) {
                arma::mat Q(bigK, bigK, arma::fill::eye);
                
                for(int cc = 0; cc < N; cc++) {
                    arma::uvec idx = shock_idx[cc];
                    int Kidx = idx.size();
                    arma::mat randMat(Kidx, Kidx), Qc(Kidx, Kidx, arma::fill::zeros);
                    
                    if(shock_cidx[cc]) {
                        // Generate random matrix using thread-local RNG
                        for(int kk = 0; kk < Kidx; kk++) {
                            for(int jj = 0; jj < Kidx; jj++) {
                                randMat(jj,kk) = rnorm(gen);
                            }
                        }
                        
                        if(nozero[cc]) {
                            arma::mat R;
                            qr(Qc, R, randMat);
                        } else {
                            for(int i = 0; i < Kidx; i++) {
                                arma::mat Ztemp = Z_cube_sorted.slice(idx[i]);
                                arma::colvec Zsum = sum(abs(Ztemp),1);
                                arma::uvec zidx = find(Zsum);
                                Ztemp = Ztemp.rows(zidx);
                                
                                arma::mat R(0,Kidx);
                                if(i == 0) {
                                    arma::mat R1 = Ztemp * irf_restr.cols(idx);
                                    R = join_cols(R,R1);
                                } else {
                                    arma::mat R2 = Qc.head_cols(i).t();
                                    if(Ztemp.n_elem != 0) {
                                        arma::mat R1 = Ztemp * irf_restr.cols(idx);
                                        R = join_cols(R, R1);
                                    }
                                    R = join_cols(R, R2);
                                }
                                
                                arma::mat Rt = R.t();
                                arma::mat NU;
                                get_nullspace(NU, Rt);
                                arma::vec x_j = randMat.col(i);
                                double div = arma::as_scalar((NU.t() * x_j).t() * (NU.t() * x_j));
                                arma::vec q_j = NU * (NU.t() * x_j / sqrt(div));
                                Qc.col(i) = q_j;
                            }
                        }
                        Q.submat(idx,idx) = Qc;
                    }
                }
                
                // Reorder Q
                for(int i = 0; i < bigK; i++) {
                    int idx = shock_order(i);
                    Q_bar.col(idx) = Q.col(i);
                }
                
                // Check IRF
                arma::mat irf_check = irf_restr * Q_bar;
                for(int kk = 0; kk < bigK; kk++) {
                    arma::vec STemp = S_cube.col(kk);
                    arma::vec PTemp = P_cube.col(kk);
                    if(sum(abs(STemp)) > 0) {
                        arma::vec prob(N_restr);
                        for(int nn = 0; nn < N_restr; nn++) {
                            prob(nn) = (PTemp(nn) > runif(gen)) ? 1 : 0;
                        }
                        arma::mat PDiag = diagmat(prob);
                        arma::vec IrfTemp = sign(irf_check.col(kk));
                        double getsum = arma::as_scalar(IrfTemp.t() * PDiag * STemp);
                        double chksum = arma::as_scalar(sum(abs(PDiag * STemp)));
                        signCheck(kk) = (getsum == chksum) ? 1 : 0;
                    } else {
                        signCheck(kk) = 1;
                    }
                }
                condall = prod(signCheck);
                icounter++;
            }
        }
        
        // Compute shock and impulse responses
        invGSigma_u = Ginv * P0G * Q_bar;
        for(int ihor = 0; ihor < nhor; ihor++) {
            irfa.slice(ihor) = PHI.slice(ihor) * invGSigma_u;
        }
        
        // Save results with proper synchronization
        #pragma omp critical
        {
            irf_output(irep) = irfa;
            if(save_rot) {
                rot_output(irep) = Q_bar;
            }
            counter_output(irep) = icounter;
            
            // Check for user interruption less frequently
            if(irep % 50 == 0) {
                Rcpp::checkUserInterrupt();
            }
        }
    }
    
    return Rcpp::List::create(
        Rcpp::Named("irf", irf_output),
        Rcpp::Named("rot", rot_output),
        Rcpp::Named("counter", counter_output)
    );
}
