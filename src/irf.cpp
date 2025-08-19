#include <RcppArmadillo.h>
#include <random>
#include <thread>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <omp.h>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
// [[Rcpp::plugins(openmp)]]

void get_PHI(arma::cube& PHI, arma::cube& Fmat, const int nhor){
  const int plag = Fmat.n_slices;
  const int bigK = Fmat.n_cols;
  
  arma::cube PHIx(bigK, bigK, plag+nhor+1, arma::fill::zeros);
  PHIx.slice(plag) = arma::mat(bigK,bigK, arma::fill::eye);
  for(int ihor=plag; ihor < plag+nhor; ihor++){
    arma::mat acc(bigK,bigK, arma::fill::zeros);
    for(int pp=0; pp<plag; pp++){
      arma::mat& Fmatslice = Fmat.slice(pp);
      arma::mat& PHIxslice = PHIx.slice(ihor-pp);
      acc = acc + Fmatslice*PHIxslice;
    }
    PHIx.slice(ihor+1) = acc;
  }
  PHI = PHIx.slices(plag,plag+nhor);
}

void get_nullspace(arma::mat& NU, arma::mat& M){
  arma::mat Q, R; qr(Q, R, M);
  arma::uword r = arma::rank(M);
  if(r == 0){
    arma::uvec set(M.n_cols); for(size_t c = 0; c < M.n_cols; c++){set(c) = c;};
    NU = Q.cols(set);
  }else{
    arma::uvec set(r); for(size_t c = 0; c < r; c++){set(c) = c;};
    NU = Q;
    NU.shed_cols(set);
  }
}

void write_irf_to_bigmat(Rcpp::XPtr<BigMatrix> xpMat, const arma::cube& irf_result,
                         const int irep, const int bigK, const int nhor) {
  // Matrix structure: nrow = bigK*bigK*nhor, ncol = thindraws
  // Each draw (irep) goes to column irep
  MatrixAccessor<double> macc(*xpMat);
  
  // Rcpp::Rcout << "BigMatrix info:" << std::endl;
  // Rcpp::Rcout << "  Type: " << xpMat->matrix_type() << std::endl;
  // Rcpp::Rcout << "  Rows: " << xpMat->nrow() << std::endl;
  // Rcpp::Rcout << "  Cols: " << xpMat->ncol() << std::endl;
  // Rcpp::Rcout << "  Current irep: " << irep << std::endl;
  
  // // Test writing to first few elements
  // for(int i = 0; i < 10; i++) {
  //   mat[i][irep] = static_cast<double>(i);
  // }

  int row_idx = 0;
  for(int h = 0; h < nhor; h++) {
    for(int j = 0; j < bigK; j++) {
      for(int i = 0; i < bigK; i++) {
        macc[irep][row_idx] = irf_result(i, j, h);
        row_idx++;
      }
    }
  }
}

void write_rot_to_bigmat(Rcpp::XPtr<BigMatrix> xpMat, const arma::mat& rot_matrix,
                         const int irep, const int bigK) {
  // Matrix structure: nrow = bigK*bigK, ncol = thindraws
  // Each draw (irep) goes to column irep
  
  MatrixAccessor<double> mat(*xpMat);
  
  int row_idx = 0;
  for(int j = 0; j < bigK; j++) {
    for(int i = 0; i < bigK; i++) {
      mat[row_idx][irep] = rot_matrix(i, j);
      row_idx++;
    }
  }
}

//' @name compute_irf
//' @noRd
// [[Rcpp::export]]
Rcpp::List compute_irf(arma::cube A_large, arma::cube S_large, arma::cube Ginv_large, const int type, const int nhor, const int thindraws, const SEXP shocklist_in, SEXP irf_bigmat, SEXP rot_bigmat, const bool verbose, const bool save_rot, const int threads) {  // input stuff
  Rcpp::List shocklist(Rcpp::clone(shocklist_in));
  //----------------------------------------------------------------------------
  // miscellaneous
  Rcpp::List shock_idx = shocklist["shock.idx"];
  Rcpp::LogicalVector shock_cidx = shocklist["shock.cidx"];
  // parameter
  const int plag = shocklist["plag"];
  const int bigK = A_large.n_rows;
  const int k = A_large.n_cols;
  const int N = shock_idx.size();
  
  // Import Rs chol function
  Rcpp::Environment base = Rcpp::Environment("package:base");
  // Rcpp::Function Rchol = base["chol"];
  
  // allocate the output matrix
  // arma::field<arma::cube> irf_output(thindraws);
  // irf_output.fill(arma::ones<arma::cube>(bigK,bigK,nhor));
  // arma::uword size_of_rot = save_rot ? thindraws : 0;
  // 
  // arma::field<arma::mat> rot_output(size_of_rot);
  // if(save_rot) {
  //   rot_output.fill(arma::ones<arma::mat>(bigK,bigK));
  // }
  // Setup BigMatrix pointers
  Rcpp::XPtr<BigMatrix> xpIrfMat(irf_bigmat);
  Rcpp::XPtr<BigMatrix> xpRotMat(rot_bigmat);

  // Verify dimensions
  const int expected_irf_rows = bigK * bigK * nhor;
  if(xpIrfMat->nrow() != expected_irf_rows || xpIrfMat->ncol() != thindraws) {
    throw std::runtime_error("IRF BigMatrix dimensions incorrect");
  }
  // Verify rotation matrix dimensions
  const int expected_rot_rows = bigK * bigK;
  if(xpRotMat->nrow() != expected_rot_rows || xpRotMat->ncol() != thindraws) {
    throw std::runtime_error("Rotation BigMatrix dimensions incorrect");
  }
  
  arma::vec counter_output(thindraws);
  
  // OpenMP Setup
  // int max_threads = omp_get_num_procs();
  // int use_threads = std::max(1, max_threads - 2);
    int use_threads = threads;
    omp_set_num_threads(use_threads);
    if (verbose) {
      Rcpp::Rcout << "Using " << use_threads << " threads for parallel computation\n";
    }
  
  // Main parallel computation
#pragma omp parallel \
  shared(A_large, S_large, Ginv_large, xpIrfMat, xpRotMat, counter_output, shock_idx, shock_cidx) \
    firstprivate(bigK, plag, nhor, N, type, save_rot)
    {
      // --- preallocate once per thread ---
      arma::cube irfa(bigK, bigK, nhor, arma::fill::zeros);
      arma::cube PHI(bigK, bigK, nhor, arma::fill::zeros);
      arma::cube Fmat(bigK, bigK, plag, arma::fill::zeros);
      arma::mat P0G(bigK, bigK, arma::fill::zeros);
      arma::mat Q_bar(bigK, bigK, arma::fill::eye);
      arma::mat invGSigma_u(bigK, bigK, arma::fill::zeros);
      
      arma::mat Amat(bigK, A_large.n_cols);
      arma::mat Smat(bigK, bigK);
      arma::mat Ginv(bigK, bigK);
      
      // Thread-local RNG
      std::mt19937 gen(std::random_device{}() + omp_get_thread_num());
      std::normal_distribution<double> rnorm(0.0, 1.0);
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      
      Rcpp::Rcout << "Allocation done, calculating";
      
      #pragma omp for schedule(dynamic)
      for(int irep = 0; irep < thindraws; irep++) {
        // reset memory
        irfa.zeros();
        PHI.zeros();
        Fmat.zeros();
        P0G.zeros();
        Q_bar.eye();
        invGSigma_u.zeros();
        
        // assign slices (no new allocation)
        const arma::mat& Amat = A_large.slice(irep);
        const arma::mat& Smat = S_large.slice(irep);
        const arma::mat& Ginv = Ginv_large.slice(irep);
        
    
    // construct Fmat
    for(int pp = 0; pp < plag; pp++){
      Fmat.slice(pp) = Amat.cols(pp*bigK,(pp+1)*bigK-1);
    }
    // create P0G
    bool chol_success = TRUE;
    for(int cc = 0; cc < N; cc++){
      arma::uvec idx = shock_idx[cc];
      if(shock_cidx[cc]){
        // fall back on pivoting if Cholesky do not work
        arma::mat Sig_chol; arma::mat sigma = Smat.submat(idx,idx);
        chol_success = chol(Sig_chol, sigma, "lower");
        if(!chol_success){
          break;
        }
        P0G.submat(idx,idx) = Sig_chol;
      }else{
        P0G.submat(idx,idx) = Smat.submat(idx,idx);
      }
    }
    // iterate through if cholesky fails
    if(!chol_success){
      counter_output(irep) = shocklist["MaxTries"]; // set to MaxTries then draw is ignored afterwards
      continue;
    }
  
    // create PHI
    get_PHI(PHI, Fmat, nhor);
    // find rotation matrix
    Q_bar = arma::eye(bigK,bigK);
  
    if(type==3){
      int MaxTries = shocklist["MaxTries"];
      arma::mat S_cube = shocklist["S.cube"];
      arma::mat P_cube = shocklist["P.cube"];
      arma::cube Z_cube = shocklist["Z.cube"];
      arma::uvec shock_order = shocklist["shock.order"];
      arma::uvec shock_horz = shocklist["shock.horz"];
      Rcpp::LogicalVector nozero = shocklist["no.zero.restr"];
      int H_restr = shock_horz.n_elem;
      int N_restr = bigK * H_restr;
      // build irf_restr
      arma::mat irf_restr(0, bigK);
      for(int hh = 0; hh < H_restr; hh++){
        int horz = shock_horz(hh);
        arma::mat irf_hh = PHI.slice(horz) * Ginv * P0G;
        irf_restr = join_cols(irf_restr, irf_hh);
      }
      
      // sort Zcube
      arma::cube Z_cube_sorted(size(Z_cube), arma::fill::zeros);
      for(int i = 0; i < bigK; i++){
        int idx = shock_order(i);
        Z_cube_sorted.slice(i) = Z_cube.slice(idx);
      }
      // set up while loop
      int icounter = 0;
      double condall = 0.0;
      // initialize stuff before the loop
      arma::vec signCheck(bigK,1);
      while( condall == 0.0 && icounter < MaxTries){
        arma::mat Q(bigK, bigK, arma::fill::eye);
        for(int cc = 0; cc < N; cc++){
          arma::uvec idx = shock_idx[cc];
          int Kidx = idx.size();
          arma::mat randMat(Kidx, Kidx), Qc(Kidx, Kidx, arma::fill::zeros);
          if(shock_cidx[cc]){
            for(int kk=0; kk<Kidx; kk++){
              for(int jj=0; jj<Kidx; jj++){
                randMat(jj,kk) = rnorm(gen);
              }
            }
            if(nozero[cc]){
              arma::mat R; qr(Qc, R, randMat);
            }else{
              for(int i = 0; i < Kidx; i++){
                arma::mat Ztemp = Z_cube_sorted.slice(idx[i]);
                arma::colvec Zsum = sum(abs(Ztemp),1); arma::uvec zidx = find(Zsum);
                Ztemp = Ztemp.rows(zidx);
                arma::mat R(0,Kidx);
                if(i == 0){
                  arma::mat R1 = Ztemp * irf_restr.cols(idx);
                  R = join_cols(R,R1);
                }else{
                  arma::mat R2 = Qc.head_cols(i).t();
                  if(Ztemp.n_elem != 0){
                    arma::mat R1 = Ztemp * irf_restr.cols(idx);
                    R = join_cols(R, R1);
                  }
                  R = join_cols(R, R2);
                }
                arma::mat Rt = R.t();
                arma::mat NU; get_nullspace(NU, Rt); 
                arma::vec x_j = randMat.col(i);
                double div = arma::as_scalar((NU.t() * x_j).t() * (NU.t() * x_j));
                arma::vec q_j = NU * ( NU.t() * x_j / sqrt(div));
                Qc.col(i) = q_j;
              } // end inner for-loop
            } // end inner if-cond
            Q.submat(idx,idx) = Qc;
          } // end outer if-cond
        } // end-for countries
        // reorder again
        for(int i = 0; i < bigK; i++){
          int idx = shock_order(i);
          Q_bar.col(idx) = Q.col(i);
        }
        // create irfcheck
        arma::mat irf_check = irf_restr * Q_bar;
        for(int kk = 0; kk < bigK; kk++){
          arma::vec STemp = S_cube.col(kk);
          arma::vec PTemp = P_cube.col(kk);
          if(sum(abs(STemp))>0){
            arma::vec prob(N_restr);
            for(int nn = 0; nn < N_restr; nn++){
              if(PTemp(nn) > runif(gen)) prob(nn) = 1; else prob(nn) = 0;
            }
            arma::mat PDiag = diagmat(prob);
            arma::vec IrfTemp = sign(irf_check.col(kk));
            double getsum = arma::as_scalar(IrfTemp.t() * PDiag * STemp);
            double chksum = arma::as_scalar(sum(abs(PDiag * STemp)));
            if(getsum == chksum) signCheck(kk) = 1; else signCheck(kk) = 0;
          }else{
            signCheck(kk) = 1;
          }
        }
        condall = prod(signCheck);
        icounter += 1;
      }
      // save counter
      counter_output(irep) = icounter;
    }else{
     counter_output(irep) = 1;
    } // end if-cond type==3
   
   // compute shock
   invGSigma_u = Ginv * P0G * Q_bar;
   
   // compute impulse responses
   for(int ihor = 0; ihor < nhor; ihor++){
   irfa.slice(ihor) = PHI.slice(ihor) * invGSigma_u;
   }
   
   // irf_output(irep) = irfa;
   // if(save_rot) {
   //   rot_output(irep) = Q_bar;
   // }
  // Write results to BigMatrix (thread-safe)
  // write_irf_to_bigmat(xpIrfMat, irfa, irep, bigK, nhor);
  // if(save_rot) {
  //   write_rot_to_bigmat(xpRotMat, Q_bar, irep, bigK);
  // }
  //if(verbose){
      // Increment progress bar
   // if (any(prog_rep_points == irep)){
        //prog.increment();
      //}
   //}
   // check user interruption
   if(irep % 100 == 0) {
   Rcpp::checkUserInterrupt();
   }
  } // end for-loop of impulse response computation
   //----------------------------------------------------------------------------
    }   
   // return Rcpp::List::create(
   // Rcpp::Named("irf", irf_output),
   // Rcpp::Named("rot", rot_output),
   // Rcpp::Named("counter", counter_output)
   // );
   
  return Rcpp::List::create(
    Rcpp::Named("counter", counter_output),
    Rcpp::Named("message", "Results written to BigMatrix objects")
  );
}