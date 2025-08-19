#include <RcppArmadillo.h>
#include <random>
#include <Rmath.h>
#include <R_ext/Random.h>
// #include <bigmemory/BigMatrix.h>
// #include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

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

// void write_irf_to_bigmat(Rcpp::XPtr<BigMatrix> xpMat, const arma::cube& irf_result,
//                          const int draw_idx, const int bigK, const int nhor) {
//   // Matrix structure: nrow = bigK*bigK*nhor, ncol = total_draws
//   // Each draw goes to column draw_idx
//   MatrixAccessor<double> macc(*xpMat);
//   
//   int row_idx = 0;
//   for(int h = 0; h < nhor; h++) {
//     for(int j = 0; j < bigK; j++) {
//       for(int i = 0; i < bigK; i++) {
//         macc[draw_idx][row_idx] = irf_result(i, j, h);
//         row_idx++;
//       }
//     }
//   }
// }
// 
// void write_rot_to_bigmat(Rcpp::XPtr<BigMatrix> xpMat, const arma::mat& rot_matrix,
//                          const int draw_idx, const int bigK) {
//   // Matrix structure: nrow = bigK*bigK, ncol = total_draws
//   // Each draw goes to column draw_idx
//   
//   MatrixAccessor<double> mat(*xpMat);
//   
//   int row_idx = 0;
//   for(int j = 0; j < bigK; j++) {
//     for(int i = 0; i < bigK; i++) {
//       mat[row_idx][draw_idx] = rot_matrix(i, j);
//       row_idx++;
//     }
//   }
// }

//' @name compute_irf
 //' @noRd
 // [[Rcpp::export]]
 Rcpp::List compute_irf(arma::mat A_draw, arma::mat S_draw, arma::mat Ginv_draw, 
                                    const int type, const int nhor, const int draw_idx,
                                    const SEXP shocklist_in, 
                                    const bool save_rot, const int seed = -1) {
   
   Rcpp::List shocklist(Rcpp::clone(shocklist_in));
   
   //----------------------------------------------------------------------------
   // miscellaneous
   Rcpp::List shock_idx = shocklist["shock.idx"];
   Rcpp::LogicalVector shock_cidx = shocklist["shock.cidx"];
   // parameter
   const int plag = shocklist["plag"];
   const int bigK = A_draw.n_rows;
   const int k = A_draw.n_cols;
   const int N = shock_idx.size();
   
   // Setup BigMatrix pointers
   // Rcpp::XPtr<BigMatrix> xpIrfMat(irf_bigmat);
   // Rcpp::XPtr<BigMatrix> xpRotMat(rot_bigmat);
   
   // Verify dimensions
   // const int expected_irf_rows = bigK * bigK * nhor;
   // if(xpIrfMat->nrow() != expected_irf_rows) {
   //   throw std::runtime_error("IRF BigMatrix row dimensions incorrect");
   // }
   // 
   // if(save_rot) {
   //   const int expected_rot_rows = bigK * bigK;
   //   if(xpRotMat->nrow() != expected_rot_rows) {
   //     throw std::runtime_error("Rotation BigMatrix row dimensions incorrect");
   //   }
   // }
   
   // Check column bounds
   // if(draw_idx >= xpIrfMat->ncol()) {
   //   throw std::runtime_error("draw_idx exceeds IRF BigMatrix column dimensions");
   // }
   // if(save_rot && draw_idx >= xpRotMat->ncol()) {
   //   throw std::runtime_error("draw_idx exceeds Rotation BigMatrix column dimensions");
   // }
   
   // Initialize RNG if seed provided
   std::mt19937 gen;
   if(seed > 0) {
     gen.seed(seed);
   } else {
     gen.seed(std::random_device{}());
   }
   std::normal_distribution<double> rnorm(0.0, 1.0);
   std::uniform_real_distribution<double> runif(0.0, 1.0);
   
   // --- Preallocate matrices ---
   arma::cube irfa(bigK, bigK, nhor, arma::fill::zeros);
   arma::cube PHI(bigK, bigK, nhor, arma::fill::zeros);
   arma::cube Fmat(bigK, bigK, plag, arma::fill::zeros);
   arma::mat P0G(bigK, bigK, arma::fill::zeros);
   arma::mat Q_bar(bigK, bigK, arma::fill::eye);
   arma::mat invGSigma_u(bigK, bigK, arma::fill::zeros);
   
   // Assign matrices
   const arma::mat& Amat = A_draw;
   const arma::mat& Smat = S_draw;
   const arma::mat& Ginv = Ginv_draw;
   
   // construct Fmat
   for(int pp = 0; pp < plag; pp++){
     Fmat.slice(pp) = Amat.cols(pp*bigK,(pp+1)*bigK-1);
   }
   
   // create P0G
   bool chol_success = true;
   for(int cc = 0; cc < N; cc++){
     arma::uvec idx = shock_idx[cc];
     if(shock_cidx[cc]){
       // fall back on pivoting if Cholesky do not work
       arma::mat Sig_chol; 
       arma::mat sigma = Smat.submat(idx,idx);
       chol_success = chol(Sig_chol, sigma, "lower");
       if(!chol_success){
         break;
       }
       P0G.submat(idx,idx) = Sig_chol;
     }else{
       P0G.submat(idx,idx) = Smat.submat(idx,idx);
     }
   }
   
   // Check if cholesky failed
   int counter_result = 1;
   if(!chol_success){
     counter_result = shocklist["MaxTries"]; // set to MaxTries then draw is ignored afterwards
     return Rcpp::List::create(
       Rcpp::Named("counter", counter_result),
       Rcpp::Named("success", false),
       Rcpp::Named("message", "Cholesky decomposition failed")
     );
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
     counter_result = icounter;
   } // end if-cond type==3
   
   // compute shock
   invGSigma_u = Ginv * P0G * Q_bar;
   
   // compute impulse responses
   for(int ihor = 0; ihor < nhor; ihor++){
     irfa.slice(ihor) = PHI.slice(ihor) * invGSigma_u;
   }
   
   // Write results to BigMatrix
   // write_irf_to_bigmat(xpIrfMat, irfa, draw_idx, bigK, nhor);
   // if(save_rot) {
   //   write_rot_to_bigmat(xpRotMat, Q_bar, draw_idx, bigK);
   // }
   // 
   // return Rcpp::List::create(
   //   Rcpp::Named("counter", counter_result),
   //   Rcpp::Named("success", true),
   //   Rcpp::Named("message", "Results written to BigMatrix objects")
   // );
   return Rcpp::List::create(
     Rcpp::Named("counter", counter_result),
     Rcpp::Named("success", true),
     Rcpp::Named("irf_result", irfa),
     Rcpp::Named("rot_result", save_rot ? Rcpp::wrap(Q_bar) : R_NilValue)
   );
 }