//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm

//[[Rcpp::export]]

//[[Rcpp::export]]
mat gp(vec y, mat x, mat V, int B, int burn, bool printProg) {
  int n = y.size();
  int p = x.n_cols; 
  mat xt = x.t();
  mat S = (xt * x + V.i()).i();
  mat SX = S*xt;
  mat beta = zeros<mat>(B+burn,p);

  for (int b=1; b<B+burn; b++) {
    beta.row(b) = reshape(mvrnorm(SX * z, S),1,p);

    if (printProg) Rcout << "\rProgress: " << b+1 << "/" << B+burn;
  }
  
  return beta.tail_rows(B);
}
