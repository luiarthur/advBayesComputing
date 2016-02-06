//[[Rcpp::depends(RcppArmadillo)]]
#include <func.cpp> // wsample, ldnorm, mvrnorm


//[[Rcpp::export]]
List gdp(vec y, mat x, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  List ret;

  for (int b=1; b<B; b++) {

    
    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  ret["beta"] = beta.tail_rows(B-burn);

  return ret;
}
