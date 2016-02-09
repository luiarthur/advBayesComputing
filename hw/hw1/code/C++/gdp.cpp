//[[Rcpp::depends(RcppArmadillo)]]
#include <func.cpp> // wsample, ldnorm, mvrnorm


//[[Rcpp::export]]
List gdp(vec y, mat x, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  vec lam = mat<vec>(B,J);
  vec tau = mat<vec>(B,J);
  double sig2 = 1;
  vec alpha = ones<vec>(B);
  vec eta = ones<vec>(B);
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  mat D = eye(J,J);
  mat S = eye(J,J);
  double m, sc;
  List ret;

  for (int b=1; b<B; b++) {
    // Update beta
    S = (xx + D.i()).i();
    beta.row(b) = reshape(mvrnorm(S*xy, sig2*S), 1, J);

    // Update lambda_j, tau_j
    for (int j=0; j<J; j++) {
      // Update lambda
      sc = 1 / ( abs(beta(b,j)) / sqrt(sig2) + eta[b-1] );
      lam(b,j) = rgamma(1,alpha[b-1]+1, sc); // shape, scale

      // Update tau
      m = sqrt(sig2) * lam(b,j) / abs(beta(b,j));
      tau(b,j) = 1 /  rInvGaussian(m, lam(b,j) * lam(b,j));
    }

    // Update alpha
    // Update eta


    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  ret["beta"] = beta.tail_rows(B-burn);
  ret["lambda"] = beta.tail_rows(B-burn);
  ret["tau"] = beta.tail_rows(B-burn);
  ret["sig2"] = beta.tail_rows(B-burn);
  ret["alpha"] = beta.tail_rows(B-burn);
  ret["eta"] = beta.tail_rows(B-burn);

  return ret;
}
