//[[Rcpp::depends(RcppArmadillo)]]
#include <func.h> // wsample, ldnorm, mvrnorm


//[[Rcpp::export]]
List blasso(vec y, mat x, double r, double delta, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  vec lam = ones<vec>(B);
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  mat outvec = zeros<mat>(B-burn,1);
  mat S = zeros<mat>(J,J);
  mat D = eye(J,J);
  vec m = zeros<vec>(J);
  mat t2 = ones<mat>(B,J);
  double sh = J+r;
  double sc;
  double l2;
  List ret;


  for (int b=1; b<B; b++) {

    // Update t2
    for (int j=0; j<J; j++) {
      l2 = lam[b-1] * lam[b-1];
      t2(b,j) = 1 / rinvGaussian(sqrt(l2/beta(b-1,j)),l2);
      D(j,j) = 1;
    }

    // Update beta
    S = (xx + D.i()).i();
    beta.row(b) = reshape(mvrnorm(S*xy,S),1,J);
    
    // Update lambda
    sc = sum(t2.row(b))/2 + delta;
    lam[b] = rgamma(1,sh,1/sc)[0]; // shape and rate
    
    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  outvec = reshape(lam,B,1);
  ret["beta"] = beta.tail_rows(B-burn);
  ret["lambda"] = outvec.tail_rows(B-burn);
  ret["t2"] = t2.tail_rows(B-burn);

  return ret;
}
