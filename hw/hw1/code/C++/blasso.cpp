//[[Rcpp::depends(RcppArmadillo)]]
#include <func.h> // wsample, ldnorm, mvrnorm


//[[Rcpp::export]]
List blasso(vec y, mat x, double r, double delta, vec t2_a, vec t2_b, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  vec lam = ones<vec>(B);
  int beta_acc = 0;
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  mat outvec = zeros<mat>(B-burn,1);
  mat S = zeros<mat>(J,J);
  mat D = zeros<mat>(J,J);
  vec m = zeros<vec>(J);
  double sh, sc;
  List ret;


  for (int b=1; b<B; b++) {

    // Update beta
    S = (xx + D.i()).i();


    
    // Update lambda
    //sh = J+r;
    //sc = sum(t2)/2 + delta;
    //lam[b] = rgamma(1,sh,1/sc); // shape and rate
    

    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  outvec = reshape(lam,B,1);
  ret["beta"] = beta.tail_rows(B-burn);
  ret["lambda"] = outvec.tail_rows(B-burn);
  ret["beta_acc"] = beta_acc;

  return ret;
}
