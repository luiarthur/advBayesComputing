//[[Rcpp::depends(RcppArmadillo)]]
#include <func.h> // wsample, ldnorm, mvrnorm


double lg(vec y, mat x, double l, vec b, vec cand) {
  vec m;
  double lg_cand, lg_curr;
  
  m = y-x*cand;
  lg_cand = (-.5 * m.t() * m -.5 * l * sum(abs(cand)))[0];
  m = y-x*b;
  lg_curr = (-.5 * m.t() * m -.5 * l * sum(abs(b)))[0];

  return lg_cand-lg_curr;

}

//[[Rcpp::export]]
List blasso(vec y, mat x, double r, double delta, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  vec lam = ones<vec>(B);
  int beta_acc = 0;
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  mat outvec = zeros<mat>(B-burn,1);
  vec cand;
  List ret;


  for (int b=1; b<B; b++) {

    // Update beta
    cand = vectorise(mvrnorm(beta.row(b-1),eye(J,J)));
    if ( lg(y,x,lam[b-1],beta.row(b-1),cand) > log(randu()) ) {
      beta_acc++;
      beta.row(b) = reshape(cand, 1, J);
    } else {
      beta.row(b) = beta.row(b-1);
    }
    

    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  outvec = reshape(lam,B,1);
  ret["beta"] = beta.tail_rows(B-burn);
  ret["lambda"] = outvec.tail_rows(B-burn);
  ret["beta_acc"] = beta_acc;

  return ret;
}
