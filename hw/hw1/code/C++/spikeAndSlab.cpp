//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include <RcppArmadilloExtensions/sample.h> // For wsample

using namespace std;
using namespace arma;
using namespace Rcpp;
const double pi  =3.141592653589793238462;

double wsample( vec x, vec prob ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(
      as<NumericVector>(wrap(x)), 1, true, 
      as<NumericVector>(wrap(prob))) ;
  return ret[0] ;
}

double pdfnorm (double x, double m, double s2) {
  return exp(-.5 * pow(x-m,2.0)/s2) / sqrt(2*pi*s2);
}

double lpj (double bj, double gamj, double tau2j, double g) {
  double prior;
  prior = (1-gamj) * pdfnorm(bj,0,tau2j) + gamj * pdfnorm(bj,0,g*tau2j);
  return log(prior);
}

double lhp(double gamj, double wj) {
  return gamj * log(wj) + (1-gamj) * log(1-wj);
}

double ll (vec y, mat x, vec beta) {
  vec m = (y - x * beta);
  double out = -sum(m.t() * m) / 2;
  return out;
}

//[[Rcpp::export]]
List spikeAndSlab(vec y, mat x, vec tau2, double g, vec w, vec cs, int B, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = zeros<mat>(B,J);
  vec beta_acc = zeros<vec>(J);
  vec gam_acc = zeros<vec>(J);
  double cand, curr;
  double lg_cand, lg_curr;
  vec bvec_curr, bvec_cand;
  mat gam = zeros<mat>(B,J);
  List ret;


  for (int b=1; b<B; b++) {
    for (int j=0; j<J; j++) {

      // Update beta
      beta(b,j) = beta(b-1,j); 
      curr = beta(b,j);
      cand = randn() * cs[j] + curr;
      bvec_curr = vectorise( beta.row(b) );
      bvec_cand = bvec_curr;
      bvec_cand(j) = cand;

      lg_curr = ll(y, x, bvec_curr) + lpj(curr, gam(b,j), tau2[j], g);
      lg_cand = ll(y, x, bvec_cand) + lpj(cand, gam(b,j), tau2[j], g);

      if (lg_cand - lg_curr > log(randu())) {
        beta(b,j) = cand;
        beta_acc[j] = beta_acc[j] + 1;
      }

      // Update gamma
      gam(b,j) = gam(b-1,j);
      curr = gam(b,j);
      cand = wsample({0.0,1.0}, {.5,.5});
      lg_curr = lpj(beta(b,j), curr, tau2[j], g) + lhp(curr,w[j]);
      lg_cand = lpj(beta(b,j), cand, tau2[j], g) + lhp(cand,w[j]);

      if (lg_cand - lg_curr > log(randu())) {
        gam(b,j) = cand;
        gam_acc[j] = gam_acc[j] + 1;
      }
    }
    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  ret["beta"] = beta;
  ret["gam"] = gam;
  ret["beta_acc"] = beta_acc;
  ret["gam_acc"] = gam_acc;

  return ret;
}
