//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include <RcppArmadilloExtensions/sample.h> // For wsample

using namespace std;
using namespace arma;
using namespace Rcpp;

double wsample( vec x, vec prob ) {
  NumericVector ret = Rcpp::RcppArmadillo::sample(
      as<NumericVector>(wrap(x)), 1, true, 
      as<NumericVector>(wrap(prob))) ;
  return ret[0] ;
}

double lp (double b, double p, double d, double c) {
  double out;
  if ( p > randu() ) {
    out = -b*b / (2*d);
  } else {
    out = -b*b / (2*c);
  }
  return out;
}

double ll (vec y, mat x, vec beta) {
  vec m = (y - x * beta);
  double out = -sum(m.t() * m) / 2;
  return out;
}

//[[Rcpp::export]]
List spikeAndSlab(vec y, mat x, double p, double d, double c, vec cs, int B, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = zeros<mat>(B,J);
  vec acc = zeros<vec>(J);
  double cand, curr;
  double lg_cand, lg_curr;
  vec bvec_curr, bvec_cand;
  List ret;


  for (int b=1; b<B; b++) {
    for (int j=0; j<J; j++) {
      beta(b,j) = beta(b-1,j); 

      curr = beta(b,j);
      cand = randn() * cs[j] + curr;
      bvec_curr = vectorise( beta.row(b) );
      bvec_cand = bvec_curr;
      bvec_cand(j) = cand;

      lg_curr = ll(y, x, bvec_curr) + lp(curr, p, d, c);
      lg_cand = ll(y, x, bvec_cand) + lp(cand, p, d, c);

      if (lg_cand - lg_curr > log(randu())) {
        beta(b,j) = cand;
        acc[j] = acc[j] + 1;
      }
    }
    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  ret["beta"] = beta;
  ret["acc"] = acc;

  return ret;
}
