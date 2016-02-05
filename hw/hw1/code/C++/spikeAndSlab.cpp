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

double ldnorm (double x, double m, double s2) {
  return -.5*log(s2) -.5 * pow(x-m, 2.0) / s2;
}

vec mvrnorm(vec m, mat S) {
  int n = m.size();
  vec e = randn(n);
  return m + chol(S).t() * e;
}

//Working!

//[[Rcpp::export]]
List spikeAndSlab(vec y, mat x, vec tau2, double g, vec w, int B, int burn, bool printProg) {
  int n = y.size();
  int J = x.n_cols;
  mat beta = ones<mat>(B,J);
  double p1, p0, lp1, lp0;
  mat gam = zeros<mat>(B,J);
  mat xx = x.t() * x;
  mat xy = x.t() * y;
  mat D2 = eye<mat>(J,J);
  int dummy;
  List ret;
  vec m;
  mat S;


  for (int b=1; b<B; b++) {

    // Update beta
    S = (xx + D2.i()).i();
    m = S * xy;
    beta.row(b) = reshape(mvrnorm(m,S),1,J);

    // Update gamma
    for (int j=0; j<J; j++) {
      lp1 = ldnorm(beta(b,j),0,g*tau2[j]) + log(w[j]);
      lp0 = ldnorm(beta(b,j),0,tau2[j]) + log(1-w[j]);

      p0 = exp( lp0 - max({lp0,lp1}) );
      p1 = exp( lp1 - max({lp0,lp1}) );

      gam(b,j) = wsample( {0,1}, {p0,p1} );
      if ( gam(b,j) == 1 ) {
        D2(j,j) = g * tau2[j];
      } else {
        D2(j,j) = tau2[j];
      }
    }

    if (printProg) Rcout << "\rProgress: " << b << "/" << B;
  }

  ret["beta"] = beta.rows(burn,B-1);
  ret["gam"] = gam.rows(burn,B-1);

  return ret;
}
