//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include <../../../02_dpmm/code/Rcpp/dp/dp.h> //wsample

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
List spikeAndSlab(mat Y, mat D, double beta_m, double beta_s2,
    double tau2_a, double tau2_b, double alpha_a, double alpha_b,
    double sig2_a, double sig2_b, double phi_a, double phi_b, int L, int B) {

  return ret;
}
