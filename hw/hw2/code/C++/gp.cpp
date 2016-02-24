//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm



mat wood_inv(double s2, mat I, mat H, mat Ht, mat Ksi) {
  mat out = I/s2 - H/s2* (Ksi + Ht*H / s2).i() * Ht/s2;
  return out;
}

double log_dinvgamma(double x, double a, double b) {
  return (-a-1) * log(x) - b/x;
}

double log_like(vec y, double s2, double phi, double tau, mat D, mat C, mat I) {
  mat Ks = tau * exp(-phi*D);
  mat Ksi = Ks.i();
  mat H = C * Ksi;
  mat Ht = H.t();
  int n = I.n_rows;
  double val, sign, ldet;

  log_det(val,sign, I/s2 + H*Ks*Ht);
  ldet = val;

  return (-.5 * ldet - .5 * y.t() * wood_inv(s2,I,H,Ht,Ksi) * y)[0];
}


//[[Rcpp::export]]
List gp(vec y, mat x, mat s, mat C, mat D, double cs_tau, double cs_phi, double cs_sig2, int B, int burn, bool printProg) {
  int n = y.size();
  mat In = eye<mat>(n,n);
  int acc_tau = 0;
  int acc_phi = 0;
  int acc_sig2 = 0;
  vec tau = ones<vec>(B+burn);
  vec phi = ones<vec>(B+burn);
  vec sig2 = ones<vec>(B+burn);
  double log_ratio;
  double cand;
  List ret;

  for (int b=1; b<B+burn; b++) {
    tau[b] = tau[b-1];
    phi[b] = phi[b-1];
    sig2[b] = sig2[b-1];

    // Update tau: IG(2,5) prior
    cand = randn() * cs_tau + tau[b];
    if (cand > 0) {
      log_ratio = log_like(y, sig2[b], phi[b], cand, D, C, In) + 
                  log_dinvgamma(cand,2,5) - 
                  log_like(y, sig2[b], phi[b], tau[b], D, C, In) -
                  log_dinvgamma(tau[b],2,5);
      if ( log_ratio > log(randu()) ) {
        tau[b] = cand;
        acc_tau++;
      }
    }

    // Update phi: Uniform(0,5) prior
    cand = randn() * cs_phi + phi[b];
    if (0 < cand && cand < 5) {
      log_ratio = log_like(y, sig2[b] , cand, tau[b], D, C, In) - 
                  log_like(y, sig2[b], phi[b], tau[b], D, C, In);
      if ( log_ratio > log(randu()) ) {
        phi[b] = cand;
        acc_phi++;
      }
    }


    // Update sig2: InvGamma(2,1) prior
    cand = randn() * cs_sig2 + sig2[b];
    if (cand > 0) {
      log_ratio = log_like(y, cand, phi[b], tau[b], D, C, In) + 
                  log_dinvgamma(cand,2,1) - 
                  log_like(y, sig2[b], phi[b], tau[b], D, C, In) -
                  log_dinvgamma(sig2[b],2,1);
      if ( log_ratio > log(randu()) ) {
        sig2[b] = cand;
        acc_sig2++;
      }
    }


    if (printProg) Rcout << "\rProgress: " << b+1 << "/" << B+burn;
  }

  vec out_tau  = tau.tail(B);   ret["tau"]  = out_tau;
  vec out_phi  = phi.tail(B);   ret["phi"]  = out_phi;
  vec out_sig2 = sig2.tail(B);  ret["sig2"] = out_sig2;

  return ret;
}
