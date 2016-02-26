//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm
#include <ctime> // timing

void time_remain(clock_t start_time, int it, int total_its, int freq) {
  clock_t end_time = clock();
  int its_remaining = total_its - it;
  if (it % freq == 0) {
    double out = ( (end_time-start_time) / 100000.0 / freq / 60) * its_remaining;
    Rcout << "\rProgress: " << it << "/" << total_its+1 <<"; Time Remaining: " << out << " mins        ";
  }
}

mat wood_inv(double s2, mat I, mat C, mat Ct, mat Ksi) {
  mat out = I/s2 - C/s2* (Ksi + Ct*C / s2).i() * Ct/s2;
  return out;
}

double wood_ldet(mat C, mat Ks, mat Ksi, mat Ct, double s2, int n) {
  double ldet_1,ldet_2, ldet_3, sign;

  log_det(ldet_1,sign, Ksi + Ct*C/s2);
  log_det(ldet_2,sign, Ks);
  ldet_3 = n*log(s2);

  return ldet_1 + ldet_2 + ldet_3;
}

double log_like(vec y, double ls2, double w, double ltau, mat D, mat C, mat I) {
  double s2 = exp(ls2);
  double phi = 5 / (exp(-w)+1);
  double tau = exp(ltau);

  mat Ks = tau * exp(-phi*D);
  mat Ksi = Ks.i();
  mat Ct = C.t();
  int n = I.n_rows;
  double ldet;
  
  ldet = wood_ldet(C,Ks,Ksi,Ct,s2,n);

  return (-.5 * ldet - .5 * y.t() * wood_inv(s2,I,C,Ct,Ksi) * y)[0];
}

double log_prior(vec param) { // s2, phi, tau
  double ls2 = param[0];
  double w = param[1];
  double ltau = param[2];

  double a_s2 = 2;
  double b_s2 = 5;
  double a_phi = 0;
  double b_phi = 5;
  double a_tau = 2;
  double b_tau = 5;

  //return lphi - a_s2*ls2 - b_s2/exp(ls2) - a_tau*ltau - b_tau/exp(ltau);
  return ( w-2*log(exp(w)+1) ) - a_s2*ls2 - b_s2/exp(ls2) - a_tau*ltau - b_tau/exp(ltau);
}

//[[Rcpp::export]]
List gp(vec y, mat x, mat s, mat C, mat D, mat cand_S, int B, int burn, bool printProg) {
  int n = y.size();
  int p = x.n_cols;
  mat In = eye<mat>(n,n);
  int acc_rate = 0;
  mat param = zeros<mat>(B+burn,3);
  double log_ratio;
  vec cand;
  List ret;
  time_t start;
  int freq = 50;

  for (int b=1; b<B+burn; b++) {
    if (b % freq == 0) start = clock();

    // Update tau: IG(2,5) prior
    cand = mvrnorm(vectorise(param.row(b-1)), cand_S); // s2, phi, tau
    log_ratio = log_like(y, cand[0], cand[1], cand[2], D, C, In) + 
                log_prior(cand) - 
                log_like(y, param(b-1,0), param(b-1,1), param(b-1,2), D, C, In) -
                log_prior(  vectorise(param.row(b-1))  );
    if ( log_ratio > log(randu()) ) {
      param.row(b) = reshape(cand,1,p);
      if (b > burn) acc_rate++;
    }

    //if (printProg) Rcout << "\rProgress: " << b << "/" << B+burn-1;
    if (printProg) time_remain(start, b, B+burn-1, freq);
  }

  param.col(0) = exp(param.col(0));
  param.col(1) = (5*exp(param.col(1))+0) / ( exp(param.col(1))+1 );
  param.col(2) = exp(param.col(2));

  ret["param"] = param.tail_rows(B);
  ret["acc_rate"] = acc_rate * 1.0 / B;

  return ret;
}
