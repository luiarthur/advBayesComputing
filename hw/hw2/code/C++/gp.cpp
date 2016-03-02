//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm
#include <time.h> // timing

void time_remain(clock_t start_time, int it, int total_its, int freq) {
  clock_t end_time;
  end_time = clock();
  int its_remaining = total_its - it;
  if (it % freq == 0) {
    double time_elapsed = ((float)end_time - (float)start_time) / CLOCKS_PER_SEC * its_remaining / freq ;
    int mins = (int) time_elapsed / 60;
    int secs = (int) time_elapsed % 60;
    Rcout << "\rProgress: " << it << "/" << total_its+1 <<"; Time Remaining: " << mins << "m " << secs << "s       ";
  }
}

mat wood_inv(double s2, mat I, mat C, mat Ct, mat Ks) {
  mat out = I/s2 - C/s2* (Ks + Ct*C / s2).i() * Ct/s2;
  return out;
}

double wood_ldet(mat C, mat Ks, mat Ct, double s2, int n) {
  double ldet_1,ldet_2, ldet_3, sign;

  log_det(ldet_1,sign, Ks + Ct*C/s2);
  log_det(ldet_2,sign, Ks);
  ldet_3 = n*log(s2);

  return ldet_1 - ldet_2 + ldet_3;
}


double log_like_plus_log_prior(vec y, vec param, mat D, mat Cd, mat I, vec priors) {
  double ls2 = param[0];
  double w = param[1];
  double ltau = param[2];

  double a_s2 = priors[0];//2;
  double b_s2 = priors[1];//5;
  double a_phi = priors[2];//0;
  double b_phi = priors[3];//5;
  double a_tau = priors[4];//2;
  double b_tau = priors[5];//5;

  double s2 = exp(ls2);
  double phi = (b_phi*exp(w)+a_phi) / (exp(w)+1);
  double tau = exp(ltau);

  mat Ks = tau * exp(-phi*D);

  mat C = tau * exp(-phi*Cd);
  mat Ct = C.t();

  int n = I.n_rows;
  double ldet;
  
  ldet = wood_ldet(C,Ks,Ct,s2,n);

  double log_prior = ( w-2*log(exp(w)+1) ) - a_s2*ls2 - b_s2*exp(-ls2) - a_tau*ltau - b_tau*exp(-ltau);
  double log_like = (-.5 * ldet - .5 * y.t() * wood_inv(s2,I,C,Ct,Ks) * y)[0];

  return log_like + log_prior;
}

//[[Rcpp::export]]
List gp(vec y, mat x, mat s, mat Cd, mat D, mat cand_S, vec init, vec priors, int B, int burn, bool printProg) {
  int n = y.size();
  int num_params = cand_S.n_rows;
  mat In = eye<mat>(n,n);
  int acc_rate = 0;
  mat param = zeros<mat>(B+burn,num_params);
  double log_ratio;
  vec cand,curr;
  List ret;
  clock_t start_time = clock();
  int freq = 50;
  param.row(0) = reshape(init,1,num_params);

  for (int b=1; b<B+burn; b++) {
    curr = vectorise(param.row(b-1));
    cand = mvrnorm(curr, cand_S); // s2, phi, tau

    log_ratio = log_like_plus_log_prior(y,cand,D,Cd,In,priors) - 
                log_like_plus_log_prior(y,curr,D,Cd,In,priors);

    if ( log_ratio > log(randu()) ) {
      param.row(b) = reshape(cand,1,num_params);
      if (b > burn) acc_rate++;
    } else {
      param.row(b) = param.row(b-1);
    }

    if (printProg) time_remain(start_time, b, B+burn-1, freq);
    if (b % freq == 0) start_time = clock();
  }
  Rcout << endl;

  param.col(0) = exp(param.col(0));
  param.col(1) = (priors[3]*exp(param.col(1))+priors[2]) / ( exp(param.col(1))+1 );
  param.col(2) = exp(param.col(2));

  ret["param"] = param.tail_rows(B);
  ret["acc_rate"] = acc_rate * 1.0 / B;
  ret["y"] = y;
  ret["x"] = x;
  ret["s"] = s;
  ret["Cd"] = Cd;
  ret["D"] = D;
  ret["cand_S"] = cand_S;

  return ret;
}
