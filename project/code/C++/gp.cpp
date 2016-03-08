//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm
#include <time.h> // timing

void time_remain(clock_t start_time, int it, int total_its, int freq) {

  if (it % freq == 0 || it==1) {
    if (it == 1) freq = 1;

    clock_t end_time;
    end_time = clock();
    int its_remaining = total_its - it;

    double time_elapsed = ((float)end_time - (float)start_time) / CLOCKS_PER_SEC * its_remaining / freq ;
    int mins = (int) time_elapsed / 60;
    int secs = (int) time_elapsed % 60;
    Rcout << "\rProgress: " << it << "/" << total_its+1 <<"; Time Remaining: " << mins << "m " << secs << "s       ";
  }
}

double log_like_plus_log_prior(vec y, vec param, mat D, mat I, vec priors) {
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
  double phi = (b_phi*exp(w) + a_phi) / (exp(w)+1); // inverse logit
  double tau = exp(ltau);

  mat K = tau * exp(-phi*D);

  double ldet_K,sign;
  mat s2I_plus_K = s2*I + K;

  ldet_K = log_det(ldet_K,sign, s2I_plus_K);
  
  double log_prior = ( w-2*log(exp(w)+1) ) - a_s2*ls2 - b_s2/s2 - a_tau*ltau - b_tau/tau;
  double log_like = (-.5 * ldet_K - .5 * y.t() * s2I_plus_K.i() * y)[0];

  return log_like + log_prior;
}

//[[Rcpp::export]]
List gp(vec y, mat x, mat D, mat cand_S, vec init, vec priors, int B, int burn, bool printProg) {
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

  Rcout << endl;
  for (int b=1; b<B+burn; b++) {
    curr = vectorise(param.row(b-1));
    cand = mvrnorm(curr, cand_S); // s2, phi, tau

    log_ratio = log_like_plus_log_prior(y,cand,D,In,priors) - 
                log_like_plus_log_prior(y,curr,D,In,priors);

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
  param.col(1) = (priors[3]*exp(param.col(1))+priors[2]) / ( exp(param.col(1))+1 );// inverse logit
  param.col(2) = exp(param.col(2));

  Rcout <<"Acceptance Rate: " << acc_rate * 1.0 / B << endl;
  Rcout <<"The parameters in $param are 's2,phi,tau'" << endl;

  ret["param"] = param.tail_rows(B); //s2, phi, tau
  ret["acc_rate"] = acc_rate * 1.0 / B;
  ret["y"] = y;
  ret["x"] = x;
  ret["D"] = D;
  ret["cand_S"] = cand_S;

  return ret;
}
