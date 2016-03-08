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

mat xDx (mat x, vec d) {
  int p = x.n_cols;
  int n = x.n_rows;
  mat G = zeros<mat>(n,n);
  double diff_ij;

  for (int i=0; i<n; i++) {
    for (int j=0; j<=i; j++) {
      for (int k=0; k<p; k++) {
        diff_ij = x(i,k) - x(j,k);
        G(i,j) = G(i,j) + diff_ij*diff_ij * d[k];
      }
      G(i,j) = sqrt( G(i,j) );
      if (i != j) G(j,i) = G(i,j);
    }
  }

  return G;
}

double log_like_plus_log_prior(vec y, mat X, vec param, mat I, vec priors) {
  double ls2 = param[0];
  double w = param[1];
  double ltau = param[2];
  vec d_vec = param.tail( param.size()-3 );
  //mat D_mat = xDx(X,d_vec);
  mat D_mat = xDx(X,d_vec % d_vec); //check this.

  double a_s2 = priors[0];//2;
  double b_s2 = priors[1];//5;
  double a_phi = priors[2];//0;
  double b_phi = priors[3];//5;
  double a_tau = priors[4];//2;
  double b_tau = priors[5];//5;
  double a_d = priors[6]; //1;
  double b_d = priors[7]; //1;

  double s2 = exp(ls2);
  double phi = (b_phi*exp(w) + a_phi) / (exp(w)+1); // inverse logit
  double tau = exp(ltau);

  mat K = tau * exp(-phi*D_mat);

  double ldet_K,sign;
  mat s2I_plus_K = s2*I + K;

  log_det(ldet_K, sign, s2I_plus_K);
  
  double sum_log_dj_prior = sum( (-b_d-1) * log(1+abs(d_vec) / (a_d * b_d)) );

  double log_prior = ( w-2*log(exp(w)+1) ) - a_s2*ls2 - b_s2/s2 - a_tau*ltau - b_tau/tau + sum_log_dj_prior;
  double log_like = (-.5 * ldet_K - .5 * y.t() * s2I_plus_K.i() * y)[0];

  return log_like + log_prior;
}

//[[Rcpp::export]]
List gp_gdp(vec y, mat X, mat cand_S, vec init, vec priors, int B, int burn, bool printProg) {
  int n = y.size();
  int num_params = cand_S.n_rows;
  mat In = eye<mat>(n,n);
  int acc_rate = 0;
  mat param = zeros<mat>(B+burn,num_params);
  double log_ratio;
  vec cand = zeros<vec>(num_params);
  vec curr = zeros<vec>(num_params);
  List ret;
  clock_t start_time = clock();
  int freq = 50;
  param.row(0) = reshape(init,1,num_params);

  Rcout << endl;
  for (int b=1; b<B+burn; b++) {
    // Update s2, phi, tau:
    curr = vectorise(param.row(b-1));
    cand = mvrnorm(curr, cand_S); // s2, phi, tau, d1,...,dp

    log_ratio = log_like_plus_log_prior(y,X,cand,In,priors) - 
                log_like_plus_log_prior(y,X,curr,In,priors);

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
  ret["X"] = X;
  ret["cand_S"] = cand_S;

  return ret;
}
