//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm

//http://www.cs.ubc.ca/~arnaud/delmoral_doucet_jasra_sequentialmontecarloforbayesiancomputation_valenciaalmostfinal.pdf
//http://luiarthur.github.io/assets/ams268/notes/smc.pdf

//R::pnorm(double x, double m, double s, bool lower.tail,bool log)

double log_like(vec y, mat x, mat beta_row) {
  vec beta = vectorise(beta_row);
  vec xb = x*beta;
  int n = y.size();
  double log_p;
  double log_q;
  vec like = zeros<vec>(n);

  for (int i=0; i<n; i++) {
    log_p = R::pnorm(xb[i], 0, 1, true, true);
    log_q = R::pnorm(xb[i], 0, 1, false, true);
    like[i] = y[i]*log_p+ (1-y[i])*log_q;
  }

  return sum(like);
}


//[[Rcpp::export]]
mat update_smc(vec y, mat x, int N, mat beta, bool printProg) {
  int T = y.size();
  int p = x.n_cols; 
  mat xt = x.t();
  mat beta_new = zeros<mat>(N,p);
  vec w_old = ones<vec>(N) / N;
  vec w_new = ones<vec>(N);
  vec z = zeros<vec>(p);
  vec ind = zeros<vec>(N);

  for (int t=1; t<T; t++) {
    for (int i=0; i<N; i++) {
      w_new[i] = w_old[i] + log_like({y[t+1]},x,beta.row(i));
    }
    if (printProg) Rcout << "\rProgress: " << t+1 << "/" << T;
  }

  w_new = exp(w_new - max(w_new));
  ind = sample_replace(linspace(0,N-1,N), N, w_new );
  for (int i=0; i<N; i++) {
    beta_new.row(i) = beta.row(ind[i]);
  }

  return beta_new;
}
