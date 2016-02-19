//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm

//[[Rcpp::export]]
vec rtruncnorm(vec a, vec b, vec mu, vec sigma) {
  int n = mu.size();
  vec draw(n);

  for (int i=0; i<n; i++) {
    double lo = Rcpp::stats::pnorm_1((a(i)-mu(i))/sigma(i), 0, true, false);
    double hi = Rcpp::stats::pnorm_1((b(i)-mu(i))/sigma(i), 0, true, false);
    double u = Rf_runif(0, 1);
    double res= mu(i) + sigma(i) * Rcpp::stats::qnorm_1(u*(lo-hi)+hi, 0, true, false);
    draw(i)=res;
  }

  return draw;
}

//[[Rcpp::export]]
vec update_z(vec y, mat x, vec b) {
  int n = y.size();
  vec z = zeros<vec>(n);

  uvec eq_0 = find(y==0);
  uvec eq_1 = find(y==1);
  int j = eq_0.size();
  int k = eq_1.size();
  double inf = 100000;
  
  mat x0 = x.rows(eq_0);
  mat x1 = x.rows(eq_1);

  z( eq_0 ) = rtruncnorm(-inf*ones<vec>(j),zeros<vec>(j),x0*b,ones<vec>(j) );
  z( eq_1 ) = rtruncnorm( zeros<vec>(k),inf*ones<vec>(k),x1*b,ones<vec>(k) );
 
  return z;
}

//[[Rcpp::export]]
mat gibbs_ac(vec y, mat x, mat V, int B, int burn, bool printProg) {
  int n = y.size();
  int p = x.n_cols; 
  mat xt = x.t();
  mat S = (xt * x + V.i()).i();
  mat SX = S*xt;

  vec z = ones<vec>(n);
  mat beta = zeros<mat>(B+burn,p);

  for (int b=1; b<B+burn; b++) {
    z = update_z(y,x,vectorise(beta.row(b-1)));
    beta.row(b) = reshape(mvrnorm(SX * z, S),1,p);

    if (printProg) Rcout << "\rProgress: " << b+1 << "/" << B+burn;
  }
  
  return beta.tail_rows(B);
}
