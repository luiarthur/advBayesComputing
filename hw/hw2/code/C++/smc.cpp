//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm

//http://www.cs.ubc.ca/~arnaud/delmoral_doucet_jasra_sequentialmontecarloforbayesiancomputation_valenciaalmostfinal.pdf

//[[Rcpp::export]]
mat update_smc(vec y, mat x, vec theta, int s, int B, int burn, bool printProg) {
  int n = y.size();
  int p = x.n_cols; 
  mat xt = x.t();
  mat S = (xt * x + V.i()).i();
  mat SX = S*xt;
  mat beta = zeros<mat>(B+burn,p);
  vec w = zeros<vec>(s);
  vec z = zeros<vec>(p);

  //vec tt = sample_replace( {1,2,3}, 5, {1,1,1} );
  //Rcout << tt << endl;

  for (int j=0; j<s; j++) {
  }


  //for (int b=1; b<B+burn; b++) {
  //  if (printProg) Rcout << "\rProgress: " << b+1 << "/" << B+burn;
  //}
  
  return beta.tail_rows(B);
}
