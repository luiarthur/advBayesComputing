//[[Rcpp::depends(RcppArmadillo)]]
#include <../../../cpp_functions/func.cpp> // wsample, ldnorm, mvrnorm
#include <time.h> // timing


//[[Rcpp::export]]
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

