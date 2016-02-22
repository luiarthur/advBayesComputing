set.seed(268)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../../../../cpp_functions/func.cpp")
sourceCpp("../C++/gp.cpp")
system("mkdir -p output")

n <- 5000
sig2 <- .5
x <- matrix(rnorm(n*3),n,3)
y <- x[,1] + ifelse(x[,2]-.5 > 0, x[,2]-.5, 0) + x[,3]^2 + rnorm(n,0,sqrt(sig2))


