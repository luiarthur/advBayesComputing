set.seed(268)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../../../../cpp_functions/func.cpp")
sourceCpp("../C++/gp.cpp")
system("mkdir -p output")

n <- 1000
sig2 <- .5
sn <- 100
s <- matrix(rnorm(sn*3),sn,3) # knots
x <- matrix(rnorm(n*3),n,3)     # data (simulated covariates)
y <- x[,1] + ifelse(x[,2]-.5 > 0, x[,2]-.5, 0) + x[,3]^2 + rnorm(n,0,sqrt(sig2)) # data (simulated responses)
C <- cov(t(rbind(x,s)))[1:n,-c(1:n)] # covariance between data and knots
D <- as.matrix(dist(s))

# y | ... ~ N(0,s^2 + K)
sourceCpp("../C++/gp.cpp")
system.time( out <- gp(y, x, s, C, D, cs_tau=1, cs_phi=1, cs_sig2=1, B=3000, burn=1000, printProg=T) )

par(mfrow=c(3,1))
  plot(out$tau,type="l")
  plot(out$phi,type="l")
  plot(out$sig2,type="l")
par(mfrow=c(1,1))
