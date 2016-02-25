set.seed(268)
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../../R_Functions/plotPost.R",chdir=T)
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
system.time( out <- gp(y, x, s, C, D, cand_S=diag(3), B=2000, burn=2000, printProg=T) )

par(mfrow=c(3,1))
  plot(out$param[,1],type="l",ylab=expression(sigma^2))
  plot(out$param[,2],type="l",ylab=expression(phi))
  plot(out$param[,3],type="l",ylab=expression(tau))
par(mfrow=c(1,1))

out$acc_rate
plot.post(out$param[,1],stay=T)
plot.post(out$param[,2],stay=T)
plot.post(out$param[,3],stay=T)

apply(out$param,2,mean)
apply(out$param,2,sd)
apply(out$param,2,quantile)


