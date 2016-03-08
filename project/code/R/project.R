set.seed(268)
library(Rcpp)
system("mkdir -p output")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
source("sim_dat.R")

cat("sourcing func.cpp...\n"); sourceCpp("../../../cpp_functions/func.cpp")
cat("sourcing gp.cpp...\n"); sourceCpp("../C++/gp.cpp")
cat("sourcing gp_gdp.cpp...\n"); sourceCpp("../C++/gp_gdp.cpp")
cat("Starting Main Program...\n")

f3 <- function(X) apply(as.matrix(X),1, function(x) x[1] + ifelse(x[2]>.5,x[2]-.5,0) + x[3]^2)
f10 <- function(X) apply(as.matrix(X),1, function(x) x[1] + ifelse(x[2]>.5,x[2]-.5,0) + x[3]^2)

dat <- sim_dat(f3,p=3,n=1000)
priors <- c(2,.5,  .9,1.1,  2,1) #s2, phi, tau
out <- gp(y=dat$y, x=dat$x, D=dat$D, cand_S=diag(dat$p)*1,
          init=rep(1,3), priors=priors, B=1000, burn=3000, printProg=TRUE)
colnames(out$param) <- c("s2","phi","tau")
plot.posts(out$param,cex.l=1.3,cex.a=1,names=colnames(out$param))
apply(out$param,2,quantile)

plot(ts(out$param))
