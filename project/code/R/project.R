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


### p = 3
f3 <- function(X) apply(as.matrix(X),1, function(x) x[1] + ifelse(x[2]>.5,x[2]-.5,0) + x[3]^2)
dat <- sim_dat(f3,p=3,n=50)
#priors <- c(4.5,17.5,  .1,5,  .1,.1) #s2, phi, tau
priors <- c(5,5,  .1,5,  .1,.1) #s2, phi, tau
out <- gp(y=dat$y, x=dat$x, D=dat$D, cand_S=diag(dat$p)*.1,
          init=rep(1,3), priors=priors, B=2000, burn=4000, printProg=TRUE)
colnames(out$param) <- c("s2","phi","tau")

### p = 500
f3.2 <- function(X) apply(as.matrix(X),1, function(x) 10*x[1] + 15*sin(x[2]) + 10*x[3]^2)
dat2 <- sim_dat(f3.2,p=500,n=100)
priors <- c(2,1,    .1,5,    2,5,    1,1) #s2, phi, tau, d_vec
cat("sourcing gp_gdp.cpp...\n"); sourceCpp("../C++/gp_gdp.cpp")
out2 <- gp_gdp(y=dat2$y, X=dat2$x, cand_S=diag(3+dat2$p)*1e-9,
               init=rep(0,3+dat2$p), priors=priors, B=2000, burn=100000, printProg=TRUE)
save(out2,file="output/out2.RData")
load("output/out2.RData")


#Plots: #####################################################3

### p=3
plot.posts(out$param,cex.l=1.3,cex.a=1,names=colnames(out$param))
plot(ts(out$param))
apply(out$param,2,summary)


### p=500
plot(ts(out2$param))
plot.posts(out2$param[,1:3],cex.l=1.3,cex.a=1)
plot( apply(out2$param[,-c(1:3)],2,mean) )
apply(out2$param,2,summary)
