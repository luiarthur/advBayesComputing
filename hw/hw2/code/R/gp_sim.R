set.seed(268)
library(Rcpp)
system("mkdir -p output")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../../R_Functions/plotPost.R",chdir=T)
cat("sourcing func.cpp...\n")
sourceCpp("../../../../cpp_functions/func.cpp")
cat("sourcing gp.cpp...\n")
sourceCpp("../C++/gp.cpp")
cat("Starting Main Program...\n")

n <- 1000
sig2 <- .5
sn <- 100

# ORIGINAL:
#p <- 3
#x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
#mu <- x[,1] + ifelse(x[,2]-.5 > 0, x[,2]-.5, 0) + x[,3]^2
#s <- matrix(rnorm(sn*p),sn,p) # knots

# TESTING
p <- 2
x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
f <- function(xx) ifelse(xx[,1]-.5 > 0, xx[,1]-.5, 0) + xx[,2]^2
mu <- f(x)
s <- matrix(rnorm(sn*p),sn,p) # knots

# Sorting
ord <- order(mu)
mu <- sort(mu)
x <- x[ord,]
# End Sorting
y <- mu + rnorm(n,0,sqrt(sig2)) # data (simulated responses)
C <- cov(t(x),t(s)) # covariance between data and knots
D <- as.matrix(dist(s))^2

# y | ... ~ N(0,s^2 + K)
prelim <- gp(y, x, s, C, D, cand_S=diag(3), init=rep(0,3), B=500, burn=500, printProg=T)
V <- cov( prelim$param )
system.time( out <- gp(y, x, s, C, D, cand_S=V, init=tail(prelim$param,1), B=2000, burn=500, printProg=T) )

colnames(out$param) <- c("s2","phi","tau")
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


onePred <- function(param,o) {
  phi <- param[2]
  tau <- param[3]
  Cs <- o$C
  Ds <- o$D
  
  Ks <- tau * exp(-phi * Ds)
  ks <- ncol(Ks)

  mu_star <- mvrnorm(rep(0,ks), Ks)
  mu <- Cs %*% solve(Ks) %*% mu_star

  mu
}

system.time( preds <- t(apply(out$param,1,function(p) onePred(p,out))) )

plot(apply(preds,2,mean),type='l',ylim=c(-30,30),col="blue",lwd=2)
lines(apply(preds,2,function(x)quantile(x,.025)),type='l',ylim=range(mu),col="blue",lwd=2)
lines(apply(preds,2,function(x)quantile(x,.975)),type='l',ylim=range(mu),col="blue",lwd=2)
lines(mu,lwd=3,col='grey')


### MLE?
M <- 3 * exp(-3*D)
ms <- mvrnorm(rep(0,ncol(M)), M)
mu_pred <- C %*% solve(M) %*% ms

plot(mu,lwd=3,col="grey",pch=20,ylim=c(-10,10))
points(mu_pred,pch=20,ylim=range(mu))
points(mu-mu_pred,pch=20,col="red")

#Map Plots for testing######################
source("plotmap.R")
col.map <- colorRampPalette(c('darkred','orange','yellow'),bias=2)(length(mu))
plotmap(mu,x,bks=c(0,3),xlim=c(-2,2),ylim=c(-3,3),col.map=col.map)


