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
sn <- 50

# ORIGINAL:
p <- 3
x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
f <- function(xx) x[,1] + ifelse(xx[,2]>.5, xx[,2]-.5, 0) + xx[,3]^2
mu <- f(x)
ord <- order(mu)
s <- matrix(runif(sn*p,range(x)[1],range(x)[2]),sn,p) # knots
#s <- matrix(rnorm(sn*p),sn,p) # knots

# TESTING
#p <- 2
#x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
#f <- function(xx) ifelse(xx[,1]>.5, xx[,1]-.5, 0) + xx[,2]^2
#mu <- f(x)
#s <- matrix(runif(sn*p,range(x)[1],range(x)[2]),sn,p) # knots

y <- rnorm(n,f(x),sqrt(sig2)) # data (simulated responses)

D <- as.matrix(dist(s))
Cd <- matrix(0,n,sn)
for (i in 1:n) for (j in 1:sn) Cd[i,j] <- sqrt(sum((x[i,] - s[j,])^2))

# y | ... ~ N(0,s^2 + K)
priors <- c(2,.5,  .1,3,  2,2) #s2, phi, tau
system.time( out <- gp(y, x, s, Cd, D, cand_S=diag(3)*.01, init=rep(0,3), priors=priors, B=300, burn=400, printProg=T) )

save(out,file="output/out.RData")

colnames(out$param) <- c("s2","phi","tau")
par(mfrow=c(3,1))
  plot(out$param[,1],type="l",ylab=expression(sigma^2))
  plot(out$param[,2],type="l",ylab=expression(phi))
  plot(out$param[,3],type="l",ylab=expression(tau))
par(mfrow=c(1,1))

out$acc_rate
plot.posts(out$param,names=c("s2","phi","tau"))


apply(out$param,2,mean)
apply(out$param,2,sd)
apply(out$param,2,quantile)


onePred_mu_star <- function(param,o,retList=F) {
  s2  <- param[1]
  phi <- param[2]
  tau <- param[3]

  Ks <- tau * exp(-phi * o$D)
  Ks.i <- solve(Ks)
  C <- tau * exp (-phi * o$Cd)
  H <- C %*% Ks.i
  Ht <- t(H)

  S.i <- solve( Ks.i + Ht%*%H / s2 )
  m <- S.i %*% Ht %*% o$y / s2

  mu_star <- mvrnorm(m,S.i)
  out <- NULL

  if (retList) out <- list("mu_star"=mu_star,"C"=C,"Ks"=Ks) 
  else out <- mu_star
  
  out
}


onePred <- function(param,o) {
  ms <- onePred_mu_star(param,o,T)
  mu <- ms$C %*% solve(ms$Ks) %*% ms$mu_star
  mu
}

system.time( preds <- t(apply(out$param,1,function(p) onePred(p,out))) )
plot(apply(preds,2,mean)[ord],type='l',col="blue",lwd=2)
lines(mu[ord],col="grey",lwd=2)

plot(0,cex=0,ylim=c(-3,6),xlim=c(1,n))
apply(preds[,ord],1,function(r) points(r,type='p',ylim=c(-3,6),col=rgb(.1,.1,.1,.1),cex=.1))
lines(apply(preds,2,mean)[ord],ylim=c(-3,6),col="blue")
lines(mu[ord],col='green',lwd=3)

plot(mu - apply(preds,2,mean))


#Map Plots for testing######################
#source("plotmap.R")
#col.map <- colorRampPalette(c('white','yellow','gold','orange','darkred'),bias=2)(n)
#col.map.s <- colorRampPalette(c('darkred','orange','yellow'),bias=2)(sn)
#col.diff <- colorRampPalette(c('darkblue','lightblue','white','yellow','darkred'))(n)
#
#
#par(mfrow=c(2,2))
#  plotmap(f(out$x),out$x, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
#  plotmap(apply(preds,2,mean),out$x, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
#  plotmap(f(out$x)-apply(preds,2,mean),out$x, bks=c(-1,1)*.5,xlim=c(-3,3),ylim=c(-3,3),col.map=col.diff,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
#  plotmap(apply(preds,2,sd),out$x, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
#par(mfrow=c(1,1))


#system.time( preds_ms <- t(apply(out$param,1,function(p) onePred_mu_star(p,out))) )
#par(mfrow=c(1,2))
#  plotmap(f(out$s),out$s,bks=c(0,3),xlim=c(-2,3),ylim=c(-3,3),col.map=col.map.s)
#  plotmap(apply(preds_ms,2,mean),out$s,bks=c(0,3),xlim=c(-2,3),ylim=c(-3,3),col.map=col.map.s)
#par(mfrow=c(1,1))
