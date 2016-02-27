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

y <- rnorm(n,f(x),sqrt(sig2)) # data (simulated responses)
C <- cov(t(x),t(s)) # covariance between data and knots
D <- as.matrix(dist(s))

# y | ... ~ N(0,s^2 + K)
priors <- c(2,5,0,5,2,5) #s2, phi, tau
#prelim <- gp(y, x, s, C, D, cand_S=.01*diag(3), init=rep(0,3), priors=priors, B=500, burn=50, printProg=T)
#V <- round( cov(prelim$param) )
#diag(V) <- ifelse(diag(V)<1,sign(diag(V)),diag(V))
#system.time( out <- gp(y, x, s, C, D, cand_S=V, init=tail(prelim$param,1), B=3000, burn=5000, printProg=T) )
system.time( out <- gp(y, x, s, C, D, cand_S=.01*diag(3), init=rep(0,3), priors=priors, B=3000, burn=5000, printProg=T) )

save(out,file="output/out.RData")

colnames(out$param) <- c("s2","phi","tau")
par(mfrow=c(3,1))
  plot(out$param[,1],type="l",ylab=expression(sigma^2))
  plot(out$param[,2],type="l",ylab=expression(phi))
  plot(out$param[,3],type="l",ylab=expression(tau))
par(mfrow=c(1,1))

out$acc_rate
plot.post(out$param[,1],stay=T,main=expression(sigma^2))
plot.post(out$param[,2],stay=T,main=expression(phi))
plot.post(out$param[,3],stay=T,main=expression(tau))

apply(out$param,2,mean)
apply(out$param,2,sd)
apply(out$param,2,quantile)

plot.posts(out$param)

onePred_mu_star <- function(param,o,retList=F) {
  phi <- param[2]
  tau <- param[3]
  Cs <- o$C
  Ds <- o$D
  
  Ks <- tau * exp(-phi * Ds)
  ks <- ncol(Ks)

  mu_star <- mvrnorm(rep(0,ks), Ks)
  out <- NULL

  if (retList) out <- list("mu_star"=mu_star,"Cs"=Cs,"Ks"=Ks) 
  else out <- mu_star
  
  out
}


onePred <- function(param,o) {
  ms <- onePred_mu_star(param,o,T)
  mu <- ms$Cs %*% solve(ms$Ks) %*% ms$mu_star
  mu
}

#onePred(out$param[1,],out)

system.time( preds <- t(apply(out$param,1,function(p) onePred(p,out))) )
plot(apply(preds,2,mean),type='l',col="blue",lwd=2)

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
col.map.s <- colorRampPalette(c('darkred','orange','yellow'),bias=2)(length(f(s)))
col.map.diff <- colorRampPalette(c('darkred','white','yellow'))(length(f(s)))
par(mfrow=c(1,3))
  plotmap(mu,x,bks=c(-.1,3.1),xlim=c(-2,2),ylim=c(-3,3),col.map=col.map)
  plotmap( y,x,bks=c(-.1,3.1),xlim=c(-2,2),ylim=c(-3,3),col.map=col.map)
  plotmap(apply(preds,2,mean),bks=c(-.1,3.1),x,col.map=col.map,ylim=c(-3,3))
par(mfrow=c(1,1))

system.time( preds_ms <- t(apply(out$param,1,function(p) onePred_mu_star(p,out))) )
par(mfrow=c(1,2))
  plotmap(f(s),s,bks=c(0,3),xlim=c(-2,2),ylim=c(-3,3),col.map=col.map.s)
  plotmap(apply(preds_ms,2,mean),s,bks=c(-1,3),xlim=c(-2,2),ylim=c(-3,3),col.map=col.map.s)
par(mfrow=c(1,1))

  plotmap(f(s)-apply(preds_ms,2,mean),s,xlim=c(-2,2),ylim=c(-3,3),col.map=col.map.diff)


#x2 <- sort( rnorm(5000,0,10) )
#X2 <- t(apply(matrix(1:500),1,function(i)rnorm(5000,g(x2),sqrt(.5))))
#plot(x2,X2[1,],col=rgb(.1,.1,.1,0),cex=1)
#apply(X2,1,function(row) points(x2,row,col=rgb(.1,.1,.1,.1),cex=.1))
#lines(x2,apply(X2,2,mean),col="blue",lwd=2)
#lines(x2,apply(X2,2,function(x) mean(x) + sd(x)),col="blue",lwd=1)
#lines(x2,apply(X2,2,function(x) mean(x) - sd(x)),col="blue",lwd=1)
