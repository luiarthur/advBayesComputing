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
system.time( out <- gp(y, x, s, Cd, D, cand_S=diag(3)*.01, init=rep(0,3), priors=priors, B=2000, burn=10000, printProg=T) )

save(out,file="output/out.RData")

colnames(out$param) <- c("s2","phi","tau")
par(mfrow=c(3,1))
  plot(out$param[,1],type="l",ylab=expression(sigma^2))
  plot(out$param[,2],type="l",ylab=expression(phi))
  plot(out$param[,3],type="l",ylab=expression(tau))
par(mfrow=c(1,1))

out$acc_rate
plot.posts(out$param,names=c(expression(sigma^2),expression(phi),expression(tau^2)),cex.a=1.5,cex.legend=1.2)


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

# mu^ and mu in 1-D
pdf("output/gpOrderedData.pdf",w=13,h=9)
  plot(apply(preds,2,mean)[ord],type='l',col="blue",lwd=.5,bty="n",ylab=expression(mu),xlab="Ordered Data")
  lines(mu[ord],col="black",lwd=3)
dev.off()

# Residuals
pdf("output/gpResid.pdf",w=13,h=9)
  plot((mu - apply(preds,2,mean))[ord],pch=20,main=expression(mu-hat(mu)),ylab="Residuals",xlab="Ordered Data")
dev.off()

source("plotmap.R")
rmse_map <- function(X,x_true) {
  rmse <- function(x,x_true) {
    sqrt(sum((x-x_true)^2)/length(x))
  }
  apply(matrix(1:ncol(X)),1,function(i) rmse(X[,i],x_true[i]))
}

col.map <- colorRampPalette(c("blue","grey90","red"))(n)
preds.mean <- apply(preds,2,mean)

pdf("output/plot3d.pdf",w=13,h=9)
  par(mfrow=c(2,2))
    plotmap3d(x[,3],x[,2],x[,1],xlab="x3",ylab="x2",zlab="x1",clim.map=c(-.9,2),pch=20,cex=1,theta=45,phi=20,col.map=col.map,cex.main=2,main=expression(hat(mu)),val=preds.mean)
    plotmap3d(x[,3],x[,2],x[,1],xlab="x3",ylab="x2",zlab="x1",clim.map=c(-.9,2),pch=20,cex=1,theta=45,phi=20,col.map=col.map,cex.main=2,main=expression(mu),val=mu)
    plotmap3d(x[,3],x[,2],x[,1],xlab="x3",ylab="x2",zlab="x1",clim.map=c(-.9,2),pch=20,cex=1,theta=45,phi=20,col.map=col.map,cex.main=2,main=expression(E~"["~(hat(mu)-mu)^2~"|"~mu~"]"),val=rmse_map(preds,mu))
    plotmap3d(x[,3],x[,2],x[,1],xlab="x3",ylab="x2",zlab="x1",clim.map=c(-.9,2),pch=20,cex=1,theta=45,phi=20,col.map=col.map,cex.main=2,main=expression(sd(hat(mu))),val=apply(preds,2,sd))
  par(mfrow=c(1,1))
dev.off()




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
