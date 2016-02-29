set.seed(268)
system("mkdir -p output")
source("plotmap.R")
source("../../../../R_Functions/plotPost.R",chdir=T)
library(spBayes)

n <- 1000
sig2 <- .5
sn <- 30

# TESTING
p <- 2
x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
f <- function(xx) ifelse(xx[,1]>.5, xx[,1]-.5, 0) + xx[,2]^2
mu <- f(x)
s <- matrix(runif(sn*p,range(x)[1],range(x)[2]),sn,p) # knots

y <- rnorm(n,f(x),sqrt(sig2)) # data (simulated responses)

D <- as.matrix(dist(s))
Cd <- matrix(0,n,sn) # distance between data and knots
for (i in 1:n) for (j in 1:sn) Cd[i,j] <- sqrt(sum((x[i,] - s[j,])^2))

# y | ... ~ N(0,s^2 + K)
#priors <- c(2,.5,  4.9,5.1,  2,2) #s2, phi, tau
priors <- list("sigma.sq.IG"=c(2,2), # the scaling for the covariance matrix. 
               "phi.Unif"=c(0.1,5),  # decay parameter
               "tau.sq.IG"=c(2,.5),  # the variance between observations. Should be .5
               "beta.norm"=list(0,.00001))
starting <- list("phi"=1,"sigma.sq"=1,"tau.sq"=1)
tuning <- list("phi"=.01,"sigma.sq"=.01,"tau.sq"=.01)
out <- spLM(y~1,coords=x,knots=s,cov.model="exponential",n.samples=12000,priors=priors,starting=starting,tuning=tuning,verbose=T,modified.pp=F)


apply(as.matrix(out$p.theta.samples),2,mean) # 16, .5, .12
apply(as.matrix(out$p.theta.samples),2,quantile) # 16, .5, .12
plot.posts(tail(as.matrix(out$p.theta.samples),2000)) # scaling for cov, decay, obs var
plot.post(tail(as.matrix(out$p.beta),2000),stay=T)
m1 <- spRecover(out,start=10001,verbose=TRUE)
preds <- m1$p.w.recover.samples

#col.map <- colorRampPalette(c('darkred','orange','yellow'),bias=2)(length(mu))
col.map <- colorRampPalette(c('white','yellow','gold','orange','darkred'),bias=2)(length(mu))
col.diff <- colorRampPalette(c('darkblue','lightblue','white','yellow','darkred'))(length(mu))

par(mfrow=c(2,2))
  plotmap(f(m1$coords),m1$coords, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(apply(preds,1,mean),m1$coords, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(f(m1$coords)-apply(preds,1,mean),m1$coords, bks=c(-1,1)*.5,xlim=c(-3,3),ylim=c(-3,3),col.map=col.diff,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(apply(preds,1,sd),m1$coords, bks=c(0,2),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
par(mfrow=c(1,1))


# How do I get wStar??? ######################
onePred_mu_star <- function(param) {
  tau <- param[1]
  s2  <- param[2]
  phi <- param[3]

  Ks <- tau * exp(-phi * D)
  Ks.i <- solve(Ks)
  C <- tau * exp(-phi * Cd)
  H <- C %*% Ks.i
  Ht <- t(H)

  S.i <- solve( Ks.i + Ht%*%H / s2 )
  m <- S.i %*% Ht %*% y / s2

  mu_star <- mvrnorm(1,m,S.i)
  mu_star
}

R <- apply(as.matrix(m1$p.theta.recover.samples),1,function(param) onePred_mu_star(param))

onePred_mu_star(m1$p.theta.recover.samples[3,])

plot(apply(m1$p.wStr,1,mean),col="blue",pch=20,cex=2)
points(apply(R,1,mean),pch=20,col="red")
