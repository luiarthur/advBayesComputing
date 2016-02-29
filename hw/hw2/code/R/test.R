set.seed(268)
system("mkdir -p output")
source("plotmap.R")
source("../../../../R_Functions/plotPost.R",chdir=T)
library(spBayes)

spLM

n <- 1000
sig2 <- .5
sn <- 300

# TESTING
p <- 2
x <- matrix(rnorm(n*p),n,p)     # data (simulated covariates)
f <- function(xx) ifelse(xx[,1]>.5, xx[,1]-.5, 0) + xx[,2]^2
mu <- f(x)
s <- matrix(runif(sn*p,range(x)[1],range(x)[2]),sn,p) # knots

y <- rnorm(n,f(x),sqrt(sig2)) # data (simulated responses)
C <- cov(t(x),t(s)) # covariance between data and knots
D <- as.matrix(dist(s))

# y | ... ~ N(0,s^2 + K)
#priors <- c(2,.5,  4.9,5.1,  2,2) #s2, phi, tau
priors <- list("sigma.sq.IG"=c(2,2), # the scaling for the covariance matrix. 
               "phi.Unif"=c(0.1,5),  # decay parameter
               "tau.sq.IG"=c(2,.5),  # the variance between observations. Should be .5
               "beta.norm"=list(0,.00001))
starting <- list("phi"=1,"sigma.sq"=1,"tau.sq"=1)
tuning <- list("phi"=.01,"sigma.sq"=.01,"tau.sq"=.01)
out <- spLM(y~1,coords=x,knots=s,cov.model="exponential",n.samples=12000,priors=priors,starting=starting,tuning=tuning,verbose=T,modified.pp=F)


plot.posts(as.matrix(out$p.theta.samples))
m1 <- spRecover(out,start=10001,verbose=TRUE)
preds <- m1$p.w.recover.samples

col.map <- colorRampPalette(c('darkred','orange','yellow'),bias=2)(length(mu))
col.diff <- colorRampPalette(c('darkblue','white','darkred'))(length(mu))

par(mfrow=c(2,2))
  plotmap(f(m1$coords),m1$coords, bks=c(0,1),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(apply(preds,1,mean),m1$coords, bks=c(0,1),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(f(m1$coords)-apply(preds,1,mean),m1$coords, bks=c(-1,1)*.2,xlim=c(-3,3),ylim=c(-3,3),col.map=col.diff,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
  plotmap(apply(preds,1,sd),m1$coords, bks=c(0,1),xlim=c(-3,3),ylim=c(-3,3),col.map=col.map,ylab="x2",xlab="x1"); abline(v=c(.5),col="grey")
par(mfrow=c(1,1))

