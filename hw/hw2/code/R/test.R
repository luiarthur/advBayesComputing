set.seed(268)
system("mkdir -p output")
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
priors <- list("sigma.sq.IG"=c(2,1/.5),
               "phi.Unif"=c(0.1,5),
               "tau.sq.IG"=c(2,1/2))
starting <- list("phi"=1,"sigma.sq"=1,"tau.sq"=1)
tuning <- list("phi"=.01,"sigma.sq"=.01,"tau.sq"=.01)
out <- spLM(y~rep(0,n),coords=x,knots=s,cov.model="exponential",n.samples=300,priors=priors,starting=starting,tuning=tuning,verbose=T,modified.pp=F)
#system.time( out <- gp(y, x, s, C, D, cand_S=.01*diag(3), init=rep(0,3), priors=priors, B=3000, burn=15000, printProg=T) )


