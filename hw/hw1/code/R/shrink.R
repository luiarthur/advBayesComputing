set.seed(268)

library(MASS) # ridge.lm
library(glmnet) #cv.glmnet(x,y): lasso
library(doMC)
registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../C++/spikeAndSlab.cpp")
sourceCpp("../C++/blasso.cpp")


mvrnorm <- function(M,S,n=nrow(S)) M + t(chol(S)) %*% rnorm(n)

Sig_Maker <- function(rho,p) {
  M <- matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      M[i,j] <- rho^{abs(i-j)}
    }
  }
  M
}

#n <- c(50, 500)
#p <- c(100, 1000)
n <- c(50, 100)
p <- c(100,200)
Sig_Makers <- list( function(p) diag(p), function(p) Sig_Maker(.1,p), function(p) Sig_Maker(.6,p) )
betas <- list( function(p) {b <- rep(0,p); b[1:5] <- 3; b <- c(0,b); b }, 
               function(p) {b <- rep(0,p); b[1:5] <- 5; b[6:10] <- -2; b[11:15] <- .5; b <- c(0,b); b },
               function(p) {b <- rep(1,p); b <- c(1,b); b })

numdat <- length(n) * length(p) * length(Sig_Makers) * length(betas)
counter <- 0
mod <- as.list(1:numdat)


# - For all theses Bayesian models, get E(\beta_j | y)
# - discuss accuracy wrt to a metric
# - Compare Lasso and S&S for variable selection
# - Compute mean(L_j: beta_j == 0)
# - posterior pred.

# Temp
#Sig <- Sig_Makers[[1]](100);
#beta <- betas[[1]](100)

oneSim <- function(n_i,p_i,S_i,beta_i) {
  param_index <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
  print( paste0( c("n: ","p: ","S: ","beta: "), unlist(param_index) ))
  x <- t(sapply( 1:(n[[n_i]]), function(x) c(1, mvrnorm(0,Sig_Makers[[S_i]](p[[p_i]]))) ))
  y <- x %*% betas[[beta_i]](p[[p_i]]) + rnorm(n[[n_i]])
  lasso.mod <- cv.glmnet(x[,-1],y)
  ridge.mod <- lm.ridge(y~x[,-1])
  spike.mod <- spikeAndSlab(y=y, x=x, tau2=rep(1e-6,ncol(x)), g=1e8, w=rep(.5,ncol(x)), B=2000, burn=500, printProg=F)
  blasso.mod <- blasso(y=y,x=x,r=1,delta=1.5,B=1200,burn=200, printProg=F)
  gdp.mod <- "not complete"
  mod <- list("lasso_mod"=lasso.mod, "ridge_mod"=ridge.mod, "spike_mod"=spike.mod, "blasso_mod"=blasso.mod, "gdp_mod"=gdp.mod, "param_index"=param_index)
  mod
}

mod.params.ind <- as.list(1:numdat); counter <- 0
for (n_i in 1:length(n)) for (p_i in 1:length(p)) for (S_i in 1:length(Sig_Makers)) for (beta_i in 1:length(betas)) {
  counter <- counter + 1
  mod.params.ind[[counter]] <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
}

mod <- foreach(mpi=mod.params.ind) %dopar% oneSim(mpi$n,mpi$p,mpi$S,mpi$beta)

# Bayesian Lasso WORKING!
#sourceCpp("../C++/blasso.cpp")
#blasso.mod <- blasso(y=y,x=x,r=1,delta=1.5,B=1200,burn=200, printProg=T)
#bl.b <- blasso.mod$beta
#plot(apply(bl.b,2,mean))
#plot(apply(blasso.mod$t2,2,mean))
#plot(blasso.mod$lam)
#hist(blasso.mod$lam)

# Spike and Slab WORKING!
#sourceCpp("../C++/spikeAndSlab.cpp")
#B <- 1500
#burn <- 500
#ss.mod <- spikeAndSlab(y=y, x=x, tau2=rep(1e-6,ncol(x)), g=1e8, w=rep(.5,ncol(x)), B=B+burn, burn=burn, printProg=T)
#ss.b <- ss.mod$beta
#ss.g <- ss.mod$gam
#(post.gam <- round(apply(ss.g,2,mean),5))
#(post.beta <- round(apply(ss.b, 2, mean),5))
#order(post.beta,decreasing=T)
#par(mfrow=c(3,1)); plot(ss.b[,1]); plot(ss.b[,3]); plot(ss.b[,10]); par(mfrow=c(1,1))
