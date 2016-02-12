set.seed(268)

library(MASS) # ridge.lm
library(glmnet) #cv.glmnet(x,y): lasso
library(doMC)
registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../C++/spikeAndSlab.cpp")
system("mkdir -p output")

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

n <- c(50, 500)
p <- c(100, 1000)
Sig_Makers <- list( function(p) diag(p), function(p) Sig_Maker(.1,p), function(p) Sig_Maker(.6,p) )
betas <- list( function(p) {b <- rep(0,p); b[1:5] <- 3; b <- c(0,b); b }, 
               function(p) {b <- rep(0,p); b[1:5] <- 5; b[6:10] <- -2; b[11:15] <- .5; b <- c(0,b); b },
               function(p) {b <- rep(1,p); b <- c(1,b); b })

numdat <- length(n) * length(p) * length(Sig_Makers) * length(betas)
counter <- 0
mod <- as.list(1:numdat)

beta_lookup <- function(param_index) { 
  bb <- param_index$beta
  pp <- param_index$p
  betas[[bb]](p[pp]) 
}

oneSim <- function(n_i,p_i,S_i,beta_i) {
  param_index <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
  print( paste0( c("n: ","p: ","S: ","beta: "), unlist(param_index) ))
  x <- t(sapply( 1:(n[n_i]), function(x) c(1, mvrnorm(0,Sig_Makers[[S_i]](p[p_i]))) ))
  y <- x %*% betas[[beta_i]](p[p_i]) + rnorm(n[n_i])
  lasso.mod <- cv.glmnet(x[,-1],y)
  ssvn.mod <- spikeAndSlab(y=y, x=x, tau2=rep(1e-6,ncol(x)), g=1e8, w=rep(.5,ncol(x)), B=2000, burn=500, printProg=F, returnHyper=T)
  mod <- list("lasso_mod"=lasso.mod, "ssvn_mod"=ssvn.mod, "param_index"=param_index)
  mod
}

mod.params.ind <- as.list(1:numdat); counter <- 0
for (n_i in 1:length(n)) for (p_i in 1:length(p)) for (S_i in 1:length(Sig_Makers)) for (beta_i in 1:length(betas)) {
  counter <- counter + 1
  mod.params.ind[[counter]] <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
}

mod <- foreach(mpi=mod.params.ind) %dopar% oneSim(mpi$n,mpi$p,mpi$S,mpi$beta)

mod_n <- 4
true_beta <- beta_lookup(mod[[mod_n]]$param_index)
postmean.gamma <- apply(mod[[mod_n]]$ssvn$gam,2,mean)
sel_lasso <- coef(mod[[mod_n]]$lasso) != 0
sel_ssvn <- postmean.gamma > .5
mean(as.numeric(sel_lasso) == (true_beta != 0))
mean(as.numeric(sel_ssvn) == (true_beta != 0))

#sourceCpp("../C++/spikeAndSlab.cpp")
#x <- t(sapply( 1:(n[1]), function(x) c(1, mvrnorm(0,Sig_Makers[[1]](p[1]))) ))
#y <- x %*% betas[[2]](p[1]) + rnorm(n[1])
#tt <- .05^2
#ssvn.mod <- spikeAndSlab(y=y, x=x, tau2=rep(tt,ncol(x)), g=100 / tt, w=rep(.5,ncol(x)), B=5000, burn=2000, printProg=T, returnHyper=T)
#round(apply(ssvn.mod$beta,2,mean),2)
#round(apply(ssvn.mod$g,2,mean),2)
#
#plot(ssvn.mod$b[,14],type='l',col='red')
#abline(h=betas[[2]](p[1])[14],col='grey')
#
#plot( round(apply(ssvn.mod$beta,2,mean),2) )
#points(betas[[2]](p[1]),col='grey',pch=20)
