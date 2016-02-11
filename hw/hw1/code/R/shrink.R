set.seed(268)

library(MASS) # ridge.lm
library(glmnet) #cv.glmnet(x,y): lasso
library(doMC)
registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../C++/spikeAndSlab.cpp")
sourceCpp("../C++/blasso.cpp")
sourceCpp("../C++/gdp.cpp")
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
#n <- c(50, 100)
#p <- c(100,200)
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
  ridge.mod <- lm.ridge(y~x[,-1])
  spike.mod <- spikeAndSlab(y=y, x=x, tau2=rep(1e-6,ncol(x)), g=1e8, w=rep(.5,ncol(x)), B=2000, burn=500, printProg=F, returnHyper=F)
  blasso.mod <- blasso(y=y,x=x,r=1,delta=1.5,B=1200,burn=200, printProg=F, returnHyper=F)
  gdp.mod <- gdp(y=y, x=x, B=2000, burn=500, printProg=F,returnHyper=F)
  mod <- list("lasso_mod"=lasso.mod, "ridge_mod"=ridge.mod, "spike_mod"=spike.mod$beta, "blasso_mod"=blasso.mod$beta, "gdp_mod"=gdp.mod$beta, "param_index"=param_index)
  mod
}

mod.params.ind <- as.list(1:numdat); counter <- 0
for (n_i in 1:length(n)) for (p_i in 1:length(p)) for (S_i in 1:length(Sig_Makers)) for (beta_i in 1:length(betas)) {
  counter <- counter + 1
  mod.params.ind[[counter]] <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
}

mod <- foreach(mpi=mod.params.ind) %dopar% oneSim(mpi$n,mpi$p,mpi$S,mpi$beta)

mod_n <- 35 
plot( beta_lookup(mod[[mod_n]]$param_index), col="grey", pch=1, ylim=c(-5,5), bty="n", cex =2, xlab="Parameter Index", ylab=expression(beta));
points( apply(mod[[mod_n]]$spike,2,mean),col='red',pch=20 );
points( apply(mod[[mod_n]]$gdp,2,mean),col='green',pch=20 );
points( apply(mod[[mod_n]]$blasso,2,mean),col='blue',pch=20 );

points( coef(mod[[mod_n]]$lasso),col='blue',pch=20 );
points( coef(mod[[mod_n]]$ridge),col='red',pch=20 );


# - For all theses Bayesian models, get E(\beta_j | y)
rmse <- function(b,b_true) sqrt(mean((b-b_true)^2))
compareBayesian <- function(model,rmse_ord="blasso",ylim_rmse=c(0,5),cex_rmse=1) {
  cc <- 0
  rmse_blasso <- NULL; rmse_spike <- NULL; rmse_gdp <- NULL;
  for (mm in model) {
    cc <- cc + 1
    mod_ind <- mm$param_index
    beta_true <- beta_lookup(mod_ind)
    rmse_blasso[cc] <- rmse(mm$blasso,beta_true)
    rmse_spike[cc] <- rmse(mm$spike,beta_true)
    rmse_gdp[cc] <- rmse(mm$gdp,beta_true)
  }

  ord <- 1:length(mod)
  if (rmse_ord == "blasso") {
    ord <- order(rmse_blasso)
  } else if (rmse_ord == "gdp") {
    ord <- order(rmse_gdp)
  } else if (rmse_ord == "spike") {
    ord <- order(rmse_spike)
  }

  par(mar=c(4,4,0,0))
  plot(rmse_blasso[ord],col="blue",pch=1,ylim=ylim_rmse,cex=cex_rmse,
       ylab="RMSE",xlab="",lwd=2,xaxt="n",bty="n")
  points(rmse_spike[ord],col="red",pch=2,cex=cex_rmse,lwd=2)
  points(rmse_gdp[ord],col="green",pch=4,cex=cex_rmse,lwd=2)

  simdat <- t(sapply(ord, function(o) mod[[o]]$param_index))
  simdat <- cbind(
    ifelse(simdat[,1] == 1, "n50", "n500"),
    ifelse(simdat[,2] == 1, "p100", "p5000"),
    ifelse(simdat[,3] == 1, "I", ifelse(simdat[,3] == 2, "S.1", "S.6")),
    ifelse(simdat[,4] == 1, "b1", ifelse(simdat[,4] == 2, "b2", "b3")))
  lab <- apply(simdat,1,function(x) paste(x,collapse="\n"))
  axis(1,at=1:36,label=lab,pos=-.4,tck=0,lty=0,cex.axis=.6)
  par(mar=c(5.1,4.1,4.1,2.1))
  abline(v=5*(1:7),col=rgb(.5,.5,.5))

  list("ord"=ord, "blasso"=rmse_blasso, "gdp"=rmse_gdp, "spike"=rmse_spike)
}

pdf("output/rmseblasso.pdf",w=16,h=9)
  rmse_blasso <- compareBayesian(mod,"blasso",ylim=c(0,10),cex=1.5)
dev.off()
pdf("output/rmsessvn.pdf",w=16,h=9)
rmse_spike  <- compareBayesian(mod,"spike",ylim=c(0,10),cex=1.5)
dev.off()
pdf("output/rmsegdp.pdf",w=16,h=9)
rmse_gdp    <- compareBayesian(mod,"gdp",ylim=c(0,10),cex=1.5)
dev.off()

sapply(rmse_gdp,mean) # blasso seems to perform a little better than gdp and a lot better than spike and slab
#      ord    blasso       gdp     spike
#18.500000  1.114592  1.124779  3.695219

# - discuss accuracy wrt to a metric
# - Compare Lasso and S&S for variable selection
# - Compute mean(L_j: beta_j == 0)
# - posterior pred.


