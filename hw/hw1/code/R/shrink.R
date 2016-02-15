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

n <- c(50,   500)
p <- c(100, 1000)
#n <- c(50,100)
#p <- c(50,100)
Sig_Makers <- list( function(p) diag(p), function(p) Sig_Maker(.1,p), function(p) Sig_Maker(.6,p) )
betas <- list( function(p) {b <- rep(0,p); b[1:5] <- 3; b <- c(0,b); b }, 
               function(p) {b <- rep(0,p); b[1:5] <- 5; b[6:10] <- -2; b[11:15] <- .5; b <- c(0,b); b },
               function(p) {b <- rep(1,p); b <- c(0,b); b })

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
  tt <- .05^2
  spike.mod <- spikeAndSlab(y=y, x=x, tau2=rep(tt,ncol(x)), g=100/tt, w=rep(.5,ncol(x)), B=3000, burn=2000, printProg=F, returnHyper=T)
  blasso.mod <- blasso(y=y,x=x,r=1,delta=1.5,B=1200,burn=200, printProg=F, returnHyper=F)
  gdp.mod <- gdp(y=y, x=x, B=2000, burn=500, printProg=F,returnHyper=F)
  mod <- list("lasso_mod"=lasso.mod, "ridge_mod"=ridge.mod, "spike_mod"=spike.mod, "blasso_mod"=blasso.mod$beta, "gdp_mod"=gdp.mod$beta, "param_index"=param_index)
  mod
}

mod.params.ind <- as.list(1:numdat); counter <- 0
for (n_i in 1:length(n)) for (p_i in 1:length(p)) for (S_i in 1:length(Sig_Makers)) for (beta_i in 1:length(betas)) {
  counter <- counter + 1
  mod.params.ind[[counter]] <- list("n"=n_i,"p"=p_i,"S"=S_i,"beta"=beta_i)
}

mod <- foreach(mpi=mod.params.ind) %dopar% oneSim(mpi$n,mpi$p,mpi$S,mpi$beta)

#mod_n <- 35 
#plot( beta_lookup(mod[[mod_n]]$param_index), col="grey", pch=1, ylim=c(-5,5), bty="n", cex =2, xlab="Parameter Index", ylab=expression(beta));
#points( apply(mod[[mod_n]]$spike,2,mean),col='red',pch=20 );
#points( apply(mod[[mod_n]]$gdp,2,mean),col='green',pch=20 );
#points( apply(mod[[mod_n]]$blasso,2,mean),col='blue',pch=20 );
#
#points( coef(mod[[mod_n]]$lasso),col='blue',pch=20 );
#points( coef(mod[[mod_n]]$ridge),col='red',pch=20 );


# - For all theses Bayesian models, get E(\beta_j | y)
rmse <- function(b,b_true) sqrt(mean((b-b_true)^2))
compareBayesian <- function(model,rmse_ord="blasso",ylim_rmse=c(0,5),cex_rmse=1,pos_rmse=0,cex.axis_rmse=.6) {
  cc <- 0
  rmse_blasso <- NULL; rmse_spike <- NULL; rmse_gdp <- NULL;
  for (mm in model) {
    cc <- cc + 1
    mod_ind <- mm$param_index
    beta_true <- beta_lookup(mod_ind)
    rmse_blasso[cc] <- rmse(apply(mm$blasso,2,mean),beta_true)
    rmse_spike[cc] <- rmse(apply(mm$spike$beta,2,mean),beta_true)
    rmse_gdp[cc] <- rmse(apply(mm$gdp,2,mean),beta_true)
  }

  ord <- 1:length(model)
  if (rmse_ord == "blasso") {
    ord <- order(rmse_blasso)
  } else if (rmse_ord == "gdp") {
    ord <- order(rmse_gdp)
  } else if (rmse_ord == "spike") {
    ord <- order(rmse_spike)
  }

  par(mar=c(4,4,0,0))
  plot(rmse_blasso[ord],col="black",pch=1,ylim=ylim_rmse,cex=cex_rmse,
       ylab="RMSE",xlab="",lwd=2,xaxt="n",bty="n")
  points(rmse_spike[ord],col="black",pch=2,cex=cex_rmse,lwd=2)
  points(rmse_gdp[ord],col="black",pch=4,cex=cex_rmse,lwd=2)

  simdat <- t(sapply(ord, function(o) mod[[o]]$param_index))
  simdat <- cbind(
    ifelse(simdat[,1] == 1, "n50", "n500"),
    ifelse(simdat[,2] == 1, "p100", "p1000"),
    ifelse(simdat[,3] == 1, "I", ifelse(simdat[,3] == 2, "S.1", "S.6")),
    ifelse(simdat[,4] == 1, "b1", ifelse(simdat[,4] == 2, "b2", "b3")))
  lab <- apply(simdat,1,function(x) paste(x,collapse="\n"))
  axis(1,at=1:numdat,label=lab,pos=pos_rmse,tck=0,lty=0,cex.axis=cex.axis_rmse)
  legend("topleft",cex=2,legend=c("blasso","ssvn","gdp"),pch=c(1,2,4),bty="n")
  par(mar=c(5.1,4.1,4.1,2.1))
  abline(v=5*(1:7),col=rgb(.5,.5,.5))

  list("ord"=ord, "blasso"=rmse_blasso, "gdp"=rmse_gdp, "spike"=rmse_spike)
}

# - discuss accuracy wrt to a metric
pdf("output/rmseblasso.pdf",w=16,h=9)
  rmse_blasso <- compareBayesian(mod,"blasso",ylim=c(0,3),pos_rmse=-.2,cex_rmse=1.5)
dev.off()
pdf("output/rmsessvn.pdf",w=16,h=9)
  rmse_spike  <- compareBayesian(mod,"spike",ylim=c(0,3),pos_r=-.2,cex_rmse=1.5)
dev.off()
pdf("output/rmsegdp.pdf",w=16,h=9)
  rmse_gdp    <- compareBayesian(mod,"gdp",ylim=c(0,3),pos_r=-.2,cex_rmse=1.5)
dev.off()

sapply(rmse_gdp[-1],mean) # blasso seems to perform a little better than gdp and a lot better than spike and slab
#      ord    blasso       gdp     spike
#18.500000     .3141     .3722     .5840

# Compare Lasso and SSVN

cms <- matrix(0,numdat,4)
colnames(cms) <- c("++lasso","++ssvn","--lasso","--ssvn")
for (mod_n in 1:numdat) {
  true_beta <- beta_lookup(mod[[mod_n]]$param_index)
  postmean.gamma <- apply(mod[[mod_n]]$spike$gam,2,mean)
  sel_lasso <- !!as.numeric(coef(mod[[mod_n]]$lasso) != 0)
  sel_ssvn <- postmean.gamma > .5
  cms[mod_n,1] <- sum(true_beta !=0 & sel_lasso) / sum(true_beta != 0)
  cms[mod_n,2] <- sum(true_beta !=0 & sel_ssvn) / sum(true_beta !=0)
  cms[mod_n,3] <- sum(true_beta == 0 & !sel_lasso) / sum(true_beta == 0)
  cms[mod_n,4] <- sum(true_beta == 0 & !sel_ssvn) / sum(true_beta == 0)
  print(mod_n)
}

compareVS <- function(c1=1,c2=2,ord=1,pos_rmse=0,cex.axis_rmse=1) {
  ordd <- order(cms[,ord])
  simdat <- t(sapply(ordd, function(o) mod[[o]]$param_index))
  simdat <- cbind(
    ifelse(simdat[,1] == 1, "n50", "n500"),
    ifelse(simdat[,2] == 1, "p100", "p1000"),
    ifelse(simdat[,3] == 1, "I", ifelse(simdat[,3] == 2, "S.1", "S.6")),
    ifelse(simdat[,4] == 1, "b1", ifelse(simdat[,4] == 2, "b2", "b3")))

  par(mar=c(4,4,0,0))
  plot(cms[ordd,c1],ylim=0:1,pch=20,lwd=2,xlab="",xaxt="n",ylab="proportion")
  points(cms[ordd,c2],pch=2,lwd=2)
  lab <- apply(simdat,1,function(x) paste(x,collapse="\n"))
  axis(1,at=1:numdat,label=lab,pos=pos_rmse,tck=0,lty=0,cex.axis=cex.axis_rmse)
  legend("right",cex=2,legend=c(colnames(cms)[c1],colnames(cms)[c2]),pch=c(20,2),bty="n")
  par(mar=c(5.1,4.1,4.1,2.1))
  abline(v=5*(1:7),col=rgb(.5,.5,.5))
}

#compareVS(1,2,1,-.05,.6)
pdf("output/lassoVssvnPP.pdf",w=16,h=9); compareVS(1,2,1,-.05,.6); dev.off() 
pdf("output/ssvnVlassoPP.pdf",w=16,h=9); compareVS(1,2,2,-.05,.6); dev.off()
pdf("output/lassoVssvnFF.pdf",w=16,h=9); compareVS(3,4,3,-.05,.6); dev.off()
pdf("output/ssvnVlassoFF.pdf",w=16,h=9); compareVS(3,4,4,-.05,.6); dev.off()

apply(cms,2,mean)
#  ++lasso    ++ssvn   --lasso    --ssvn
#0.7220000 0.7208333 0.6391333 0.8962522

# - Compute mean(L_j: beta_j == 0)
m_zero <- function(b_post,b_true,zero=T) {
  ind_zero <- which(b_true == 0)
  out <- NULL
  if (zero) 
    out <- mean(apply(as.matrix(b_post[,ind_zero]),2, function(x) diff(quantile(x,c(.025,.975)))))
  else 
    out <- mean(apply(as.matrix(b_post[,-ind_zero]),2, function(x) diff(quantile(x,c(.025,.975)))))
  out
}

compareL <- function(zero=T,L_ord="blasso",ylim_L=c(0,5),cex_L=1,pos_L=0,cex.axis_L=.6) {
  cc <- 0
  L_blasso <- NULL; L_spike <- NULL; L_gdp <- NULL;
  for (mm in mod) {
    cc <- cc + 1
    mod_ind <- mm$param_index
    beta_true <- beta_lookup(mod_ind)
    L_blasso[cc] <- m_zero(mod[[cc]]$blasso, beta_true,zero)
    L_spike[cc] <- m_zero(mod[[cc]]$spike$b, beta_true,zero)
    L_gdp[cc] <- m_zero(mod[[cc]]$gdp, beta_true,zero)
    print(cc)
  }

  ord <- 1:length(mod)
  if (L_ord == "blasso") {
    ord <- order(L_blasso)
  } else if (L_ord == "gdp") {
    ord <- order(L_gdp)
  } else if (L_ord == "ssvn") {
    ord <- order(L_spike)
  }

  par(mar=c(4,4,0,0))
  plot(L_blasso[ord],col="black",pch=1,ylim=ylim_L,cex=cex_L,
       ylab=ifelse(zero,"M_zero","M_nonzero"),xlab="",lwd=2,xaxt="n",bty="n")
  points(L_spike[ord],col="black",pch=2,cex=cex_L,lwd=2)
  points(L_gdp[ord],col="black",pch=4,cex=cex_L,lwd=2)

  simdat <- t(sapply(ord, function(o) mod[[o]]$param_index))
  simdat <- cbind(
    ifelse(simdat[,1] == 1, "n50", "n500"),
    ifelse(simdat[,2] == 1, "p100", "p1000"),
    ifelse(simdat[,3] == 1, "I", ifelse(simdat[,3] == 2, "S.1", "S.6")),
    ifelse(simdat[,4] == 1, "b1", ifelse(simdat[,4] == 2, "b2", "b3")))
  lab <- apply(simdat,1,function(x) paste(x,collapse="\n"))
  axis(1,at=1:numdat,label=lab,pos=pos_L,tck=0,lty=0,cex.axis=cex.axis_L)
  legend("topleft",cex=2,legend=c("blasso","ssvn","gdp"),pch=c(1,2,4),bty="n")
  par(mar=c(5.1,4.1,4.1,2.1))
  abline(v=5*(1:7),col=rgb(.5,.5,.5))

  list("ord"=ord, "blasso"=L_blasso, "gdp"=L_gdp, "ssvn"=L_spike)
}

pdf("output/Lblasso.pdf",w=16,h=9)
  L_blasso <- compareL(zero=T,L_ord="blasso",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
pdf("output/Lssvn.pdf",w=16,h=9)
  L_ssvn <- compareL(zero=T,L_ord="ssvn",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
pdf("output/Lgdp.pdf",w=16,h=9)
  L_gdp <- compareL(zero=T,L_ord="gdp",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
#L_blasso
#compareL(zero=T,L_ord="ssvn",ylim_L=c(0,50),cex_L=1,pos_L=-.3,cex.axis_L=.6)

pdf("output/L1blasso.pdf",w=16,h=9)
  L1_blasso <- compareL(zero=F,L_ord="blasso",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
pdf("output/L1ssvn.pdf",w=16,h=9)
  L1_ssvn <- compareL(zero=F,L_ord="ssvn",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
pdf("output/L1gdp.pdf",w=16,h=9)
  L1_gdp <- compareL(zero=F,L_ord="gdp",ylim_L=c(0,5),cex_L=1,pos_L=-.3,cex.axis_L=.6)
dev.off()
#L1_blasso


# - posterior pred.  mod[[7]]

sim.post.pred <- function(B=2000,burn=1000) {
  x <- t(sapply( 1:(50), function(x) c(1, mvrnorm(0,Sig_Maker(.6,100) ))))
  y <- x %*% betas[[1]](100) + rnorm(50)

  tt <- .05^2
  cat("\nssvs: \n")
  ssvs_mod <- spikeAndSlab(y=y, x=x, tau2=rep(tt,ncol(x)), g=100/tt, w=rep(.5,ncol(x)), B=B+burn, burn=burn, printProg=T, returnHyper=T)
  cat("\nblasso: \n")
  blasso_mod <- blasso(y=y,x=x,r=1,delta=1.5,B=B+burn,burn=burn, printProg=T, returnHyper=F)
  cat("\ngdp: \n")
  gdp_mod <- gdp(y=y, x=x, B=B+burn, burn=burn, printProg=T,returnHyper=F)

  get.post.pred <- function(b,g=rep(1,ncol(x))) {
    G <- diag(g)
    y.pred <- x %*% G %*% t(b)
    apply(t(y.pred),2,mean)
  }

  blasso.pp <- get.post.pred(blasso_mod$beta)
  gdp.pp <- get.post.pred(gdp_mod$beta)
  print(1)
  ssvn.pp <- get.post.pred(ssvs_mod$beta, ifelse(apply(ssvs_mod$gam,2,mean) > .5, 1,0))

  list("y"=y,"x"=x,"blasso"=blasso.pp,"gdp"=gdp.pp,"ssvs"=ssvs.pp)
}

sim.pp <- sim.post.pred(B=100,burn=10)
