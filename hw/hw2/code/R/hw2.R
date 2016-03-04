set.seed(268)

#library(doMC)
#registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../../../../cpp_functions/func.cpp")
sourceCpp("../C++/albertChib.cpp")
sourceCpp("../C++/smc.cpp")
system("mkdir -p output")
"%bw%" <- function(x,rng) x<rng[2] & x>rng[1]

n <- 1000
p <- 10#400

system.time( X <- t(sapply(1:n,function(x) mvrnorm(rep(0,p),diag(p))))) #54s
beta <- c(1.3, 4, -1, 1.6, 5, -2, rep(0,p-6))
y <- rbinom(n,1,pnorm(X%*%beta))

# Albert and Chib Probit Gibbs
b.post  <- gibbs_ac(y,X,diag(ncol(X)),B=2000,burn=1000,printProg=T);

# Plotting Functions
add.errbar <- function(ci,...) {
  x <- 1:nrow(ci)
  segments(x,ci[,1],x,ci[,2],...)
}
plot.coverage <- function(b_post,b_true=beta,...) {
  ci <- t(apply(b_post,2,function(x) quantile(x,c(.025,.975))))
  M <- cbind(b_true[1:6],ci[1:6,]); M <- cbind(M,t(apply(M,1,function(x) c( diff(c(x[2],x[3])), x[1] %bw% c(x[2],x[3])) )))
  colnames(M) <- c("Truth","CI lower","CI upper","Length","Coverage")
  rng <- range(M[,1:3])

  plot(M[,1],ylim=rng,bty="n",xaxt="n",...)
  axis(1,at=1:6,label=paste0("b",1:6))
  add.errbar(ci,lwd=4,col="blue")
}

# Tables & Plots
ci <- t(apply(b.post,2,function(x) quantile(x,c(.025,.975))))
M <- cbind(beta[1:6],ci[1:6,]); M <- cbind(M,t(apply(M,1,function(x) c( diff(c(x[2],x[3])), x[1] %bw% c(x[2],x[3])) )))
colnames(M) <- c("Truth","CI lower","CI upper","Length","Coverage"); M
plot.coverage(b.post,pch=20,col="grey",cex=2,ylab="",xlab=expression(beta))



# SMC idk...
sourceCpp("../C++/smc.cpp")
N <- 500
out_2 <- update_smc(y=y, x=X, N=N, beta=matrix(rnorm(N*p),N,p), printProg=TRUE)
plot(beta,pch=20,col="grey",cex=3)
lines(apply(out_2,2,mean),col='blue',lwd=3)
#update_smc(y, X[,1:3], diag(3), 10, 4, F)
#update_smc(y, X[,1:3], diag(3), 10, 4, F)

