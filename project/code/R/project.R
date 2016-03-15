set.seed(268)
library(Rcpp)
system("mkdir -p output")
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
source("sim_dat.R")

#cat("sourcing func.cpp...\n"); sourceCpp("../../../cpp_functions/func.cpp")
#cat("sourcing gp.cpp...\n"); sourceCpp("../C++/gp.cpp")
cat("sourcing gp_gdp.cpp...\n"); sourceCpp("../C++/gp_gdp.cpp")
cat("Starting Main Program...\n")


### p = 3
#f3 <- function(X) apply(as.matrix(X),1, function(x) x[1] + ifelse(x[2]>.5,x[2]-.5,0) + x[3]^2)
#dat <- sim_dat(f3,p=3,n=100)
#priors <- c(5,5,  .1,5,  .1,.1) #s2, phi, tau
#out <- gp(y=scale(dat$y), x=dat$x, D=dat$D, cand_S=diag(dat$p)*.1,
#          init=rep(1,3), priors=priors, B=2000, burn=4000, printProg=TRUE)
#colnames(out$param) <- c("s2","phi","tau")
#plot.posts(out$param,cex.l=1.3,cex.a=1,names=colnames(out$param))
#plot(ts(out$param))
#apply(out$param,2,summary)


cat("sourcing gp_gdp.cpp...\n"); sourceCpp("../C++/gp_gdp.cpp")
priors <- c(2,1,    .1,5,    3,5,    1,1) #s2, phi, tau, d_vec

f1 <- function(X) apply(as.matrix(X),1, function(x) -.2*x[1] + sin(x[2]) + .3*x[3]^2)
f2 <- function(X) apply(as.matrix(X),1, function(x) x[1] + sin(5*x[2]) + sin(x[3]))
f3 <- function(X) apply(as.matrix(X),1, function(x) 3*x[1] + -.3*x[2] + x[3])

fn <- list(f1,f2,f3)
dat <- lapply(fn,function(f) sim_dat(f,p=10,n=100))

one_sim <- function(d,B=2000,burn=2000) {
  gp_gdp(y=scale(d$y), X=d$x, cand_S=diag(c(rep(.1,3),rep(1e-3,d$p))),
         init=c(rep(0,3),rep(1,d$p)), priors=priors, B=B, burn=burn, printProg=TRUE)
}

plot_d <- function(o,...) {
  d_vec <- t(apply(o$param[,-c(1:3)],1,function(x) x / sum(x)))
  plot( apply(d_vec,2,mean) ,col=c(rep("red",3),rep("grey",dat$p)),pch=20,cex=5,ylab="d",...)
  add.errbar(t(apply(d_vec,2,function(x) quantile(x,c(.025,.975)))),lwd=4,col="black");abline(h=0,lty=2)
}

out <- as.list(1:3)
out[[1]] <- one_sim(dat[[1]],2000,4000)
out[[2]] <- one_sim(dat[[2]],2000,4000)
out[[3]] <- one_sim(dat[[3]],2000,4000)

for (i in 1:3) {
  pdf(paste0("output/post",i,".pdf"),w=13,h=9)
    plot.posts(out[[i]]$param[,1:3],cex.l=.01,cex.a=1.7,names=c("s2","phi","tau"),tck.dig=3)
  dev.off()
  pdf(paste0("output/trace",i,".pdf"),w=13,h=9)
    plot(ts(out[[i]]$param[,-c(1:3)]),main="")
  dev.off()
  pdf(paste0("output/d",i,".pdf"),w=13,h=9)
    plot_d(out[[i]],ylim=c(-1,1),cex.axis=2)
  dev.off()
}
pdf(paste0("output/d",2,".pdf"),w=13,h=9)
  plot_d(out[[2]],ylim=c(-100,100),cex.axis=2)
dev.off()

