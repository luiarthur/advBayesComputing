set.seed(268)

#library(doMC)
#registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
sourceCpp("../../../../cpp_functions/func.cpp")
sourceCpp("../C++/albertChib.cpp")
system("mkdir -p output")
"%bw%" <- function(x,rng) x<rng[2] & x>rng[1]

n <- 1000
p <- 400

system.time( X <- t(sapply(1:n,function(x) mvrnorm(rep(0,p),diag(p))))) #12s
beta <- c(1.3, 4, -1, 1.6, 5, -2, rep(0,p-6))
y <- rbinom(n,1,pnorm(X%*%beta))

b.post  <- gibbs_ac(y,X,diag(ncol(X)),B=2000,burn=1000,printProg=T);
b.post.mean <- apply(b.post,2,mean)
#plot(b.post[,1])
#plot(b.post[,7])

plot(b.post.mean,pch=20,col='blue')
points(beta,pch=20,col='grey')

ci <- t(apply(b.post,2,function(x) quantile(x,c(.025,.975))))
ind <- 7; beta[ind] %bw% ci[ind,] # used to determine coverage
diff(ci[ind,]) # length of the ci
