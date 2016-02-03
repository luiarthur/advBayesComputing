library(MASS) # ridge.lm

library(doMC)
registerDoMC(ncore <- as.numeric(system("nproc",intern=T)))

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
betas <- list( function(p) {b <- rep(0,p); b[1:5] <- 3; b <- c(1,b); b }, 
               function(p) {b <- rep(0,p); b[1:5] <- 5; b[6:10] <- -2; b[11:15] <- .5; b <- c(1,b); b },
               function(p) {b <- rep(1,p); b <- c(1,b); b })

numdat <- length(n) * length(p) * length(Sig_Makers) * length(betas)
counter <- 0
dat <- as.list(1:numdat)


for (nn in n) {
  for (pp in p) {
    for (ss in Sig_Makers) {
      for (bb in betas) {
        counter <- counter + 1

        # 1) Simulate the data
        Sig <- ss(pp)
        x <- foreach(i=1:nn,.combine=rbind) %dopar% c( 1 ,c(mvrnorm(0,Sig)) )
        beta <- bb(pp)
        y <- x %*% beta + rnorm(nn)

        # 2a)
        # Run Lasso

        # Run Ridge
        mod <- lm.ridge(y~x[,-1])

        # 2b)
        # Run Bayesian spike and slab (S&S)
        # Run Bayesian Lasso
        # Run Bayesian Generalized Double Pareto
        # - For all theses Bayesian models, get E(\beta_j | y)
        # - discuss accuracy wrt to a metric
        # - Compare Lasso and S&S for variable selection
        # - Compute mean(L_j: beta_j == 0)
        # - posterior pred.


        #dat[[counter]] <- list("y"=y, "x"=x, "S"=s, "n"=nn, "p"=pp)

        cat("\r n: ", nn, ";  p: ", pp, "; counter: ", counter)
      }
    }
  }
}


