sim_dat <- function(f,p,n,m=round(sqrt(n)),s2=.5,v=1) { # function, num of covariates, num of observations
  x <- matrix(rnorm(n*p,0,v),n,p)
  y <- f(x) + rnorm(n,0,sqrt(s2))

  s <- matrix(rnorm(m*p),m,p)

  list("y"=y,"f"=f,"x"=x,"D"=as.matrix(dist(x)),"s"=s,"p"=p,"s2"=s2)
}

# Plotting Functions
add.errbar <- function(ci,...) {
  x <- 1:nrow(ci)
  segments(x,ci[,1],x,ci[,2],...)
}

