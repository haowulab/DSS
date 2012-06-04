### some utility functions
## variance of each row
rowVars <- function(x) {
  n0 <- ncol(x)
  EX <- rowMeans(x, na.rm=TRUE)
  EX2 <- rowMeans(x^2, na.rm=TRUE)
      vv <- (EX2-EX^2) * n0 / (n0-1)
  vv
}

### generate negative binomial rv, given mu and phi (over-dispersion)
rnegbinom <- function (n, mu =1, phi=0.01){
  rpois(n, rgamma(n, shape=1/phi,scale=mu*phi))
}

