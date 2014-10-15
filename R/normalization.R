### function to compute normalization factors for an object of SeqCountSet
### I'm using a naive method now.
estNormFactors <- function(seqData, method=c("lr", "total", "quantile", "median")) {
  if(!is(seqData, "SeqCountSet"))
    stop("Input must be an object of SeqCountSet class!")

  method=match.arg(method)
  X <- exprs(seqData)
  if(method=="lr")
      k=estNormFactors.lr(X)
  if(method=="quantile")
      k=estNormFactors.quantile(X)
  else if(method=="total")
      k=estNormFactors.total(X)
  else if(method=="median")
      k=estNormFactors.median(X)

  normalizationFactor(seqData)=k
  seqData
}

estNormFactors.quantile <- function(X) {
  k=apply(X, 2, function(x) quantile(x[x>0], 0.75))
  k/min(k)
}

estNormFactors.total <- function(X) {
  k=colSums(X)
  k/min(k)
}

estNormFactors.median <- function(X) {
  k=apply(X, 2, median)
  k/min(k)
}

estNormFactors.lr <- function(X) {
    ss <- colSums(X)
    ix <- which.min(ss)## use this one as basis
    X0 <- X[,ix]
    k <- rep(1, ncol(X))

    for(i in 1:ncol(X)) {
        ii <- X0>0 & X[,i]>0
        lr <- log(X0[ii] / X[ii,i])
        ## trim 5%
        qq <- quantile(lr, c(0.05, 0.95))
        lr2 <- lr[lr>qq[1] & lr<qq[2]]
        k[i] <- exp(-median(lr2))
    }
    k[ix] <- 1
    k
}
