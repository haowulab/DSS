### function to compute normalization factors for an object of SeqCountSet
### I'm using a naive method now.
estNormFactors <- function(seqData, method=c("quantile", "total", "median")) {
  if(!is(seqData, "SeqCountSet"))
    stop("Input must be an object of SeqCountSet class!")
  
  method=match.arg(method)
  X <- exprs(seqData)
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
  k=apply(X, 2, quantile, 0.75)
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

