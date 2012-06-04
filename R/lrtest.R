## likihood ratio test.
## This works for two group comparison only!!
lrTest <- function(seqData) {
  ## compute EY under null
  Y=exprs(seqData)
  k=normalizationFactor(seqData)
  if(is.matrix(k))
    Y2=Y/k
  else
    Y2=sweep(Y, 2, k, FUN="/")
  muY0=sweep(matrix(rep(rowMeans(Y2),ncol(Y)),ncol=ncol(Y)),  2, k, FUN="*")
  
  ## EY under alternative
  design=as.factor(pData(phenoData(seqData))$designs)
  k=normalizationFactor(seqData)
  X=makeDesign(design, k)
  muY1 = calc.expY(Y+0.5, X)

  ## likelihoods
  phi=dispersion(seqData)
  l0=rowSums(dnbinom(Y, size=1/phi, mu=muY0, log=TRUE))
  l1=rowSums(dnbinom(Y, size=1/phi, mu=muY1, log=TRUE))
  lr=l1-l0 ## some are less than 0??
  
  pval=1-pchisq(lr, df=1)
  list(lr=lr, pval=pval)
}


                    
