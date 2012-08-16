########################################################################
## some utility functions for estimating and shrinking dispersions
########################################################################

### function to create design matrix based on design and size factor.
## Works for single factor now.
makeDesign <- function(design, k) {
  allfactors <- levels(design)
  X=matrix(0, nrow=length(k), ncol=length(allfactors))
  for(i in 1:length(allfactors)) {
    idx=design==allfactors[i]
    X[idx,i]=1
  }
  if(is.matrix(k))
    stop("Matrix normalization factor is still under development!")
  else
    X=sweep(X, 1, k, FUN="*")
  X
}

## Estimate gene specific over dispersion using MOM.
## This considers the experimental design
est.dispersion <- function(Y, k, design) {
  ## adjust by size factor
  if(is.matrix(k))
    Y2=Y/k
  else
    Y2=sweep(Y, 2, k, FUN="/")
  
  ## transform
  z=Y^2-Y
  z=sweep(z, 2, k^2, FUN="/")

  ## first find group labels and compute group means.  works only for two groups now
  ll=levels(design)
  idxA=design==ll[1]; muA=rowMeans(Y2[,idxA])
  idxB=design==ll[2]; muB=rowMeans(Y2[,idxB])
  mus=matrix(0, nrow=nrow(Y), ncol=ncol(Y))
  mus[,idxA]=muA; mus[,idxB]=muB
  phi=rowSums(z) / rowSums(mus^2) - 1
##   phiA=rowMeans(z[,idxA])/(muA^2) - 1
##   phiB=rowMeans(z[,idxB])/(muB^2) - 1
##   phi=(phiA+phiB)/2
  phi[phi<1e-8]=NA
  phi
}

## Estimate gene specific over dispersion using MOM.
## This assumes all data are from the same distribution
est.phi0 <- function(Y) {
  k=colSums(Y) ## naive estimation of size factor
  k=k/min(k)
  Y2=sweep(Y, 2, k, FUN="/")
  mm0=rowMeans(Y2)
  idx0=mm0==0
  res=rep(0, nrow(Y))
  res[idx0]=NA
  
  ix=mm0>0
  Y2=Y2[ix,]
  m=rowMeans(Y2)
  v=apply(Y2, 1, var)
  phi.g = (v-m)/m^2
  res[!idx0]=phi.g
  res[res<1e-5]=1e-5
  
  ##n0=rowMeans(Y==0)
  ##ix.valid=mm0>10 & n0<=0.2
  ##list(phi=res, ix.valid=ix.valid)
  res
}


### function to compute expected value for Y
calc.expY <- function(Y, X) {
  b=(Y %*% X) %*% solve(t(X)%*%X)
  muY=b %*% t(X)
  muY
}
