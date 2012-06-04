### functions to estimate and shrink gene specific over dispersion
## this works for single factor design only at this time!!!
estDispersion <- function(seqData) {
  if(!is(seqData, "SeqCountSet"))
    stop("Input must be an object of SeqCountSet class!")

  design=as.factor(pData(phenoData(seqData))$designs)
  k=normalizationFactor(seqData)

  ## make design matrix
  X=makeDesign(design, k)
  
  ## compute the expected Y
  Y=exprs(seqData)+0.5
  muY = calc.expY(Y, X)
  
  ## estimate gene specific phi
  df=ncol(Y)-ncol(X)
  phi.g = est.dispersion(Y, k, design)
  
  ## shrink phi
  phi.hat = shrink.dispersion(phi.g, Y, muY, k, design)
  dispersion(seqData) = phi.hat

  seqData
}

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

## Estimate gene specific over dispersion using MOM
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


### function to compute expected value for Y
calc.expY <- function(Y, X) {
  b=(Y %*% X) %*% solve(t(X)%*%X)
  muY=b %*% t(X)
  muY
}
  
### function to shrink dispersion
## Y is normalized by size factors and designs. 
shrink.dispersion <- function(phi.g, Y, muY, k, design) {
  nsamples=ncol(Y)
  ngenes=nrow(Y)
  phi.hat=rep(0,ngenes)
  
  ## estimate hyper parameters - this is tricky!!
  s=rowMeans(Y)
  ii=s>=5
  phi.g0=phi.g[ii]
  lphi.g0=log(phi.g0)
  ## It seems phi.g is a little under estimated. 
  mu=median(lphi.g0, na.rm=TRUE) ## it seems mu is under estimated when var(OD) is small
  sigma2.mar=(IQR(lphi.g0, na.rm=TRUE) / 1.349)^2
  ## estimate the part being over estimated and subtract.
  ## It seems I'm subtracting too much when OD is big!!
  sigma2.base=compute.baseSigma(exp(mu), Y[ii,], muY[ii,], k, design)
  sigma=sqrt(max(sigma2.mar-sigma2.base, 1e-2))
  ## The way to estimate mu and sigma needs some thinking.
  ## However it doesn't seem to make much differences. 
  ## NR procedure to do shrinkage
  max.value = max(lphi.g0,10, na.rm=TRUE)
  ## objective function (penalized likelihood)
  get.phi <- function(dat){
    y=dat[1:nsamples]
    Ey=dat[nsamples+(1:nsamples)]
    obj=function(phi) {
      alpha=1/phi
      tmp1=1/(1+Ey*phi)
      tmp2=1-tmp1
      -(sum(lgamma(alpha+y)) - nsamples*lgamma(alpha) + alpha*sum(log(tmp1)) + sum(y*log(tmp2)) -
        ((log(phi) - mu)^2) / (2*(sigma^2)) + log(phi) - log(sigma))
    }
    return(optimize(obj, interval=c(0.01, max.value))$minimum)
  }
  tmp=cbind(Y, muY)
  phi.hat=apply(tmp,1,get.phi)
  phi.hat
}

### compute hyperparameter sigma, when all genes have the same phi.
## this seems over estimated that. Need to think!!!
compute.baseSigma <- function(phi0, Y, muY,k, design) {
  n=length(phi0)
  m=ncol(Y)
  nsim=10
  res=rep(0, nsim)
  ## sample
  for(i in 1:nsim) {
    Ysim <- matrix(rnegbinom(length(Y), muY, phi=phi0), ncol=m)
    lag.sim=log(est.dispersion(Ysim, k, design))
    sigma=IQR(lag.sim, na.rm=TRUE) / 1.349
    res[i]=sigma^2
  }
  mean(res)
}


  

