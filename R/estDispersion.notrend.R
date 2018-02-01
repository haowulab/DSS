##### Functions to estimate and shrink dispersions
## when there's no trend between dispersion and expression
estDispersion.notrend <- function(seqData) {
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
  phi.hat = shrink.dispersion.notrend(phi.g, Y, muY, k, design)
  dispersion(seqData) = phi.hat

  seqData
}

### function to shrink dispersion
## Y is normalized by size factors and designs.
shrink.dispersion.notrend <- function(phi.g, Y, muY, k, design) {
  nsamples=ncol(Y)
  ngenes=nrow(Y)
  phi.hat=rep(0,ngenes)

  ## estimate hyper parameters - this is tricky!!
  s=rowMeans(Y)
  ii=s>=5
  phi.g0=phi.g[ii]
  lphi.g0=log(phi.g0)
  ## It seems phi.g is a little under estimated.
  m0=median(lphi.g0, na.rm=TRUE) ## it seems m0 is under estimated when var(OD) is small
  sigma2.mar=(IQR(lphi.g0, na.rm=TRUE) / 1.349)^2
  ## estimate the part being over estimated and subtract.
  ## It seems I'm subtracting too much when OD is big!!
  sigma2.base=compute.baseSigma.notrend(exp(m0), Y[ii,], muY[ii,], k, design)
  sigma=sqrt(max(sigma2.mar-sigma2.base, 1e-2))

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
        ((log(phi) - m0)^2) / (2*(sigma^2)) - log(phi) - log(sigma))
    }
    return(optimize(obj, interval=c(0.01, max.value))$minimum)
  }
  tmp=cbind(Y, muY)
  phi.hat=apply(tmp,1,get.phi)
  phi.hat
}

### compute hyperparameter sigma, when all genes have the same phi.
## this seems over estimated that. Need to think!!!
compute.baseSigma.notrend <- function(phi0, Y, muY,k, design) {
  n=length(phi0)
  m=ncol(Y)
  nsim=10
  res=rep(0, nsim)
  ## sample, ignore size factor. This is okay.
  k=rep(1, m)
  for(i in 1:nsim) {
    Ysim <- matrix(rnegbinom(length(Y), muY, phi=phi0), ncol=m)
    lag.sim=log(est.dispersion(Ysim, k, design))
    sigma=IQR(lag.sim, na.rm=TRUE) / 1.349
    res[i]=sigma^2
  }
  mean(res)
}




