### dispersion shrinkage with multiple factor design
## Assuming dispersion is not dependent on means at this time
estDispersion.multiFactor <- function(seqData) {
  Y = exprs(seqData)
  k = normalizationFactor(seqData)
  
  ## design matrix
  design = as.matrix(pData(seqData))
  ## incorporate size factors
  Y2 = sweep(Y, 2, k, FUN="/")

  ## 1. Run olr, using log(counts) as outcomes, obtain E[Y]
  Y2=log(Y2+0.5)
  muY = exp(calc.expY(Y2, design))
  
  ## 2.  estimate gene specific phi
  z=Y^2-Y
  z=sweep(z, 2, k^2, FUN="/")
  phi.g=rowSums(z)/rowSums(muY^2)-1
  phi.g[phi.g<1e-8]=NA

  ## 3. shrinkage 
  ## first estimate hyper parameters
  s=rowMeans(Y)
  ii=s>=5
  phi.g0=phi.g[ii]
  lphi.g0=log(phi.g0)
  mu=median(lphi.g0, na.rm=TRUE) ## it seems mu is under estimated when var(OD) is small
  tau=sqrt((IQR(lphi.g0, na.rm=TRUE) / 1.349)^2)
  ## scale tau down a little bit. This needs more investigation,
  ## but it seems this doesn't affect the results much.
  tau = tau * 0.7
  
  ## NR procedure to do shrinkage
  max.value = max(lphi.g0,10, na.rm=TRUE)
  ## objective function (penalized likelihood). Shrink in the log scale to geometric means
  nsamples = nrow(design)
  get.phi <- function(dat){
    y=dat[1:nsamples]
    Ey=dat[nsamples+(1:nsamples)]
    
    obj=function(phi) {
      alpha=1/phi
      tmp1=1/(1+Ey*phi)
      tmp2=1-tmp1
      -(sum(lgamma(alpha+y)) - nsamples*lgamma(alpha) + alpha*sum(log(tmp1)) + sum(y*log(tmp2)) +
        dnorm(log(phi), mean=mu, sd=tau, log=TRUE))
    }
    return(optimize(obj, interval=c(0.01, max.value))$minimum)
  }
  tmp=cbind(Y, muY)
  phi.hat=apply(tmp,1,get.phi)

  dispersion(seqData) = phi.hat
  seqData
}
