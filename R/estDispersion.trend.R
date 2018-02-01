##### Functions to estimate and shrink dispersions
## when there IS trend between dispersion and expression
estDispersion.trend <- function(seqData) {
  design=as.factor(pData(phenoData(seqData))$designs)
  k=normalizationFactor(seqData)

  ## make design matrix
  X=makeDesign(design, k)

  ## compute the expected Y. Should this be under null???
  Y=exprs(seqData)+0.5
  muY = calc.expY(Y, X)

  ## estimate gene specific phi
  df=ncol(Y)-ncol(X)
  phi.g = est.dispersion(Y, k, design)

  ## estimate phi~expr trends
##   lexpr1=rowMeans(log(sweep(Y+0.5,2,k,FUN="/")))
##   lOD=log(phi.g)
##   ix=!is.na(lOD)
##   fit1=smooth.spline(lOD[ix]~lexpr1[ix], df=3, nknots=round(sum(ix)/50))
##   lOD.pred=predict(fit1,lexpr1)$y
  lOD.pred=est.trend(Y)
  ## shrink phi
  phi.hat = shrink.dispersion.trend(phi.g, lOD.pred, Y, muY, k, design)
  dispersion(seqData) = phi.hat

  seqData
}


### function to shrink dispersion
## Y is normalized by size factors and designs.
## lOD.pred is predicted dispersion in log scale.
shrink.dispersion.trend <- function(phi.g, lOD.pred, Y, muY, k, design) {
  nsamples=ncol(Y)
  ngenes=nrow(Y)
  phi.hat=rep(0,ngenes)

  ## estimate hyper parameters - this is tricky!!
  ## lOD.pred is the mean
  ## need to estimate tau^2
  lphi.g0=log(phi.g)
  lexpr1=rowMeans(log(sweep(Y,2,k,FUN="/")))
  lexpr.cut=2
  ix=lexpr1>lexpr.cut
  sigma2.mar=(IQR(lphi.g0[ix]-lOD.pred[ix], na.rm=TRUE) / 1.349)^2

  ## remove the amount of over-estimation
  sigma2.base=compute.baseSigma.trend(lOD.pred, Y, muY, k, design)
  sigma=sqrt(max(sigma2.mar-sigma2.base, 1e-2))

  ## The way to estimate mu and sigma needs some thinking.
  ## However it doesn't seem to make much differences.
  ## NR procedure to do shrinkage
  max.value = max(lphi.g0,10, na.rm=TRUE)
  ## objective function (penalized likelihood)
  get.phi <- function(dat){
    y=dat[1:nsamples]
    Ey=dat[nsamples+(1:nsamples)]
    mu.phi=dat[length(dat)]
    obj=function(phi) {
      alpha=1/phi
      tmp1=1/(1+Ey*phi)
      tmp2=1-tmp1
      -(sum(lgamma(alpha+y)) - nsamples*lgamma(alpha) + alpha*sum(log(tmp1)) + sum(y*log(tmp2)) -
        ((log(phi) - mu.phi)^2) / (2*(sigma^2)) - log(phi) - log(sigma))
    }
    return(optimize(obj, interval=c(0.01, max.value))$minimum)
  }

  tmp=cbind(Y, muY, lOD.pred)
  phi.hat=apply(tmp,1,get.phi)
  phi.hat
}

### compute hyperparameter sigma, when all genes have the same phi.
## this seems over estimated that. Need to think!!!
compute.baseSigma.trend <- function(lOD.pred, Y, muY,k, design) {
  phi0=exp(lOD.pred)
  n=length(phi0)
  m=ncol(Y)
  nsim=5
  res=rep(0, nsim)
  ## sample
  k=rep(1,m)
  lexpr.cut=2
  for(i in 1:nsim) {
    Ysim <- matrix(rnegbinom(length(Y), muY, phi=phi0), ncol=m)
    lexpr.sim=rowMeans(log(sweep(Ysim+0.5,2,k,FUN="/")))
    phi.sim = est.dispersion(Ysim, k, design)
    lOD.sim=log(phi.sim)
    lOD.pred.sim=est.trend(Ysim)
    ix=lexpr.sim>lexpr.cut
    sigma=IQR(lOD.sim[ix]-lOD.pred.sim[ix], na.rm=TRUE) / 1.349
    res[i]=sigma^2
  }
  mean(res)
}

###############################################
## functions for estimating trends
###############################################
scv <- function(x) var(x)/mean(x)^2

## trend estimation. Results are trend in log scale
est.trend <- function(X,xmin=1){
  ss=colSums(X); ss=ss/min(ss)
  X= sweep(X,2,ss,FUN="/")
  phi.hat1 = est.phi0(X)
  lexpr=rowMeans(log(X))
  lexpr.cut=2
  phi0a=median(phi.hat1[lexpr>lexpr.cut])
  SCV=apply(X,1,scv)
  ##phi0 is the m0, average phi0 assuming no trend
  binsize=0.2
  xmax=quantile(lexpr, 0.99, na.rm=TRUE)
  trend1=tapply(phi.hat1-phi0a,cut(lexpr,c(seq(xmin,xmax,binsize),Inf)),median)
  trend1=-isoreg(-trend1)$yf
  fit1=smooth.spline(trend1~seq(xmin,xmax,binsize),df=3)
  pred1=predict(fit1,lexpr)$y
  pred1[lexpr>xmax]=predict(fit1,xmax)$y
  pred1[lexpr<xmin]=pmax(predict(fit1,xmin)$y,pred1[lexpr<xmin])
  m0.adj=pred1+phi0a
  log(m0.adj)
}



