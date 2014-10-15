#######################################################
## some utility functions for BS-seq data
#######################################################

###### make an object of BSseq given count data from several replicates
## The input is a list of  data frames with columns: chr, pos, Ntotal, Nmethy.
makeBSseqData <- function(dat, sampleNames) {
  n0 <- length(dat)

  if(missing(sampleNames))
    sampleNames <- paste("sample", 1:n0, sep="")

  alldat <- dat[[1]]
  colnames(alldat)[3:4] <- c("N1", "X1")
  if(any(alldat[,3] < alldat[,4], na.rm=TRUE))
      stop("Some methylation counts are greater than coverage.\n")

  if(n0 > 1) { ## multiple replicates, merge data
    for(i in 2:n0) {
      thisdat <- dat[[i]]
      if(any(alldat[,3] < alldat[,4], na.rm=TRUE))
          stop("Some methylation counts are greater than coverage.\n")
      colnames(thisdat)[3:4] <- paste(c("N", "X"),i, sep="")
      alldat <- merge(alldat, thisdat, all=TRUE)
    }
  }

  ## make BSseq object
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  M <- as.matrix(alldat[,ix.X, drop=FALSE])
  Cov <- as.matrix(alldat[,ix.N, drop=FALSE])
  colnames(M) <- colnames(Cov) <- sampleNames

  result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M, Cov=Cov)

  result
}


########################################################################
## A function to estimate prior parameters, assume log-normal prior.
## It takes X and N, and only use the sites with big coverages,
## then return the mean and sd of prior distribution.
## Potentially, this should take mean if smoothing is allowed.
########################################################################
est.prior.logN <- function(X, N) {
  ## keep sites with large coverage and no missing data
  ix=rowMeans(N>10)==1 & rowSums(N==0)==0
  if(sum(ix) == 0) {
    warning("The coverages are too low. Cannot get good estimations of prior. Use arbitrary prior N(-3,1).")
    return(c(-3, 1))
  }
  X=X[ix,]; N=N[ix,]
  ## compute sample mean/var
  p=X/N
  mm=rowMeans(p)
  mm[mm==0]=1e-5
  mm[mm==1]=1-1e-5
  vv=rowVars(p)
  phi=vv/mm/(1-mm)
  ## exclude those with vv==0. Those are sites with unobservable phis.
  ## But this will over estimate the prior.
  ## What will be the consequences????
  phi=phi[vv>0]
  lphi=log(phi[phi>0])
  prior.mean=median(lphi, na.rm=TRUE)
  prior.sd=IQR(lphi, na.rm=TRUE) /1.39

  ## It seems this over-estimates the truth. Need to use the tricks in
  ## my biostat paper to remove the over-estimation. To be done later.
  c(prior.mean, prior.sd)
}

####################################################
## naive estimates of dispersion, give X and N
####################################################
est.phi.naive <- function(X, N) {
  p=X/N
  mm=rowMeans(p)
  mm[mm==0]=1e-5
  mm[mm==1]=1-1e-5
  vv=rowVars(p)
  phi=vv/mm/(1-mm)
  phi[phi>=1]=1-1e-5
  phi[phi==0]=1e-5
  phi
}

########################################
## Merge two sets of counts.
## This is the version without smoothing
########################################
mergeData.counts <- function(BS1, BS2) {
  n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
  n2 <-getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
  gr1 <- getBSseq(BS1,"gr");  gr2 <- getBSseq(BS2,"gr")
  allchr1 <- seqnames(gr1); allpos1 <- start(gr1)
  allchr2 <- seqnames(gr2); allpos2 <- start(gr2)

  colnames(x1) <- paste(paste("X", 1:ncol(x1), sep=""), "cond1", sep=".")
  colnames(n1) <- paste(paste("N", 1:ncol(n1), sep=""), "cond1", sep=".")
  colnames(x2) <- paste(paste("X", 1:ncol(x2), sep=""), "cond2", sep=".")
  colnames(n2) <- paste(paste("N", 1:ncol(n2), sep=""), "cond2", sep=".")

  dat1 <- data.frame(chr=as.character(seqnames(gr1)), pos=start(gr1), x1, n1)
  dat2 <- data.frame(chr=as.character(seqnames(gr2)), pos=start(gr2), x2, n2)
  alldat <- merge(dat1, dat2, all=FALSE) ## only keep the ones with data in both samples
  alldat[is.na(alldat)] <- 0

  return(alldat)
}

########################################
## Merge two sets of counts.
## This is the version with smoothing
########################################
mergeData.counts.smooth <- function(BS1, BS2, mu1.smooth, mu2.smooth) {
  n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
  n2 <-getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
  gr1 <- getBSseq(BS1,"gr");  gr2 <- getBSseq(BS2,"gr")
  allchr1 <- seqnames(gr1); allpos1 <- start(gr1)
  allchr2 <- seqnames(gr2); allpos2 <- start(gr2)

  colnames(x1) <- paste(paste("X", 1:ncol(x1), sep=""), "cond1", sep=".")
  colnames(n1) <- paste(paste("N", 1:ncol(n1), sep=""), "cond1", sep=".")
  colnames(x2) <- paste(paste("X", 1:ncol(x2), sep=""), "cond2", sep=".")
  colnames(n2) <- paste(paste("N", 1:ncol(n2), sep=""), "cond2", sep=".")

  dat1 <- data.frame(chr=as.character(seqnames(gr1)), pos=start(gr1), x1, n1, mu1=mu1.smooth)
  dat2 <- data.frame(chr=as.character(seqnames(gr2)), pos=start(gr2), x2, n2, mu2=mu2.smooth)
  alldat <- merge(dat1, dat2, all=FALSE) ## only keep the ones with data in both samples
  alldat[is.na(alldat)] <- 0

  return(alldat)
}

########################################
## the rowVars function
########################################
rowVars <- function (x, center = NULL, ...) {
  n <- !is.na(x)
  n <- rowSums(n)
  n[n <= 1] <- NA
  if (is.null(center)) {
    center <- rowMeans(x, ...)
  }
  x <- x - center
  x <- x * x
  x <- rowSums(x, ...)
  x <- x/(n - 1)
  x
}


######################################################################################
## function to compute means when there's no smoothing.
## Just compute the percentage of methylation.
######################################################################################
compute.mean.noSmooth <- function(X, N) {
  p <- (rowSums(X)+0.5)/(rowSums(N)+1)
  res <- matrix(rep(p, ncol(X)), ncol = ncol(X))
  return(res)
}

######################################################################################
## function to compute means when there's smoothing.
## Currently it doesn't do anything except calling smooth.collapse.
## Might add things later.
######################################################################################
compute.mean.Smooth <- function(BS, method=c("ma", "BSmooth"), ws=500, ...) {
  method <- match.arg(method)

  ## collapse the replicates and smoothing
  p <- smooth.collapse(BS, method, ws, ...)
  return(p)
}

######################################################################################
## Function to perform collapsed smoothing,
## that is, to collapse counts on replicates and then BSmooth
## This works for one object
######################################################################################
smooth.collapse <- function(BS, method, ws, ...) {
  ## grab counts and collapse
  n1 <- getBSseq(BS,"Cov")
  x1 <- getBSseq(BS,"M")
  nreps <-  ncol(n1)

  if(nreps>1) { ## multiple replicate, collapse them. Add a small constant to bound away from 0/1
    n1 <- rowSums(n1)+1
    x1 <- rowSums(x1)+0.5
    BS <- BSseq(chr=seqnames(BS), pos=start(BS), M=as.matrix(x1), Cov=as.matrix(n1))
  } else {
    n1 <- n1+0.4
    x1 <- x1+0.2
  }

  ## smooth, and then return smoothed values
  if(method == "BSmooth")  {
    res <- BSmooth(BS, h=ws, ...)
    a <- getMeth(res)[,1]
  } else { ## moving avarege
    X.sm = smooth.chr(as.double(x1), ws, as.character(seqnames(BS)), start(BS))
    N.sm = smooth.chr(as.double(n1), ws, as.character(seqnames(BS)), start(BS))
    a = X.sm / N.sm
  }

  a
}


########################################################################
## Dispersion shrinkage based on log-normal penalized likelihood.
## Takes X, N, estimated mean and prior.
##
## The shrinakge is done in log scale. So data will be shrink to the
## logarithmic means.
########################################################################
dispersion.shrinkage <- function(X, N, prior, estprob) {
  ## penalized likelihood function
  plik.logN <- function(size, X,mu,m0,tau,phi) {
    -(sum(dbb(size, X, mu, exp(phi))) + dnorm(phi, mean=m0, sd=tau, log=TRUE))
  }
  ## for CG sites with no coverage, use prior
  shrk.phi=exp(rep(prior[1],nrow(N)))

  ## deal with estprob, make it a matrix if not.
  if(!is.matrix(estprob))
    estprob <- as.matrix(estprob)

  ## skip those without coverage
  ix <- rowSums(N>0) > 0
  X2 <- X[ix, ,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; estprob2 <- estprob[ix,,drop=FALSE]
  shrk.phi2 <- rep(0, nrow(X2))
  for(i in 1:nrow(X2)) {
    ## I can keep the 0's with calculation. They don't make any difference.
    shrk.one=optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
      interval=c(-5, log(0.99)),tol=1e-3)
    shrk.phi2[i]=exp(shrk.one$minimum)
  }
  shrk.phi[ix] <- shrk.phi2

  return(shrk.phi)
}

#########################################################
## beta-binomial (BB) density function.
## The BB distribution is parametrized by mean and dispersion.
#########################################################
dbb <- function (size, x, mu, phi, log=TRUE)  {
  ## first convert mu/phi to alpha/beta
  tmp=1/phi-1
  alpha=mu*tmp
  beta=tmp - alpha
  v=lchoose(size,x)-lbeta(beta, alpha)+lbeta(size-x + beta,x+alpha)
  if(!log)
    return(exp(v))
  else return(v)
}

###########################################
## sort according to chr and pos.
## Return the index
###########################################
sortPos <- function(chr, pos) {
  do.call(order, list(chr, pos))
}



