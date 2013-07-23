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

  
  if(n0 > 1) { ## multiple replicates, merge data
    for(i in 2:n0) {
      thisdat <- dat[[i]]
      colnames(thisdat)[3:4] <- paste(colnames(thisdat)[3:4], i, sep="")
      alldat <- merge(alldat, thisdat, all=TRUE)
    }
  }

  ## make BSseq object
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=as.matrix(alldat[,ix.X, drop=FALSE]),
                  Cov=as.matrix(alldat[,ix.N, drop=FALSE]))
  
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
## Merge two sets of counts
########################################
mergeData.counts <- function(n1, x1, gr1, n2, x2, gr2) {
  allchr1 <- seqnames(gr1); allpos1 <- start(gr1)
  allchr2 <- seqnames(gr2); allpos2 <- start(gr2)

  colnames(x1) <- paste(paste("X", 1:ncol(x1), sep=""), "cond1", sep=".")
  colnames(n1) <- paste(paste("N", 1:ncol(n1), sep=""), "cond1", sep=".")
  colnames(x2) <- paste(paste("X", 1:ncol(x2), sep=""), "cond2", sep=".")
  colnames(n2) <- paste(paste("N", 1:ncol(n2), sep=""), "cond2", sep=".")

  
  
  dat1 <- data.frame(chr=as.character(seqnames(gr1)), pos=start(gr1), x1, n1)
  dat2 <- data.frame(chr=as.character(seqnames(gr2)), pos=start(gr2), x2, n2)
  alldat <- merge(dat1, dat2, all=TRUE)
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
## function to compute means. If it's based on each CG site, just take an average.
## If include smoothing, use BSmoooth.
######################################################################################
compute.mean <- function(X, N) {
  p <- (rowSums(X)+1)/(rowSums(N)+2)
  res <- matrix(rep(p, ncol(X)), ncol=ncol(X))
  return(p)
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
  shrk.phi=rep(NA,nrow(N))

  ## deal with estprob, make it a matrix
  estprob <- as.matrix(estprob)
  
  ## skip those without coverage, or no replicates.
  ix <- rowSums(N>0) > 1
  X2 <- X[ix, ,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; estprob2 <- estprob[ix,,drop=FALSE]
  
  for(i in 1:nrow(X2)) {
    ## I can keep the 0's with calculation. They don't make any difference.
    shrk.one=optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
      interval=c(-5, log(0.99)),tol=1e-4)
    shrk.phi[i]=exp(shrk.one$minimum)
  }
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
    
