#######################################################
## some utility functions for BS-seq data
#######################################################

###### make an object of BSseq given count data from several replicates
## The input is a list of  data frames with columns: chr, pos, N, X.
makeBSseqData <- function(dat, sampleNames) {
    n0 <- length(dat)

    if(missing(sampleNames))
        sampleNames <- paste("sample", 1:n0, sep="")

    alldat <- dat[[1]]
    if(any(alldat[,"N"] < alldat[,"X"], na.rm=TRUE))
        stop("Some methylation counts are greater than coverage.\n")


    ix.X <- which(colnames(alldat) == "X")
    ix.N <- which(colnames(alldat) == "N")
    colnames(alldat)[ix.X] <- "X1"
    colnames(alldat)[ix.N] <- "N1"


    if(n0 > 1) { ## multiple replicates, merge data
        for(i in 2:n0) {
            thisdat <- dat[[i]]
            if(any(thisdat[,"N"] < thisdat[,"X"], na.rm=TRUE))
                stop("Some methylation counts are greater than coverage.\n")

            ix.X <- which(colnames(thisdat) == "X")
            ix.N <- which(colnames(thisdat) == "N")
            colnames(thisdat)[c(ix.X,ix.N)] <- paste(c("X", "N"),i, sep="")
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

    ## order CG sites according to positions
    ## This doesn't seem to be necessary since merge will automatically order the sites.
##     idx <- split(1:length(alldat$chr), alldat$chr)
##     M.ordered <- M
##     Cov.ordered <- Cov

##     for( i in seq(along=idx) ) {
##         thisidx = idx[[i]]
##         thispos = alldat$pos[ thisidx ]
##         dd = diff(thispos)
##         if( min(dd)<0 ) { # not ordered
##             warning( paste0("CG positions in chromosome ",  names(idx)[i], " is not ordered. Reorder CG sites.\n") )
##             iii = order(thispos)
##             M.ordered[thisidx, ] <- M[thisidx, ][iii,]
##             Cov.ordered[thisidx, ] <- Cov[thisidx, ][iii,]
##         }
##     }

##     result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M.ordered, Cov=Cov.ordered)

    result <- BSseq(chr=alldat$chr, pos=alldat$pos, M=M, Cov=Cov)

    result
}


########################################################################
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
## adding a small constant could bring trouble when there's no coverage!!!
######################################################################################
compute.mean.noSmooth <- function(X, N) {
    p <- X/N
    ##rowSums <- DelayedArray::rowSums
    const <- mean(p, na.rm=TRUE)
    p <- (rowSums(X)+const)/(rowSums(N)+1)

    nreps <-  ncol(N)
    res <- matrix(rep(p, nreps), ncol = nreps)
    return(res)
}

######################################################################################
## function to compute means when there's smoothing.
## Currently it doesn't do anything except calling smooth.collapse.
## Might add things later.
######################################################################################
compute.mean.Smooth <- function(X, N, allchr, allpos, ws=500) {
    ## rowSums <- DelayedArray::rowSums

    ## collapse the replicates and smoothing
    nreps <-  ncol(N)
    p1 <- X / N
    const <-  mean(p1, na.rm=TRUE) ##  small constant to be added
    if(nreps>1) { ## multiple replicate, collapse them. Add a small constant to bound away from 0/1
        N <- rowSums(N)+1
        X <- rowSums(X)+const
    } else {
        N <- N+0.4
        X <- X+0.4*const
    }
    ## smooth and compute p
    X.sm <- smooth.chr(as.double(X), ws, allchr, allpos)
    N.sm <- smooth.chr(as.double(N), ws, allchr, allpos)
    p <- X.sm / N.sm

    res <- matrix(rep(p, nreps), ncol = nreps)
    return(res)
}

######################################################################################
## Function to perform collapsed smoothing,
## that is, to collapse counts on replicates and then smooth
## This works for one object.
##
## 2/5/2015:
## Need a little more thinking here. How to deal with regions with no coverage???
## adding a small constant could bring trouble when there's no coverage!!!
##
## 3/2/2015: This function is not in use now. I'll remove it later.
######################################################################################
smooth.collapse <- function(BS, method, ws, ...) {
    ## grab counts and collapse
    n1 <- getBSseq(BS,"Cov")
    x1 <- getBSseq(BS,"M")
    nreps <-  ncol(n1)
    p1 <- x1 / n1
    const <-  mean(p1, na.rm=TRUE) ##  small constant to be added

    if(nreps>1) { ## multiple replicate, collapse them. Add a small constant to bound away from 0/1
        n1 <- rowSums(n1)+1
        x1 <- rowSums(x1)+const
        BS <- BSseq(chr=seqnames(BS), pos=start(BS), M=as.matrix(x1), Cov=as.matrix(n1))
    } else {
        n1 <- n1+0.4
        x1 <- x1+0.4*const
    }

    ## smooth, and then return smoothed values
    if(method == "BSmooth")  { ## use BSmooth method for smoothing. This is very slow
        res <- BSmooth(BS, h=ws, ...)
        a <- getMeth(res)[,1]
    } else { ## moving avarege. This is fast and will be default
        X.sm = smooth.chr(as.double(x1), ws, as.character(seqnames(BS)), start(BS))
        N.sm = smooth.chr(as.double(n1), ws, as.character(seqnames(BS)), start(BS))
        a = X.sm / N.sm
    }

    a
}

########################################################################
## estimate dispersion for BS-seq data, given means
########################################################################
est.dispersion.BSseq <- function(X, N, estprob) {
    prior <- est.prior.BSseq.logN(X, N)
    dispersion.shrinkage.BSseq(X, N, prior, estprob)
}

########################################################################
## A function to estimate dipersion prior for BS-seq, assuming log-normal prior.
## It takes X and N, and only use the sites with big coverages,
## then return the mean and sd of prior distribution.
##
## For single rep data: will use logN(-3,1) as prior.
########################################################################
est.prior.BSseq.logN <- function(X, N) {
    ## rowMeans = DelayedArray::rowMeans
    ## rowSums = DelayedArray::rowSums

    if(ncol(X) == 1) ## single rep
        return(c(-3, 1))

    ## keep sites with large coverage and no missing data
    ix=rowMeans(N>10)==1 & rowSums(N==0)==0
    if(sum(ix) < 50) {
        warning("The coverages are too low. Cannot get good estimations of prior. Use arbitrary prior N(-3,1).")
        return(c(-3, 1))
    }

    X=X[ix,,drop=FALSE]; N=N[ix,,drop=FALSE]
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

########################################################################
## Dispersion shrinkage based on log-normal penalized likelihood.
## Takes X, N, estimated mean and prior.
##
## The shrinakge is done in log scale. So data will be shrink to the
## logarithmic means.
########################################################################
dispersion.shrinkage.BSseq <- function(X, N, prior, estprob) {
    ## penalized likelihood function
    plik.logN <- function(size, X,mu,m0,tau,phi)
        -(sum(dbb(size, X, mu, exp(phi))) + dnorm(phi, mean=m0, sd=tau, log=TRUE))

    ## for CG sites with no coverage, use prior
    shrk.phi=exp(rep(prior[1],nrow(N)))

    ## deal with estprob, make it a matrix if not.
    if(!is.matrix(estprob))
        estprob <- as.matrix(estprob)

    ## skip those without coverage
    ix <- rowSums(N>0) > 0
    X2 <- X[ix, ,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; estprob2 <- estprob[ix,,drop=FALSE]
    shrk.phi2 <- rep(0, nrow(X2))

    ## setup a progress bar
    nCG.pb = round(nrow(X2)/100)
    pb <- txtProgressBar(style = 3)
    for(i in 1:nrow(X2)) {
        ## print a progress bar
        if((i %% nCG.pb) == 0)
            setTxtProgressBar(pb, i/nrow(X2))
        ## I can keep the 0's with calculation. They don't make any difference.
        shrk.one=optimize(f=plik.logN, size=N2[i,], X=X2[i,], mu=estprob2[i,], m0=prior[1], tau=prior[2],
        interval=c(-5, log(0.99)),tol=1e-3)
        shrk.phi2[i]=exp(shrk.one$minimum)
    }
    setTxtProgressBar(pb, 1)
    cat("\n")
    shrk.phi[ix] <- shrk.phi2

    return(shrk.phi)
}

#########################################################
## beta-binomial (BB) density function.
## The BB distribution is parametrized by mean and dispersion.
#########################################################
dbb <- function (size, x, mu, phi, log=TRUE)  {
    ## 'size' and/or 'x' could be DelayedArray objects so turn them into
    ## ordinary arrays
    size=as.array(size)
    x=as.array(x)
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



