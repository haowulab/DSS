######################################################################
##
## a list of functions for testing DML from Bisulfite seq data
##
######################################################################

######################################
## wrapper function for DML test
######################################
DMLtest <- function(BSobj, group1, group2, equal.disp=FALSE, smoothing=FALSE, smoothing.span=500) {
    ## grab two group data
    tmp <- getBSseqIndex(sampleNames(BSobj), group1, group2)
    BS1 <- BSobj[,tmp$group1]
    BS2 <- BSobj[,tmp$group2]

    ## remove loci with all 0 coverages in a condition
    ## It's not required that all replicates have coverage.
    ## But there must be some coverage from at least one replicate.
    n1 <- getBSseq(BS1, "Cov")
    n2 <- getBSseq(BS2, "Cov")
    ## remove if a rather long strech (like 200 bps) of regions have no coverage
    allpos <- start(BSobj)
    ix1 <- hasCoverage(n1, allpos)
    ix2 <- hasCoverage(n2, allpos)
    ix <- ix1 & ix2
    BS1 <- BS1[ix]
    BS2 <- BS2[ix]

    ## Check the consistence of inputs
    nreps1 <- dim(BS1)[2]
    nreps2<- dim(BS2)[2]
    if( (nreps1==1 | nreps2==1) & !equal.disp ) { ## singel replicate case, and unequal dispersion in two groups
        if( !smoothing )
            stop("There is no biological replicates in at least one condition. Please set smoothing=TRUE or equal.disp=TRUE and retry.")
    }

    if(!smoothing) { ## no smoothing.
        dmls <- DMLtest.noSmooth(BS1, BS2, equal.disp)
    } else { ## smoothing version
        dmls <- DMLtest.Smooth(BS1, BS2, equal.disp, smoothing.span)
    }
    dmls
}

##################################################
## determine sample index from a BSseq object
##################################################
getBSseqIndex <- function(sName, group1, group2) {
    check <- function(group, id) {
        thisGrp = paste("group", id, sep="")

        if (is.character(group)) {
            if(!all(group %in% sName))
                stop("Some sample names not found in", thisGrp)
            group <- match(group, sName)
        }
        if (is.numeric(group)) {
            if(min(group) < 1 | max(group) > length(sName))
                stop("Some group indices are wrong in", thisGrp)
        }
        else stop(paste("problems with argument", thisGrp))
        group
    }

    group1 <- check(group1, 1)
    group2 <- check(group2, 2)

    if(length(intersect(group1, group2)) > 0)
        stop("Two groups have common sample.")

    if(length(group1) <= 0)
        stop("group1 not found.")
    if(length(group2) <= 0)
        stop("group2 not found.")

    list(group1=group1, group2=group2)
}


######################################
## test DML without smoothing
######################################
DMLtest.noSmooth <- function(BS1, BS2, equal.disp) {
    ## grab counts
    x1 <- getCoverage(BS1, type="M")
    n1 <-getCoverage(BS1, type="Cov")
    x2 <-getCoverage(BS2, type="M")
    n2 <-getCoverage(BS2, type="Cov")
    nreps1 <- ncol(x1)
    nreps2 <- ncol(x2)

    ## estimate means
    estprob1 <- compute.mean.noSmooth(x1, n1)
    estprob2 <- compute.mean.noSmooth(x2, n2)

    ## estimate dispersion
    ## - this part is slow. Could be computed parallely. Will implement later.
    cat("Estimating dispersion for each CpG site, this will take a while ...\n")
    if(equal.disp | nreps1==1 | nreps2==1) {
        phi1 <- phi2 <- est.dispersion.BSseq(cbind(x1,x2), cbind(n1,n2), cbind(estprob1, estprob2))
    } else {
        phi1 <- est.dispersion.BSseq(x1, n1, estprob1)
        phi2 <- est.dispersion.BSseq(x2, n2, estprob2)
    }

    ## weight the counts
    wt1 <- 1 / (1+(n1-1)*phi1);    wt1 <- wt1 / mean(wt1)
    wt2 <- 1 / (1+(n2-1)*phi2);    wt2 <- wt2 / mean(wt2)
    x1.wt <- x1*wt1
    n1.wt <- n1*wt1
    x2.wt <- x2*wt2
    n2.wt <- n2*wt2

    ## re-estimate means
    estprob1 <- compute.mean.noSmooth(x1.wt, n1.wt)
    estprob2 <- compute.mean.noSmooth(x2.wt, n2.wt)

    ## perform Wald test
    allchr <- as.character(seqnames(BS1))
    allpos <- start(BS1)
    wald <- waldTest.DML(x1.wt, n1.wt, estprob1, phi1, x2.wt, n2.wt, estprob2, phi2,
                         smoothing=FALSE, allchr=allchr, allpos=allpos)

    return(wald)

}

######################################
## test DML with smoothing
######################################
DMLtest.Smooth <- function(BS1, BS2, equal.disp, smoothing.span) {
    ## grab counts
    x1 <- getCoverage(BS1, type="M")
    n1 <- getCoverage(BS1, type="Cov")
    x2 <- getCoverage(BS2, type="M")
    n2 <- getCoverage(BS2, type="Cov")
    nreps1 <- ncol(x1)
    nreps2 <- ncol(x2)
    allchr <- as.character(seqnames(BS1))
    allpos <- start(BS1)

    ## Smoothing
    cat("Smoothing ...\n")
    estprob1 <- compute.mean.Smooth(x1, n1, allchr, allpos, smoothing.span)
    estprob2 <- compute.mean.Smooth(x2, n2, allchr, allpos, smoothing.span)

    ## estimate priors from counts - need to rethink the single rep case.
    cat("Estimating dispersion for each CpG site, this will take a while ...\n")
    if(equal.disp) {
        phi1 <- phi2 <- est.dispersion.BSseq(cbind(x1,x2), cbind(n1,n2), cbind(estprob1, estprob2))
    } else {
        phi1 <- est.dispersion.BSseq(x1, n1, estprob1)
        phi2 <- est.dispersion.BSseq(x2, n2, estprob2)
    }

    ## update counts - weight by dispersion
    wt1 <- 1 / (1+(n1-1)*phi1); wt1 <- wt1 / mean(wt1)
    wt2 <- 1 / (1+(n2-1)*phi2); wt2 <- wt2 / mean(wt2)
    x1.wt <- x1*wt1
    n1.wt <- n1*wt1
    x2.wt <- x2*wt2
    n2.wt <- n2*wt2

    ## re-estimate means
    estprob1 <- compute.mean.Smooth(x1.wt, n1.wt, allchr, allpos, smoothing.span)
    estprob2 <- compute.mean.Smooth(x2.wt, n2.wt, allchr, allpos, smoothing.span)

    wald <- waldTest.DML(x1.wt, n1.wt, estprob1, phi1, x2.wt, n2.wt, estprob2, phi2,
                         smoothing=TRUE, smoothing.span, allchr=allchr, allpos=allpos)
##     wald <- waldTest.DML(x1, n1, estprob1, phi1, x2, n2, estprob2, phi2,
##                          smoothing=TRUE, smoothing.span, allchr=allchr, allpos=allpos)

    return(wald)
}


###############################################################################
## Perform Wald tests for calling DML
###############################################################################
waldTest.DML <- function(x1,n1,estprob1, phi1, x2,n2, estprob2, phi2, smoothing,
                         smoothing.span, allchr, allpos) {

    ## Wald test
    if(smoothing) {
        wald <- compute.waldStat.Smooth(estprob1[,1], estprob2[,1], n1, n2, phi1, phi2,
                                        smoothing.span, allchr, allpos)
    } else {
        wald <- compute.waldStat.noSmooth(estprob1[,1], estprob2[,1], n1, n2, phi1, phi2)
    }

    ## combine with chr/pos and output
    result <- data.frame(chr=allchr, pos=allpos, wald)
    ## remove NA entries - Maybe I should keep them so result have the same dimension as inputs???
    ii <- !is.na(result$stat)
    result <- result[ii,]

    ## sort result according to chr and pos - maybe this is not important??
    ix <- sortPos(result$chr, result$pos)
    result <- result[ix,]
    return(result)
}

###########################################################
## compute Wald test statistics when there's no smoothing
###########################################################
compute.waldStat.noSmooth <- function(estprob1, estprob2, n1, n2, phi1, phi2) {
    dif <- estprob1 - estprob2
    n1m <- rowSums(n1);    n2m <- rowSums(n2)
    var1 <- rowSums(n1*estprob1*(1-estprob1)*(1+(n1-1)*phi1)) / (n1m)^2
    var2 <- rowSums(n2*estprob2*(1-estprob2)*(1+(n2-1)*phi2)) / (n2m)^2
    ##vv <- var1/ncol1+var2/ncol2
    vv <- var1 + var2
    ## bound vv a little bit??
    vv[vv<1e-5] <- 1e-5
    se <- sqrt(vv)
    stat <- dif/se
    pval <- 2 * (1 - pnorm(abs(stat))) ## p-value for hypothesis testing
    fdr <- p.adjust(pval, method="fdr")

    data.frame(mu1=estprob1, mu2=estprob2, diff=dif, diff.se=se, stat=stat,
               phi1=phi1, phi2=phi2, pval=pval, fdr=fdr)
}

###########################################################
## compute Wald test statistics when there is smoothing
## Note that the variance computation is different when there's smoothing!!!
## The variances will be smaller in the CG dense regions.
## This is reasonable because there're more data points in smoothing.
## But does it make sense biologically???
###########################################################
compute.waldStat.Smooth <- function(estprob1, estprob2, n1, n2, phi1, phi2, smoothing.span,
                                    allchr, allpos) {
    dif <- estprob1 - estprob2
    ## compute variances for moving average values. This is tricky!
    var1 <- compute.var.Smooth(estprob1, n1, phi1, smoothing.span, allchr, allpos)
    var2 <- compute.var.Smooth(estprob2, n2, phi2, smoothing.span, allchr, allpos)

    ##  var1 <- compute.var.Smooth.old(estprob1, n1, phi1, smoothing.span, allchr, allpos)
    ##  var2 <- compute.var.Smooth.old(estprob2, n2, phi2, smoothing.span, allchr, allpos)

    vv <- var1 + var2
    ## bound vv a little bit
    vv[vv<1e-5] <- 1e-4
    vv[vv> 1- 1e-5] <- 1-1e-4

    se <- sqrt(vv)
    stat <- dif/se
    pval <- 2 * (1 - pnorm(abs(stat))) ## p-value for hypothesis testing
    fdr <- p.adjust(pval, method="fdr")

    data.frame(mu1=estprob1, mu2=estprob2, diff=dif, diff.se=se, stat=stat,
               phi1=phi1, phi2=phi2, pval=pval, fdr=fdr)
}


################################################################
## function to compute the variance for moving average values
################################################################
compute.var.Smooth <- function(estprob1, n1, phi1, smoothing.span, allchr, allpos) {
    nreps <- ncol(n1)
    ## estimate distance-dependent autocorrelations
    ## rhos <- est.rho(estprob1, allchr, allpos)
    ## autocorrelation.
    rho <- 0.8
    ## rho <- acf(estprob1, 1, plot=FALSE)$acf[2]  ## lag-1 autocorrelation, use 0.8 fixed
    ## compute vars for each replicate
    vars.rep <- matrix(0, nrow=nrow(n1), ncol=nreps)
    idx <- split(1:length(allchr), allchr)

    ## note the smoothed value were calculated with a small constant added to N and X.
    ## Need to account for that!!!
    const <- 1/ncol(n1)
    n1 <- n1 + const

    for(irep in 1:nreps) { ## loop on replicates
        ## variances at each position
        vars <- n1[,irep] * estprob1*(1-estprob1) * (1+(n1[,irep]-1)*phi1)
        ## 'vars' could be a DelayedArray object so turn it into an ordinary
        ## array
        vars <- as.vector(vars)
        ## compute covariances - do by chr
        tmp1 <- estprob1*(1-estprob1)*phi1
        vars.smooth <- rep(0, length(allchr))
        for(i in seq(along=idx)) { ## loop on chromosomes
            thisidx <- idx[[i]]
            vars.smooth[idx[[i]]]=.Call( "compute_var_smooth", vars[thisidx], tmp1[thisidx],
                       as.double(n1[thisidx,irep]), as.integer(allpos[thisidx]),
                       as.integer(smoothing.span), as.double(rho) )
        }
        vars.rep[,irep] <- vars.smooth
    }
    ## compute denominators
    flag <- "sum"
    n1.sm <- as.array(n1)
    for(i in 1:nreps)
        n1.sm[,i] <- smooth.chr(as.double(n1[,i]), smoothing.span, allchr, allpos, flag)
    denom <- rowSums(n1.sm) ^ 2
    ## results
    vars <- rowSums(vars.rep) / denom
    vars
}



################################################################
## old function to compute the variance for moving average values.
## This uses an approximation that will speed up the calculation.
################################################################
compute.var.Smooth.old <- function(estprob1, n1, phi1, smoothing.span, allchr, allpos) {
    ## variance of X at each position for each replicate
    var1.X <- (n1*estprob1*(1-estprob1)*(1+(n1-1)*phi1))

    ## Consider the smoothing effect
    n1.sm <- n1
    var1.X.sm <- var1.X
    nCG1 <- n1

    flag <- "sum"
    for(i in 1:ncol(n1)) {
        n1.sm[,i] <- smooth.chr(as.double(n1[,i]), smoothing.span, allchr, allpos, flag)
        var1.X.sm[,i] <- smooth.chr(as.double(var1.X[,i]), smoothing.span, allchr, allpos, flag)
        nCG1[,i] <- smooth.chr(rep(1.0, nrow(n1)), smoothing.span, allchr, allpos, flag)
    }
    ## adjust for correlation - this is very crude, but serve the purpose
    ## Need to work this out !!!!
    var1.X.sm2 <- var1.X.sm #* nCG1 * 0.1
    n1m <- rowSums(n1.sm)
    var1 <- rowSums(var1.X.sm2) / n1m^2
    var1
}

################################################
## wrapper function for calling DML
################################################
callDML <- function(DMLresult, delta=0.1, p.threshold=1e-5) {

    ## obtain posterior probability that the differnce of two means are greater than a threshold
    if( delta > 0 ) {
        p1 <- pnorm(DMLresult$diff-delta, sd=DMLresult$diff.se) ## Pr(delta.mu > delta)
        p2 <- pnorm(DMLresult$diff+delta, sd=DMLresult$diff.se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
        postprob.overThreshold <- p1 + p2
        DMLresult <- data.frame(DMLresult, postprob.overThreshold=postprob.overThreshold)
        scores <- 1 - postprob.overThreshold
    } else {
        scores <- DMLresult$pval
    }

    ix <- scores < p.threshold
    result <- DMLresult[ix,]

    ## sort by score
    ii <- sort(scores[ix], index.return=TRUE)$ix
    result[ii,]

}


####################################################################
## function to determine what loci to keep, based on coverage depth
####################################################################
hasCoverage <- function(nn, allpos, thresh=2) {
    nn2 <- rowSums(nn)
    ws <- 200
    flag <- 0
    nn.sm <- .Call("windowFilter", as.double(nn2), as.integer(allpos), as.integer(ws), as.integer(flag))
    nn.sm > thresh*ncol(nn)
}
