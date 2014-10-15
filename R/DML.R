######################################################################
## a list of functions for testing DML from Bisulfite seq data
##
## These functions are subject to change because they might not work
## for real data when there are too many CG sites.
## Maybe I should only report significant CG sites.
##
######################################################################

######################################
## wrapper function for DML test
######################################
DMLtest <- function(BSobj, group1, group2, equal.disp=FALSE, smoothing=FALSE,
                    smoothing.method=c("ma", "BSmooth"), smoothing.span=500, ...) {
  ## grab two group data
  tmp <- getBSseqIndex(sampleNames(BSobj), group1, group2)
  BS1 <- BSobj[,tmp$group1]
  BS2 <- BSobj[,tmp$group2]

  ## Check the consistence of inputs
  nreps1 <- dim(BS1)[2]
  nreps2<- dim(BS2)[2]
  if( (nreps1==1 | nreps2==1) ) { ## singel replicate case
    if(!smoothing )
      stop("There is no biological replicates. Please set smoothing=TRUE and retry.")
  }

  if(!smoothing) { ## no smoothing.
    dmls <- DMLtest.noSmooth(BS1, BS2, equal.disp=FALSE)
  } else { ## smoothing version
    smoothing.method <- match.arg(smoothing.method)
    dmls <- DMLtest.Smooth(BS1, BS2, equal.disp=FALSE, smoothing.method, smoothing.span, ...)
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
DMLtest.noSmooth <- function(BS1, BS2, equal.disp=FALSE) {
  ## grab counts
  x1 <- getCoverage(BS1, type="M")
  n1 <-getCoverage(BS1, type="Cov")
  x2 <-getCoverage(BS2, type="M")
  n2 <-getCoverage(BS2, type="Cov")
  estprob1 <- compute.mean.noSmooth(x1, n1)
  estprob2 <- compute.mean.noSmooth(x2, n2)

  nreps1 <- ncol(x1)
  nreps2 <- ncol(x2)
  ## estimate priors from counts.
  if(equal.disp | nreps1==1 | nreps2==1) {
    ## assume equal dispersion in two conditions. Combine counts from two conditions and estimate dispersions.
    ## Should keep only those sites didn't show much differences??
    ## I'll ignore that part for now, but presumably there shouldn't be too many of those.
    ## IF there are single replicate I'll estimate prior assuming equal dispersion.
    prior1 <- est.prior.logN(cbind(x1,x2), cbind(n1,n2))
    prior2 <- 0
  } else { # different prior for two conditions
    prior1 <- est.prior.logN(x1, n1)
    prior2 <- est.prior.logN(x2, n2)
  }

  ## perform Wald test
  allchr <- as.character(seqnames(BS1))
  allpos <- start(BS1)
  wald <- waldTest.DML(x1, n1, estprob1, x2, n2, estprob2,
                       prior1, prior2, equal.disp=equal.disp,
                       smoothing=FALSE, allchr=allchr, allpos=allpos)

  return(wald)

}

######################################
## test DML with smoothing
######################################
DMLtest.Smooth <- function(BS1, BS2, equal.disp=FALSE,
                           smoothing.method, smoothing.span, ...) {
  ## grab counts
  x1 <- getCoverage(BS1, type="M")
  n1 <- getCoverage(BS1, type="Cov")
  x2 <- getCoverage(BS2, type="M")
  n2 <- getCoverage(BS2, type="Cov")
  nreps1 <- ncol(x1)
  nreps2 <- ncol(x2)

  ## Smoothing
  cat("Smoothing, this will take a while ...\n")
  mu1 <- compute.mean.Smooth(BS1, smoothing.method, smoothing.span, ...)
  mu2 <- compute.mean.Smooth(BS2, smoothing.method, smoothing.span, ...)
  estprob1 <- matrix(rep(mu1, nreps1), ncol=nreps1)
  estprob2 <- matrix(rep(mu2, nreps2), ncol=nreps2)

  ## estimate priors from counts -
  ## Maybe I should modify the way to estimate prior because this is the smoothing version.
  if(equal.disp | nreps1==1 | nreps2==1) {  ## assume equal dispersion in two conditions
    prior1 <- est.prior.logN(cbind(x1,x2), cbind(n1,n2))
    prior2 <- prior1
  } else { # different prior for two conditions
    prior1 <- est.prior.logN(x1, n1)
    prior2 <- est.prior.logN(x2, n2)
  }

  allchr <- as.character(seqnames(BS1))
  allpos <- start(BS1)
  ## perform Wald test
  wald <- waldTest.DML(x1, n1, estprob1, x2, n2, estprob2,
                       prior1, prior2, equal.disp=equal.disp,
                       smoothing=TRUE, smoothing.span, allchr, allpos)

  return(wald)
}


###############################################################################
## Perform Wald tests for calling DML
###############################################################################
waldTest.DML <- function(x1,n1,estprob1, x2,n2, estprob2, prior1, prior2,
                         equal.disp=equal.disp, smoothing,
                         smoothing.span, allchr, allpos) {

  if(equal.disp)
    prior2 <- prior1

  ## estimated shrunk dispersions
  ## - this part is slow. Should be computed parallely. Will implement later.
  cat("Estimating dispersion for each CpG site, this will take a while ...\n")

  if(equal.disp) { ## equal dispersion. Combine two groups and shrink
    x <- cbind(x1, x2); n <- cbind(n1, n2)
    estprob <- cbind(estprob1, estprob2)
    shrk.phi1 <- shrk.phi2 <- dispersion.shrinkage(x, n, prior1, estprob)

  } else { ## shrink two groups separately
    shrk.phi1 <- dispersion.shrinkage(x1, n1, prior1, estprob1)
    shrk.phi2 <- dispersion.shrinkage(x2, n2, prior2, estprob2)
  }

  ## Wald test
  if(smoothing) {
    wald <- compute.waldStat.Smooth(estprob1[,1], estprob2[,1], n1, n2, shrk.phi1, shrk.phi2,
                                    smoothing.span, allchr, allpos)
  } else
  wald <- compute.waldStat.noSmooth(estprob1[,1], estprob2[,1], n1, n2, shrk.phi1, shrk.phi2)

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
  n1m <- rowSums(n1);    n2m <- rowSums(n2)
  vx1 <- rowSums(n1*estprob1*(1-estprob1)*(1+(n1-1)*phi1))
  vx2 <- rowSums(n2*estprob2*(1-estprob2)*(1+(n2-1)*phi2))
  flag <- "sum"
  n1m.sm <- smooth.chr(as.double(n1m), smoothing.span, allchr, allpos, flag)
  n2m.sm <- smooth.chr(as.double(n2m), smoothing.span, allchr, allpos, flag)
  vx1.sm <- smooth.chr(as.double(vx1), smoothing.span, allchr, allpos, flag)
  vx2.sm <- smooth.chr(as.double(vx2), smoothing.span, allchr, allpos, flag)
  var1 <- vx1.sm / n1m.sm^2
  var2 <- vx2.sm / n2m.sm^2
  vv <- var1 + var2
  ## bound vv a little bit
  vv[vv<1e-5] <- 1e-5
  se <- sqrt(vv)
  stat <- dif/se
  pval <- 2 * (1 - pnorm(abs(stat))) ## p-value for hypothesis testing
  fdr <- p.adjust(pval, method="fdr")

  data.frame(mu1=estprob1, mu2=estprob2, diff=dif, diff.se=se, stat=stat,
             phi1=phi1, phi2=phi2, pval=pval, fdr=fdr)
}

################
## call DML
################
callDML <- function(DMLresult, delta=0, p.threshold=1e-5) {

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


