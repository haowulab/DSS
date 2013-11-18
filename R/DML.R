######################################################################
## a list of functions for DML/DMR detections from Bisulfite seq data
######################################################################
require(bsseq)

######################################
## wrapper function to calling DML. 
## This is without smoothing.
######################################
callDML <- function(BS1, BS2, equal.disp=FALSE, threshold=0) {
  ## first merge data from two conditions according to chr and pos
  n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
  n2 <-getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
  gr1 <- getBSseq(BS1,"gr");  gr2 <- getBSseq(BS2,"gr")
  alldata <- mergeData.counts(n1, x1, gr1, n2, x2, gr2)
  
  
  ## estimate priors from counts.
  if(equal.disp) {
    ## assume equal dispersion in two conditions. Combine counts from two conditions and estimate dispersions.
    ## Should keep only those sites didn't show much differences??
    ## I'll ignore that part for now. 
    ix.X <- grep("X", colnames(alldata))
    x1 <- alldata[,ix.X]
    ix.N <- grep("N", colnames(alldata))
    n1 <- alldata[,ix.N]
    prior1 <- est.prior.logN(x1, n1)
    prior2 <- 0
  } else { # different prior for two conditions
    n1 <- getBSseq(BS1,"Cov"); x1 <- getBSseq(BS1,"M")
    prior1 <- est.prior.logN(x1, n1)
    n2 <- getBSseq(BS2,"Cov"); x2 <- getBSseq(BS2,"M")
    prior2 <- est.prior.logN(x2, n2)
  }
  
  ## grab data 
  cc <- colnames(alldata)
  ix.X1 <- grep("X.*cond1", cc);
  ix.N1 <- grep("N.*cond1", cc)
  ix.X2 <- grep("X.*cond2", cc);
  ix.N2 <- grep("N.*cond2", cc)
  ncol1 <- length(ix.X1); ncol2 <- length(ix.X2)
  x1 <- as.matrix(alldata[,ix.X1]); n1 <- as.matrix(alldata[,ix.N1])
  x2 <- as.matrix(alldata[,ix.X2]); n2 <- as.matrix(alldata[,ix.N2])

  ## compute means. Spatial correlations are ignored at this time
  ## the means need to be of the same dimension as X and N
  estprob1 <- compute.mean(x1, n1);   estprob2 <- compute.mean(x2, n2)

  ## perform Wald test 
  wald <- waldTest.DML(x1, n1, estprob1, x2, n2, estprob2,
                       prior1, prior2, threshold, equal.disp=equal.disp)

  ## combine with chr/pos and output
  result <- data.frame(chr=alldata$chr, pos=alldata$pos, wald)

  ## sort result according to chr and pos
  ix <- sortPos(alldata$chr, alldata$pos)
  
  return(result[ix,])
}


###############################################################################
## Perform Wald tests for calling DML.
###############################################################################
waldTest.DML <- function(x1,n1,estprob1, x2,n2, estprob2, prior1, prior2,
                         threshold, equal.disp=equal.disp) {

  if(equal.disp)
    prior2 <- prior1
  
  ## estimated shrunk dispersions
  ## - this part is slow. Should be computed parallely
  if(equal.disp) { ## equal dispersion. Combine two groups and shrink
    x <- cbind(x1, x2); n <- cbind(n1, n2)
#    estprob <- cbind(matrix(rep(estprob1,ncol(x1)), ncol=ncol(x1)),
#                     matrix(rep(estprob1,ncol(x2)), ncol=ncol(x2)))
    estprob <- cbind(estprob1, estprob2)
    shrk.phi1 <- shrk.phi2 <- dispersion.shrinkage(x, n, prior1, estprob)
    
  } else { ## shrink two groups separately 
    shrk.phi1 <- dispersion.shrinkage(x1, n1, prior1, estprob1)
    shrk.phi2 <- dispersion.shrinkage(x2, n2, prior2, estprob2)
  }

  ## Wald test
  wald <- compute.waldStat(estprob1[,1], estprob2[,1], n1, n2, shrk.phi1, shrk.phi2)

  ## obtain posterior probability that the differnce of two means are greater than a threshold
  if( threshold>0 ) {
    p1 <- pnorm(wald$diff-threshold, sd=wald$diff.se) ## Pr(delta.mu > threshold)
    p2 <- pnorm(wald$diff+threshold, sd=wald$diff.se, lower.tail=FALSE) ## Pr(-delta.mu < -threshold)
    pp.diff <- p1 + p2
    wald <- data.frame(wald, postprob.overThreshold=pp.diff)
  }
  
  return(wald)
}

###############################################
## compute Wald test statistics
## Need to deal with missing data
## Currently work for two conditions.
## Condition means are assumed to be the same among replicates.
###############################################
compute.waldStat <- function(estprob1, estprob2, n1, n2, phi1, phi2) {
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
  
  data.frame(mu1=estprob1, mu2=estprob2, diff=dif, diff.se=se, stat=stat, pval=pval, fdr=fdr)
}
