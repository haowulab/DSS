########################################################
## function to call DMR from the DML test results
########################################################
#### function to call DMR
callDMR <- function(DMLresult, p.threshold, minlen=50, minCG=3, dis.merge=50, pct.sig=0.5) {

  if(dis.merge > minlen)
    dis.merge = minlen
  dmrs = findBumps(DMLresult$chr, DMLresult$pos, DMLresult$pval,
    cutoff=p.threshold, sep=1000,dis.merge=dis.merge,
    pct.sig=pct.sig)
  
  ## compute average methylation levels in two groups
  if(nrow(dmrs)==0) {
    warning("No DMR found!\n")
    return(NULL)
  }

  nCG=meanMethy1=meanMethy2=rep(0, nrow(dmrs))
  for(i in 1:nrow(dmrs)) {
    ii=dmrs[i,"idx.start.global"]:dmrs[i,"idx.end.global"]
    nCG[i]=length(ii)
    meanMethy1[i] = mean(DMLresult[ii,"mu1"])
    meanMethy2[i] = mean(DMLresult[ii,"mu2"])
  }
  result= data.frame(dmrs[,1:4], nCG=nCG, meanMethy1, meanMethy2, diff.Methy=meanMethy1-meanMethy2)

  ## discard really short ones
  ix.good <- result$length>minlen & result$nCG>minCG
  result[ix.good,]
}


### seperate data into blocks according to chr and pos
findRegion <- function(chr, pos, sep=1000) {
  pos.diff <- abs(c(as.integer(0), diff(pos)))
  idx.jump <- which(pos.diff>sep)
  regions <- rbind(c(1, idx.jump), c(idx.jump-1, length(pos)))
  regions
}

## Bump finding, given score and cutoff.
findBumps <- function(chr, pos, x, cutoff, sep=1000,dis.merge=100, pct.sig=0.3) {
  flag <- as.numeric(x<cutoff)
  flag[is.na(flag)]=FALSE
  ## find regions
  regions <- findRegion(chr, pos, sep)
  ## loop on regions
  initn <- 50000
  result <- data.frame(chr=rep("chr1",initn), start=rep(0,initn),
                       end=rep(0, initn), length=rep(0, initn),
                       idx.start.global=rep(0, initn),
                       idx.end.global=rep(0, initn))
  levels(result[,1]) <- unique(chr)
  result.idx <- 0
  for(i in 1:ncol(regions)) {
    idx <- regions[1,i]:regions[2,i]
    pos.region <- pos[idx]
    nn <- length(idx)
    flag.region <- flag[idx]
    ## get start/end position
    startidx <- which(flag.region[-nn]==0 & flag.region[-1]==1)+1
    if(flag.region[1]==1)
      startidx <- c(1, startidx)
    if(length(startidx)==0)
      next
    endidx <- which(flag.region[-nn]==1 & flag.region[-1]==0)
    if(flag.region[nn]==1)
      endidx <- c(endidx, nn)
    
    ## remove if there are less than minCount probes
##     idx.keep <- (endidx-startidx+1)>=minCount
##     startidx <- startidx[idx.keep]
##                 endidx <- endidx[idx.keep]
##     if(length(endidx)==0) next

    ## merge if they are really close
    nbump <- length(startidx)
    if(nbump>1) {
      bumppos <- cbind(pos[idx][startidx], pos[idx][endidx])
      dis <- bumppos[-1,1]>(bumppos[-nbump,2]+dis.merge)
      idx.start <- which(c(1,dis)==1)
      idx.end <- which(c(dis,1)==1)
      ## merged
      startidx <- startidx[idx.start]
      endidx <- endidx[idx.end]
    }
    
    ## after merging, make sure certain % of CG sites being significant is enough.
    ## this is not easy!!!
    ix.good <- NULL
    x.thisregion <- x[idx]
    for(ibump in 1:length(startidx)) {
      ii <- startidx[ibump]:endidx[ibump]
      pp <- x.thisregion[ii] < cutoff
      if(mean(pp) > pct.sig)
        ix.good <- c(ix.good, ibump)
    }
    if(length(ix.good) > 1) {
      startidx <- startidx[ix.good]
      endidx <- endidx[ix.good]
    } else next
    
    nbump <- length(startidx)
    ll <- pos.region[endidx] - pos.region[startidx] + 1
    tmpn <- length(ll)

    ## make result
    result[result.idx+(1:tmpn),] <- data.frame(chr=as.character(chr[idx][startidx]),
                                               start=pos[idx][startidx],
                                               end=pos[idx][endidx], length=ll,
                                               idx.start.global=idx[startidx],
                                               idx.end.global=idx[endidx])
    result.idx <- result.idx + tmpn
  }
  
  result <- result[1:result.idx,]
  ## remove really short ones
##  result <- result[result[,4]>minlen,]
  result
}


