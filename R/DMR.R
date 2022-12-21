########################################################
##
## function to call DMRs from the DML test results
##
########################################################

callDMR <- function(DMLresult, delta=0, p.threshold=1e-5,
                    minlen=50, minCG=3, dis.merge=100, pct.sig=0.5) {
    ## remove the NA entries
    ix.keep = !is.na(DMLresult$stat)
    if(mean(ix.keep) < 1)  ## with NA entries in the results
        DMLresult = DMLresult[ix.keep,]

    flag.multifactor = FALSE
    if(inherits(DMLresult, "DMLtest.multiFactor"))
        flag.multifactor = TRUE

    if(dis.merge > minlen)
        dis.merge = minlen

    ## deal with delta
    if( delta > 0 ) {
        if(flag.multifactor) { # multifactor, doesn't support delta
            stop("The test results is based on multifactor design, 'delta' is not supported")
        }
        p1 <- pnorm(DMLresult$diff-delta, sd=DMLresult$diff.se) ## Pr(delta.mu > delta)
        p2 <- pnorm(DMLresult$diff+delta, sd=DMLresult$diff.se, lower.tail=FALSE) ## Pr(-delta.mu < -delta)
        postprob.overThreshold <- p1 + p2
        DMLresult <- data.frame(DMLresult, postprob.overThreshold=postprob.overThreshold)
        scores <- 1 - postprob.overThreshold
        ## modify the test statistics according to postprob - this seems make the results worse
        ##        DMLresult[,"stat"] <- qnorm(scores) * sign(DMLresult[,"stat"])
    } else {
        scores <- DMLresult$pval
    }

    ## bump finding
    dmrs <- findBumps(DMLresult$chr, DMLresult$pos, scores,
                      cutoff=p.threshold, sep=5000, dis.merge=dis.merge,
                      pct.sig=pct.sig, minCG=minCG)

    ## compute average methylation levels in two groups
    if(is.null(dmrs)) {
        warning("No DMR found! Please use less stringent criteria. \n")
        return(NULL)
    }

    ## Looping seems to be very slow.
    nCG <- dmrs[,"idx.end.global"] - dmrs[,"idx.start.global"] + 1
    ix.good <- dmrs$length>minlen & nCG>minCG
    if(sum(ix.good) == 0) {
        warning("No DMR found! Please use less stringent criteria. \n")
        return(NULL)
    }
    dmrs <- dmrs[ix.good,]
    nCG <- dmrs[,"idx.end.global"] - dmrs[,"idx.start.global"] + 1

    ## create final result data frame
    if(flag.multifactor) { ## multifactor
        areaStat = rep(0, nrow(dmrs))
        for(i in 1:nrow(dmrs)) { ## this part is kind of slow when number of DMRs is large
            ii=dmrs[i,"idx.start.global"]:dmrs[i,"idx.end.global"]
            areaStat[i] = sum(DMLresult[ii,"stat"])
        }
        result <- data.frame(dmrs[,1:4], nCG=nCG, areaStat=areaStat)

    } else { ## single factor
        meanMethy1 = meanMethy2 = areaStat = rep(0, nrow(dmrs))
        for(i in 1:nrow(dmrs)) { ## this part is kind of slow when number of DMRs is large
            ii=dmrs[i,"idx.start.global"]:dmrs[i,"idx.end.global"]
            meanMethy1[i] = mean(DMLresult[ii,"mu1"])
            meanMethy2[i] = mean(DMLresult[ii,"mu2"])
            areaStat[i] = sum(DMLresult[ii,"stat"])
        }

        result <- data.frame(dmrs[,1:4], nCG=nCG, meanMethy1, meanMethy2,
                             diff.Methy=meanMethy1-meanMethy2, areaStat=areaStat)
    }


    ## sort by areaStat
    ix = sort(abs(result$areaStat), decreasing=TRUE, index.return=TRUE)$ix
    result[ix,]
}


#######################################################
## seperate data into blocks according to chr and pos
#######################################################
findRegion <- function(chr, pos, sep=1000) {
    pos.diff <- abs(c(as.integer(0), diff(pos)))
    idx.jump <- which(pos.diff>sep)
    regions <- rbind(c(1, idx.jump), c(idx.jump-1, length(pos)))
    regions
}

#######################################################
## Bump finding, given score and cutoff.
## This is slow. Need to rewrite in C.
#######################################################
findBumps <- function(chr, pos, x, cutoff, sep=1000, dis.merge=200, pct.sig=0.3, minCG) {
    flag <- as.numeric(x<cutoff)
    if(sum(flag) == 0) ## none
        return(NULL)

    flag[is.na(flag)]=FALSE
    ## find regions
    regions <- findRegion(chr, pos, sep)
    ## loop on regions
    initn <- 100000 ## initialze number of DMRs. Allocate enough rows to start with.
    result <- data.frame(chr=rep("chr1",initn), start=rep(0,initn),
                         end=rep(0, initn), length=rep(0, initn),
                         idx.start.global=rep(0, initn),
                         idx.end.global=rep(0, initn))
    levels(result[,1]) <- unique(chr)
    result.idx <- 0
    for(i in 1:ncol(regions)) {
        idx <- regions[1,i]:regions[2,i]
        if(length(idx) <= minCG) next
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

        ## skip if there are less than minCount probes
##         idx.keep <- (endidx-startidx+1)>=minCG
##         startidx <- startidx[idx.keep]
##         endidx <- endidx[idx.keep]
##         if(length(endidx)==0) next

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
        if(length(ix.good) == 0) next
        else {
            startidx <- startidx[ix.good]
            endidx <- endidx[ix.good]
        }

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

    if(result.idx >= 1)
        result <- result[1:result.idx,]
    else return(NULL)
    result
}


