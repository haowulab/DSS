###################################
## some plotting functions for DMR
###################################


###################################
## plot one DMR.
###################################

showOneDMR <- function(OneDMR, BSobj, ext=500, ylim=c(0,1)) {
    ## get chr, position and counts - could this be slow if data is huge??
    allchr = as.character(seqnames(BSobj))
    allpos = start(BSobj)
    X = getBSseq(BSobj, "M");  N = getBSseq(BSobj, "Cov")
    ## locate the data for plotting
    chr = as.character(OneDMR$chr)
    ix.chr = which(allchr==chr)
    thispos = allpos[ix.chr]
    thisN = N[ix.chr,]
    thisX = X[ix.chr,]
    xlim = c(OneDMR$start-ext, OneDMR$end+ext)
    ix1 = which(thispos<=xlim[2] & thispos>=xlim[1])
    ## plot it
    nSample = ncol(X)
    if(nSample > 2) {
        y.cex=0.66
    } else  y.cex=1

    sNames = sampleNames(BSobj)
    par(mfrow=c(nSample, 1), mar = c(2.5, 2.5, 1.6, 2.5), mgp = c(1.5, 0.5, 0))
    thisP = thisX/thisN

    for(i  in 1:ncol(X)) {
        ## blue foreground, methylated, normalized to ylim
        ## plot(thispos[ix1], thisP[ix1,i], type="h", col="blue", lwd=1, xlab=chr, ylim=ylim, xlim=xlim,
        ## ylab="methyl%", main=sNames[i], col.lab="blue")
        plot(thispos[ix1], thisP[ix1,i], type="h", col="blue", axes=F, lwd=1.5,
             xlab='', ylab='', ylim=ylim, xlim=xlim,
             main=sNames[i])
        box(col="black")
        axis(1,)
        axis(2, col="blue", col.axis="blue")
        mtext(chr, side=1, line=1.33, cex=y.cex)
        mtext("methyl%", side=2, line=1.33, col="blue", cex=y.cex)
        ## black curve, total number, normalized to ylim
        thisN.norm = thisN[ix1,i]/max(thisN[ix1,])*ylim[2]
        lines(thispos[ix1], thisN.norm, type="l", col="gray", lwd=1.5)
        axis(side=4, at=seq(0,ylim[2],length.out=5),
             labels=round(seq(0, max(thisN[ix1,]), length.out=5)) )
        mtext("read depth", side=4, line=1.33, cex=y.cex)
        ## plot a shaded region for DMR
        rect(OneDMR$start, ylim[1], OneDMR$end, ylim[2], col="#FF00001A", border = NA)
    }
}

