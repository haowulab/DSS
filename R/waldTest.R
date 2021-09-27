## perform wald test
## This works for two group comparison only!
waldTest <- function(seqData, sampleA, sampleB, equal.var=FALSE,
                     fdr.method=c("BH", "locfdr"))  {

    fdr.method = match.arg(fdr.method)
    Y0 = exprs(seqData)
    idx1 = rowSums(Y0)>0 ## genes with some counts

    ## compute two group means
    Y = exprs(seqData) ## +0.1
    design = pData(phenoData(seqData))$designs
    k = normalizationFactor(seqData)
    if(is.matrix(k))
        Y2 = Y/k
    else
        Y2 = sweep(Y, 2, k, FUN="/")

    ## difference
    idxA = design==sampleA
    nA = sum(idxA)
    muA = rowMeans(Y2[,idxA,drop=FALSE])
    idxB = design==sampleB
    nB = sum(idxB)
    muB = rowMeans(Y2[,idxB, drop=FALSE])
    d = muA-muB
    ## variance
    phi.hat = dispersion(seqData)
    ## use common mu:
    mu0 = rowMeans(Y2)
    n = length(k)
    if(equal.var) {
        std = sqrt((mu0*sum(1/k[idxA])+ nA*mu0^2*phi.hat)/nA^2+
                     (mu0*sum(1/k[idxB])+nB*mu0^2*phi.hat)/nB^2)
    } else  {
        std = sqrt((muA*sum(1/k[idxA])+ nA*muA^2*phi.hat)/nA^2+ 
                     (muB*sum(1/k[idxB])+nB*muB^2*phi.hat)/nB^2)
    }

    stat = d / std
    stat[is.na(stat)] = 0
    pval = 2*(1-pnorm(abs(stat)))

    ## generate a data frame of reports
    lfc = log((muA+0.5)/(muB+0.5)) ## log fold change
    difExpr = muA-muB ## difference in expressions
    genes = featureData(seqData)
    if(!is.null(genes))
      genes=as(genes, "data.frame")
    result = data.frame(geneIndex=1:nrow(Y), muA=muA, muB=muB, 
                        lfc=lfc, difExpr=difExpr, stats=stat, 
                        pval=pval, genes)
    
    if(fdr.method == "BH") {
      fdr0 = p.adjust(pval[idx1], method="BH")
      fdr = pval
      fdr[idx1] = fdr0
      result = data.frame(result, fdr=fdr)
    }
    
    if(fdr.method == "locfdr") {
        ## FDR esimation using locfdr. This is only for genes with non-zero counts.
        ## Need to work on this more, since it's not very good.
        stat1 = stat[idx1]
        normstat = (stat1 - median(stat1)) / (IQR(stat1) / 1.349)
        fdrres = locfdr(normstat, plot=0)
        fdr.global = numeric(length(normstat))
        xx = fdrres$mat[, "x"]
        leftbreaks = xx - (xx[2]-xx[1])/2; leftbreaks[1] = min(normstat)
        rightbreaks = xx + (xx[2]-xx[1])/2; rightbreaks[length(rightbreaks)] = max(normstat)
        for (i in 1:length(normstat)) {
            ind = ((leftbreaks <= (-1)*abs(normstat[i])) | (rightbreaks >= abs(normstat[i])))
            F1l = sum(fdrres$mat[ind,"p1f1"])
            Fl = sum(fdrres$mat[ind, "f"])
            fdr.global[i] = 1 - F1l/Fl
        }

        ## go back to all genes. Genes with all 0 counts will have local and global FDR NA.
        fdr0 = rep(NA, nrow(Y))
        fdr.global0 = rep(NA, nrow(Y))
        fdr0[idx1] = fdrres$fdr
        fdr.global0[idx1] = fdr.global
        result = data.frame(result, local.fdr=fdr0, fdr=fdr.global0)
        ## fix the global FDR for first gene
        result[1,"fdr"] = result[1,"local.fdr"]
    }
    
    ## sort the results by test statistics
    ix = sort(abs(result[,"stats"] ), decreasing=TRUE, index.return=TRUE)$ix
    result = result[ix,]

    return(result)
}
