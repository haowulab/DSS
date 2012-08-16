## perform wald test
## This works for two group comparison only!
waldTest <- function(seqData, sampleA, sampleB, equal.var=FALSE) {
  ## compute two group means
  Y=exprs(seqData)+0.1
  design=pData(phenoData(seqData))$designs
  k=normalizationFactor(seqData)
  if(is.matrix(k))
    Y2=Y/k
  else
    Y2=sweep(Y, 2, k, FUN="/")

  ## difference
  idxA=design==sampleA
  nA=sum(idxA)
  muA=rowMeans(Y2[,idxA])
  idxB=design==sampleB
  nB=sum(idxB)
  muB=rowMeans(Y2[,idxB])
  d=muA-muB
  ## variance
  phi.hat=dispersion(seqData)
  ## use common mu:
  mu0=rowMeans(Y2)
  n=length(k)
  if(equal.var) {
    std=sqrt((mu0*sum(1/k[idxA])+ nA*mu0^2*phi.hat)/nA^2+
      (mu0*sum(1/k[idxB])+nB*mu0^2*phi.hat)/nB^2)
  } else  {
    std=sqrt((muA*sum(1/k[idxA])+ nA*muA^2*phi.hat)/nA^2+
      (muB*sum(1/k[idxB])+nB*muB^2*phi.hat)/nB^2)
  }
  
  stat=d / std
  pval=2*(1-pnorm(abs(stat)))

  ## FDR
  normstat = (stat - median(stat)) / (IQR(stat) / 1.349)
  fdrres = locfdr(normstat, plot=0)
  fdr = fdrres$fdr
  ix = sort(abs(stat), decreasing = TRUE, index.return = TRUE)$ix
  tmp = cumsum(fdr[ix])/(1:length(fdr))
  fdr.avg = tmp
  fdr.avg[ix] = tmp
  fdr.global = numeric(length(normstat))
  for (i in 1:length(normstat)){
    ind = (abs(fdrres$mat[,"x"]) >= abs(normstat[i]))
    F0l = sum(fdrres$mat[ind,"f0"] * fdrres$fp0["mlest", 3])
    Fl = sum(fdrres$mat[ind,"f"])
    fdr.global[i] = F0l / Fl
  }
  ## some times (e.g., when there's no true DE) localfdr isn't stable and will
  ## generate p-values bigger than one. Fix it.
  fdr.global[fdr.global>1]=1
  
  ## generate a data frame of reports
  lfc=log((muA+0.5)/(muB+0.5))
  difExpr=muA-muB
  genes=featureData(seqData)
  if(!is.null(genes))
    genes=as(genes, "data.frame")
  result=data.frame(geneIndex=1:nrow(Y), muA=muA, muB=muB, lfc=lfc, difExpr=difExpr, stats=stat,
    pval=pval, local.fdr=fdr, fdr=fdr.global, genes)
  ix=sort(abs(result[,"stats"] ), decreasing=TRUE, index.return=TRUE)$ix
  result=result[ix,]
  ## fix the global FDR for first gene
  result[1,"fdr"]=result[1,"local.fdr"]
  result
}
