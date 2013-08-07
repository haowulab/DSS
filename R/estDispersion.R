### functions to estimate and shrink gene specific over dispersion
## this works for single factor design only at this time!!!
estDispersion <- function(seqData, trend=FALSE) {
  if(!is(seqData, "SeqCountSet"))
    stop("Input must be an object of SeqCountSet class!")
  
  ## design matrix
  design <- pData(seqData)
  if(ncol(design) == 1) { ## single factor design
    if(trend) 
      seqData=estDispersion.trend(seqData)
    else
      seqData=estDispersion.notrend(seqData)
  } else { ## multiple factor
    seqData=estDispersion.multiFactor(seqData)
  }
  
  seqData
}

