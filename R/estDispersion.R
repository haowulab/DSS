### functions to estimate and shrink gene specific over dispersion
## this works for single factor design only at this time!!!
estDispersion <- function(seqData, trend=FALSE) {
  if(!is(seqData, "SeqCountSet"))
    stop("Input must be an object of SeqCountSet class!")
  if(trend) 
    seqData=estDispersion.trend(seqData)
  else
    seqData=estDispersion.notrend(seqData)
  seqData
}

