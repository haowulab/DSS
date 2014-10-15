################################################
## A wrapper function for performing RNA-seq DE.
## Use all default parameter values.
################################################
DSS.DE <- function(counts, design) {
    ## some error checkings
    if(! all(unique(design) %in% c(0,1)) )
        stop("design must be a vector of 0 and 1.")
    if(length(design) != ncol(counts))
        stop("length of design doesn't match the number of columns of count matrix.")

    seqData = newSeqCountSet(counts, design)
    seqData = estNormFactors(seqData)
    seqData = estDispersion(seqData)
    result = waldTest(seqData, 0, 1)
    result
}


