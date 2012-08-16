setClass( "SeqCountSet",
         contains = "ExpressionSet",
         representation = representation(
           normalizationFactor="numeric",
           dispersion="numeric"
           )
         )

newSeqCountSet <- function(counts, designs, normalizationFactor=NULL,featureData = NULL ) {
  ## some error checking stuff

  if( is.null( featureData ) )
    featureData <- annotatedDataFrameFrom( counts, byrow=TRUE )

  ## check input normalizationFactor
  if(is.null(normalizationFactor))
    normalizationFactor <- rep(1, ncol(counts))
  else { # with normalizationFactor provided
    if(is(normalizationFactor, "numeric")) { ## input is vector
      if(length(normalizationFactor) != ncol(counts))
        stop("Input normalizationFactor has wrong length.")
    } else if(is.matrix(normalizationFactor)) { ## matrix normalization factor
      ## check dimension
      if(any(dim(counts) != dim(normalizationFactor)))
        stop("Input normalizationFactor has wrong dimension.")
      ## convert to vector
      normalizationFactor=as.numeric(normalizationFactor)
    } else { ## other things
      stop("Wrong input normalizationFactor, must be a numeric vector or matrix.")
    }
  }
  
  ## work on input design matrix. Carefully about the input design.
  if(is(designs, "vector") ) { # single factor design
    designs <- as.data.frame(designs )
    colnames(designs)="designs"
  } else  if( is( designs, "matrix" ) ) { ## multiple factor design
    if(ncol(designs) == 1) { ## still single factor
      designs <- as.data.frame(designs )
      colnames(designs)="designs"
    } else { ## input is a data frame
      designs <- as.data.frame(designs )
    }
  }
  ## if there's no row names in designs, assume it has the same sequence as the columns in count matrix
  if( is.null(rownames(designs)) & !is.null(colnames(counts)) ) {
    rownames(designs)=colnames(counts)
  }

  if( is( designs, "data.frame" ) || is( designs, "AnnotatedDataFrame" ) ) {
    stopifnot( nrow( designs ) == ncol( counts ) )
    designs <- as( designs, "AnnotatedDataFrame" )
  }

  
  res <- new("SeqCountSet",
             exprs = counts, 
             phenoData = designs,
             featureData = featureData,
             normalizationFactor=normalizationFactor,
             dispersion=numeric(0)
             )
  res

}

### normalization factor
setGeneric("normalizationFactor", function(object) standardGeneric("normalizationFactor"))
setGeneric("normalizationFactor<-", function(object, value) standardGeneric("normalizationFactor<-"))

setMethod("normalizationFactor", signature(object="SeqCountSet"),
          function(object) {
            f=object@normalizationFactor
            n=ncol(exprs(object))
            if(length(f)==n) { ## size factor as a vector
              res=f
            } else { ## a matrix
              res=matrix(f, ncol=n)
            }
            res
          })

setReplaceMethod("normalizationFactor", signature(object="SeqCountSet", value="numeric"), 
                 function(object, value ) {
                   if(length(value) != ncol(exprs(object)))
                     stop("wrong length for normalization factor vector!")
                   object@normalizationFactor <- value
                   object
                 })

setReplaceMethod("normalizationFactor", signature(object="SeqCountSet", value="matrix"), 
                 function(object, value ) {
                   if(any(dim(value) != dim(exprs(object))))
                     stop("wrong dimension for normalization factor matrix!")
                   object@normalizationFactor <- as.numeric(value)
                   object
                 })

## for gene specific over-dispersion
setGeneric("dispersion", function(object) standardGeneric("dispersion"))
setGeneric("dispersion<-", function(object, value) standardGeneric("dispersion<-"))

setMethod("dispersion", signature(object="SeqCountSet"),
          function(object) {
            object@dispersion
          })

setReplaceMethod("dispersion", signature(object="SeqCountSet", value="numeric"), 
                 function(object, value ) {
                   if(length(value) != nrow(exprs(object)))
                     stop("wrong length for dispersion vector!")
                   object@dispersion <- value
                   object
                 })

