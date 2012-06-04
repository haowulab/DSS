\name{normalizationFactor}
\alias{normalizationFactor,SeqCountSet-method}
\alias{normalizationFactor<-,SeqCountSet,numeric-method}
\alias{normalizationFactor<-,SeqCountSet,matrix-method}
\alias{normalizationFactor}
\alias{normalizationFactor<-}

\title{
  Accessor functions for the 'normalizationFactor' slot in a
  SeqCountData object.
  
}
\description{
  The normalization factors are used to adjust for technical or
  biological biases in the sequencing experiments. The factors 
  can either be (1) a vector with length equals
  to the number of columns of the count data; or (2) a matrix with the same
  dimension of the count data.
}
\usage{
\S4method{normalizationFactor}{SeqCountSet}(object)
\S4method{normalizationFactor}{SeqCountSet,numeric}(object) <- value
\S4method{normalizationFactor}{SeqCountSet,matrix}(object) <- value
}

\arguments{
  \item{object}{A SeqCountData object.}
  \item{value}{A numeric vector or matrix. If it is a vector it must
    have length equals to the number of columns of the count data. For
    matrix it must have the same dimension of the count data.}
}

\details{
  The vector normalization factors are used mostly to
  correct for sequencing depth from different datasets. The matrix
  factor applies a different normalizing constant for each gene at each
  sample to adjust for a broader range of artifacts such as GC content. 
}
  
\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
  dispersion
}
\examples{
data(seqData)
## obtain nomalization factor
seqData=estNormFactors(seqData, "quantile")
normalizationFactor(seqData)

## assign as vector
normalizationFactor(seqData)=rep(1, ncol(exprs(seqData))) ## getan error here

## or assign as a matrix
f=matrix(1, nrow=nrow(exprs(seqData)), ncol=ncol(exprs(seqData)))
normalizationFactor(seqData)=f
}

