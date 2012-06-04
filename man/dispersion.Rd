\name{dispersion}
\alias{dispersion,SeqCountSet-method}
\alias{dispersion<-,SeqCountSet,numeric-method}
\alias{dispersion}
\alias{dispersion<-}
\title{
  Acessor functions for the 'dispersion' slot in a
  SeqCountData object.
}
\description{
  Dispersion parameter for a gene represents its coefficient of
  variation of expressions. It characterizes the biological
  variations. 
}
\usage{
\S4method{dispersion}{SeqCountSet}(object)
\S4method{dispersion}{SeqCountSet,numeric}(object) <- value
}
\arguments{
  \item{object}{A SeqCountData object.}
  \item{value}{A numeric vector with the same length as number of genes.}
}

\details{
  If the counts from biological replicates are modeled as negative
  binomial distribution, the variance (v) and mean (m) should hold
  following relationship: v=m+m^2*phi, where phi is the
  dispersion. Another interpretation is that phi represents the
  biological variations among replicates when underlying expressions are
  modeled as a Gamma distribution.
}
\author{
  Hao Wu <hao.wu@emory.edu>
}
\seealso{
  normalizationFactor
}
\examples{
data(seqData)
## obtain
seqData=estNormFactors(seqData, "quantile")
seqData=estDispersion(seqData)
dispersion(seqData)

## assign
dispersion(seqData)=rep(0.1, nrow(exprs(seqData)))

}
