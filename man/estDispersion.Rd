\name{estDispersion}
\alias{estDispersion,SeqCountSet-method}
\alias{estDispersion}

\title{
  Estimate and shrink tag-specific dipsersions 
}
\description{
  This function first estimate tag-specific dipsersions using a method
  of moment estimator. Then the dipsersions are shrunk based a penalized
  likelihood approach. 
}
\usage{
\S4method{estDispersion}{SeqCountSet}(seqData, trend=FALSE)
}
\arguments{
  \item{seqData}{An object of SeqCountSet class.}
  \item{trend}{A binary indicator for modeling the dispersion~expression
  trend.} 
}
\details{
  The function takes and object of seqCountData class and return the
  same oject with ``dispersion'' field filled.

  With ``trend=TRUE'' the dependence of dispersion on mean expressions
  will be modeled. In that case the shrinkage will be performed
  conditional on mean expressions. 
}

\examples{
data(seqData)
seqData=estNormFactors(seqData)
seqData=estDispersion(seqData)
head(dispersion(seqData))
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\keyword{RNA-seq}
