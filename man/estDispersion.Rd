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
\S4method{estDispersion}{SeqCountSet}(seqData)
}
\arguments{
  \item{seqData}{An object of SeqCountSet class.}
}
\details{
  The function takes and object of seqCountData class and return the
  same oject with ``dispersion'' field filled.
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
