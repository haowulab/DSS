\name{estNormFactors}
\alias{estNormFactors,SeqCountSet-method}
\alias{estNormFactors}

\title{
  Estimate normalization factors
}
\description{
  This function estimates normalization factors for the input
  'seqCountSet' object and return the same object with
  normalizationFactor field filled or replaced.
}
\usage{
\S4method{estNormFactors}{SeqCountSet}(seqData, method=c("lr", "quantile", "total", "median"))
}

\arguments{
  \item{seqData}{An object of "SeqCountSet" class. }
  \item{method}{Methods to be used in computing normalization
    factors. Currently available options only include methods to compute
    normalization factor to adjust for sequencing depths.
    Available options use (1) "lr": using median of logratio of
	counts. Similar to the TMM method. (2) "quantile"
    (default):  75th quantile, (3) "total":  total counts, or (4)
    "median":  median counts to constuct the normalization
    factors. From all methods the normalization factor will be a vector
    with same length as number of columns for input counts.}
}

\value{
  The same "SeqCountSet" object with normalizationFactor field filled or
  replaced. 
}

\examples{
data(seqData)
## compare different methods
seqData=estNormFactors(seqData, "lr")
k1=normalizationFactor(seqData)
seqData=estNormFactors(seqData, "quantile")
k2=normalizationFactor(seqData)
seqData=estNormFactors(seqData, "total")
k3=normalizationFactor(seqData)
cor(cbind(k1,k2,k3))

## assign size factor
normalizationFactor(seqData)=k1

## or normalization factor can be a matrix
dd=exprs(seqData)
f=matrix(runif(length(dd), 1,10), nrow=nrow(dd), ncol=ncol(dd))
normalizationFactor(seqData)=f
head(normalizationFactor(seqData))
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\keyword{normalization}
