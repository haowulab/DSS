\name{estDispersion}
\alias{estDispersion,SeqCountSet-method}
\alias{estDispersion}

\title{
  Estimate and shrink tag-specific dipsersions 
}
\description{
  This function first estimate tag-specific dipsersions using a method
  of moment estimator. Then the dipsersions are shrunk based a penalized
  likelihood approach. The function works for general experimental
  designs. 
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

  The function works for multiple factor designs. But option
  ``trend=TRUE'' only applicable for single factor experiment. 
}

\examples{
data(seqData)
seqData=estNormFactors(seqData)
seqData=estDispersion(seqData)
head(dispersion(seqData))

## For multiple factor design
data(seqData)
Y=exprs(seqData)
design=data.frame(gender=c(rep("M",4), rep("F",4)), strain=rep(c("WT", "Mutant"),4))
X=as.data.frame(model.matrix(~gender+strain, data=design))
seqData=newSeqCountSet(Y, X)
seqData=estDispersion(seqData)
head(dispersion(seqData))

## the hypothesis testing for multifactor experiments can be performed
## using edgeR function, with DSS estimated dispersions
\dontrun{
library(edgeR)
fit.edgeR <- glmFit(Y, X, lib.size=normalizationFactor(seqData), dispersion=dispersion(seqData))
lrt.edgeR <- glmLRT(fit.edgeR, coef=2)
head(lrt.edgeR$table)
}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\keyword{RNA-seq}
