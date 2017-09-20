\name{DMLfit.multiFactor}
\alias{DMLfit.multiFactor}

\title{
  Fit a linear model for BS-seq data from general experimental design
}

\description{
  This function takes a BSseq object, a data frame for experimental
  design and a model formula and then fit a linear model. 
}

\usage{
DMLfit.multiFactor(BSobj, design, formula, smoothing=FALSE, smoothing.span=500)
}

\arguments{
  \item{BSobj}{An object of BSseq class for the BS-seq data.}
  \item{design}{A data frame for experimental design. Number of rows
	must match the number of columns of the counts in BSobj.}
  \item{formula}{A formula for the linear model.}
  \item{smoothing}{A flag to indicate whether to apply smoothing. When
	true, the counts will be smoothed by a simple moving average method.}
  \item{smoothing.span}{The size of smoothing window, in
	basepairs. Default is 500.} 
}

\details{
  The lineear model fitting is done through ordinary least square on the
  arscine transformed methylation percentages. The estimated standard
  errors are computed with consideration of the data (count)
  distribution and transformation. This function is extremely
  efficient. The computation takes around 20 minutes for 4 million CpG
  sites.
}

\value{
  A list with following components
  \item{gr}{An object of 'GRanges' for locations of the CpG sites.}
  \item{design}{The input data frame for experimental design.}
  \item{formula}{The input formula for the model.}
  \item{X}{The design matrix used in regression. It is created based on
  design and formula.}
  \item{fit}{The model fitting results. This is a list itself, with two
  components: 'beta' - the estimated coefficients; 'var.beta' - estimated
  variance/covariance matrices for beta. Note that var.beta for a CpG
  site should be a ncol(X) x ncol(X) matrix, but is flattend to a vector
  so that the matrices for all CpG sites can be saved as a matrix.}
}

\author{
  Hao Wu<hao.wu@emory.edu>
}

\seealso{
  DMLtest.multiFactor, DMLtest
}

\examples{
data(RRBS)
## model fitting
DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell+case:cell)

## with smoothing:
DMLfit.sm = DMLfit.multiFactor(RRBS, design, ~case+cell+case:cell, smoothing=TRUE)

## hypothesis testing
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=3)

## look at distributions of test statistics and p-values
par(mfrow=c(1,2))
hist(DMLtest.cell$stat, 100, main="test statistics")
hist(DMLtest.cell$pvals, 100, main="P values")
}

