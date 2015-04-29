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
DMLfit.multiFactor(BSobj, design, formula)
}

\arguments{
  \item{BSobj}{An object of BSseq class for the BS-seq data.}
  \item{design}{A data frame for experimental design. Number of rows
  must match the number of columns of the counts in BSobj.}
  \item{formula}{A formula for the linear model.}
}

\details{
  The line model fitting is done through ordinary least square on the
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
  \item{fit}{The model fitting results. This is a list itself, with two
  components: 'beta' - the estimated coefficients; 'se.beta' - estimated
  standard error for beta.}
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

## hypothesis testing
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=3)

## look at distributions of test statistics and p-values
par(mfrow=c(1,2))
hist(DMLtest.cell$stat, 100, main="test statistics")
hist(DMLtest.cell$pvals, 100, main="P values")
}

