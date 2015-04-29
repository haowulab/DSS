\name{DMLtest.multiFactor}
\alias{DMLtest.multiFactor}

\title{
  Perform statistical test for BS-seq data from general experimental design
}

\description{
  This function takes the linar model fitting results and performs Wald
  test at each CpG site, then return test statistics, p-values and FDR.
}

\usage{
DMLtest.multiFactor(DMLfit, coef = 2)
}

\arguments{
  \item{DMLfit}{Result object returned from 'DMLfit.multiFactor' function.}
  \item{coef}{Integer to indicate which coefficient in the linear model
	is be tested for being zero. Be careful of intercept. If the model
	contains intercept, coef=2 indicate testing the first factor in the
	formula. If the model has no intercept, testing first factor should
	use coef=1. }
}

\details{
  A Wald test is performed using the estimated coefficient and standard
  error from 'DMLfit.multiFactor'. P-values are derived from test
  statistics based on normal distribution. FDR is computed using cannonical
  Benjamini-Hochberg procedure. 
}

\value{
  A data frame with following columns: chr, pos, stat, pvals, fdr. Each
  row is for a CpG site. Note that the CpG sites are sorted by
  chromosome and position. 
}

\author{
  Hao Wu<hao.wu@emory.edu>
}

\seealso{
  DMLfit.multiFactor, DMLtest
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


