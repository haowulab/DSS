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
DMLtest.multiFactor(DMLfit, coef = 2, term, Contrast)
}

\arguments{
  \item{DMLfit}{Result object returned from 'DMLfit.multiFactor' function.}
  \item{coef}{It can be an integer to indicate which coefficient in the linear model
	is be tested for being zero. Be careful of intercept. If the model
	contains intercept, coef=2 indicate testing the first factor in the
	formula. If the model has no intercept, testing first factor should
	use coef=1.

	It can also be a character for the terms to be tested. In that case it
	must match one of the column names in the design matrix. One can
	look at colnames(DMLfit$X) to obtain the column names. 
  }
  \item{term}{The term(s) to be tested, as one or a vector of
	characters. Can be multiple terms. See 'Hypothesis test' section for
	details.} 
  \item{Contrast}{A contrast matrix for hypothesis testing. The number
	of rows must equal to the number of columns in the design matrix
	ncol(DMLfit$X). See 'Hypothesis test' section for details.}																
}

\section{Hypothesis test}{
  User can specify one of the following parameter for testing:
  'coef', 'term', or 'Contrast'.

  When specifying 'coef', it tests *one* parameter in the model, which
  corresponds to one column in the design matrix. In this case,   
  A Wald test is performed using the estimated coefficient and standard
  error from 'DMLfit.multiFactor'. P-values are derived from test
  statistics based on normal distribution.

  When specifying 'term', it tests the whole term in the model. If the
  term is continuous or a categorical variable with only two levels (one
  degree of freedom), it is equivalent to specifying 'coef' because it
  tests only one parameter in the model. However, when the term is a
  categorical variable with more than two levels, it will test multiple
  parameters at the same time so it's a compound hypothesis test, and a
  F-test will be performed. 
  
  Specifying 'Contrast' matrix provides the most flexible test
  procedure. It can test any linear combination of the parameter, and
  F-test will be performed. Let L be the contrast matrix. The hypothesis
  test performed is H0: L^T * beta = 0. Thus the number of rows in L must
  equal to the number of items of beta (which is the number of columns
  in the design matrix). 

  Using 'term' or 'Contrast' will be slower especially when there are a
  lot of CpG sites, because the computation cannot be vectorized (each
  CpG site has a different variance/covariance matrix for the estimated
  coefficients).
  
  FDR is computed using cannonical Benjamini-Hochberg procedure.
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
\dontrun{
data(RRBS)
## model fitting
DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell+case:cell)

## hypothesis testing - following two lines do the same thing
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=3)
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="cellrN")

## this doesn't work
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="cell")

## look at distributions of test statistics and p-values
par(mfrow=c(1,2))
hist(DMLtest.cell$stat, 100, main="test statistics")
hist(DMLtest.cell$pvals, 100, main="P values")

## Using term or Contrast
DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell)

## following 4 tests should produce the same results,
## since 'case' only has two levels.
## However the p-values from F-tests (using term or Contrast) are
## slightly different, due to normal approximation in Wald test.
test1 = DMLtest.multiFactor(DMLfit, coef=2)
test2 = DMLtest.multiFactor(DMLfit, coef="caseSLE")
test3 = DMLtest.multiFactor(DMLfit, term="case")
Contrast = matrix(c(0,1,0,0), ncol=1)
test4 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)
cor(cbind(test1$pval, test2$pval, test3$pval, test4$pval))


## note the different usage of term and coef.
## 'term' has to be in the formula, whereas 'coef' has to be in colnames
## of the design matrix.
DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell+case:cell)
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="cellrN")
DMLtest.cell = DMLtest.multiFactor(DMLfit, term="cell")
DMLtest.int = DMLtest.multiFactor(DMLfit, coef="caseSLE:cellrN")
DMLtest.int = DMLtest.multiFactor(DMLfit, term="case:cell")

}
}


