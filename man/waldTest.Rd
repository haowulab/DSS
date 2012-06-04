\name{waldTest}
\alias{waldTest,SeqCountSet-method}
\alias{waldTest}
\title{
  Perform gene-wise Wald test for two group comparisons for sequencing
  count data. 
}
\description{
  The counts from two groups are modeled as negative binomial random
  variables with means and dispersions estimated. Wald statistics will
  be constructed. P-values will be obtained based on Gaussian assumption.
}
\usage{
\S4method{waldTest}{SeqCountSet}(seqData, sampleA, sampleB, equal.var)
}
\arguments{
  \item{seqData}{An object of SeqCountSet class.}
  \item{sampleA}{The sample labels for the first sample to be
    compared in two-group comparison.}
  \item{sampleB}{The sample labels for the second sample to be
    compared in two-group comparison.}
  \item{equal.var}{A boolean to indicate whether to use the same or
    different means in two groups for computing variances in Wald
    test. Default is FALSE.}
}

\details{
  The input seqCountData object Must have normalizationFactor and
  dispersion fields filled, e.g., estNormFactors and estDispersion need
  to be called prior to this. With group means and shrunk dispersions
  ready, the variances for difference in group means will be constructed
  based on Negative Binomial distribution. P-values will be obtained
  under the assumption that the Wald test statistics are normally
  distributed. Genes with 0 counts in both groups will be assigned 0 for
  test statistics and 1 for p-values.
}

\value{
  A data frame with each row corresponding to a gene. Rows are sorted
  according to wald test statistics. The columns are:
  \item{gene Index}{index for input gene orders, integers from 1 to the
    number of genes.}
  \item{muA}{sample mean (after normalization) for sample A.}
  \item{muB}{sample mean (after normalization) for sample B.}
  \item{lfc}{log fold change of expressions between two groups.}
  \item{difExpr}{differences in expressions between two groups.}
  \item{stats}{Wald test statistics.}
  \item{pval}{p-values.}
  \item{others}{input gene annotations supplied as AnnotatedDataFrame
    when constructed the SeqCountData object.}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\examples{
data(seqData)
seqData=estNormFactors(seqData)
seqData=estDispersion(seqData)
result=waldTest(seqData, 0, 1)
head(result)
}


