\name{DSS.DE}
\alias{DSS.DE}
\title{
  Perform RNA-seq differential expression analysis in two-group comparison
}
\description{
  This is the top level wrapper function for RNA-seq differential expression
  analysis in a two-group comparison. Users only need to provide the
  count matrix and a vector for design, and obtain DE test results. 
}
\usage{
DSS.DE(counts, design)
}
\arguments{
  \item{counts}{A matrix of integers with rows corresponding to
    genes and columns for samples.}
  \item{designs}{A vector representing the treatment groups. It must be
	a vector of 0 and 1. The length of the vector must match the number
	of columns of input count matrix.}
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
counts = matrix(rpois(600, 10), ncol=6)
designs = c(0,0,0,1,1,1)
result = DSS.DE(counts, designs)
head(result)
}


