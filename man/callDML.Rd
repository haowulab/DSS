\name{callDML}
\alias{callDML}
\title{
  Function to detect differntially methylated loci (DML) for two group
  comparisons of bisulfite sequencing (BS-seq) data.
}
\description{
  The replicated BS-seq data are modeled as Beta-Binomial
  distribution. Similar to that in Gamma-Poisson distributions, the
  biological variations are captured by the dispersion parameter. This
  function takes BS-seq data from two conditions, and estimate the
  dispersion parameters through a Bayesian hierarchical model. Then a
  Wald test is constructed to test the differntial methylation.
}
\usage{
callDML(BS1, BS2, equal.disp = FALSE, threshold = 0)
}

\arguments{
  \item{BS1}{An object of BSseq class for the BS-seq data from the first
  condition.} 
\item{BS2}{Anobject of BSseqclass for the BS-seq data from the second
  condition.}
  \item{equal.disp}{A flag to indicate whether the dispersion in two
  groups are deemed equal. Default is FALSE, and the dispersion
  shrinkages are performed on two conditions independently.}
  \item{threshold}{A threshold for defining DML. See details for more description.}
}

\details{
  By default, a hypothesis testing that the two groups means are equal
  is conducted at each CpG site. The spatial correlations among the CpG
  sites are ignored in this function, and the tests are conducted
  independently at each CpG site. If 'threshold' is specified, the
  function will also report the posterior probability that the
  difference of the means are greater than the threshold.

  Due to the differences in coverages, some CpG sites are not covered
  at all replicates. For a CpG sites, if it's only covered by one
  replicate, the test cannot be performed because the within group
  variances cannot be estimated. For those loci the results will be NA.
}

\value{A data frame with each row corresponding to a CpG site. Rows are
  sorted by chromosome number and genomic coordinates. The columns
  include:
  \item{chr}{Chromosome number.}
  \item{pos}{Genomic coordinates.}
  \item{mu1, mu2}{Mean methylations of two groups.}
  \item{diff}{Difference of mean methylations of two groups.}
  \item{diff.se}{Standard error of the methylation difference.}
  \item{stat}{Wald statistics.}
  \item{pval}{P-values. This is obtained from normal distribution.}
  \item{fdr}{False discovery rate.}
}

\author{
  Hao Wu <hao.wu@emory.edu>}

\seealso{
makeBSseqData
}

\examples{
\dontrun{
require(bsseq)

## first read in methylation data.
path <- file.path(system.file(package="DSS"), "extdata")
dat1.1 <- read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 <- read.table(file.path(path, "cond1_2.txt"), header=TRUE)
dat2.1 <- read.table(file.path(path, "cond2_1.txt"), header=TRUE)
dat2.2 <- read.table(file.path(path, "cond2_2.txt"), header=TRUE)

## make BSseq objects
BS1 <- makeBSseqData( list(dat1.1, dat1.2), paste("cond1",1:2,sep=".") )
BS2 <- makeBSseqData( list(dat2.1, dat2.2), paste("cond2",1:2,sep=".") )

## call DML 
result <- callDML(BS1, BS2)
head(result)
}
}
