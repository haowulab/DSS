\name{callDML}
\alias{callDML}

\title{
  Function to detect differntially methylated loci (DML) from bisulfite
  sequencing (BS-seq) data. 
}
\description{
  This function takes the results from DML testing procedure ('DMLtest'
  function) and calls DMLs. Regions will CpG sites being statistically
  significant are deemed as DMLs. 
}
\usage{
callDML(DMLresult, delta=0.1, p.threshold=1e-5)
}

\arguments{
  \item{DMLresult}{
    A data frame representing the results for DML detection. This should
    be the result returned from 'DMLtest' function.}
  \item{delta}{A threshold for defining DML. In DML testing procedure, 
    hypothesis test that the two groups means are equalis is conducted at each CpG
    site. Here if 'delta' is specified, the function will compute
    the posterior probability that the difference of the means are
    greater than delta,and then call DML based on that. This only works
    when the test results are from 'DMLtest', which is for two-group
    comparison. For general design, this has to be set to 0.}
  \item{p.threshold}{
    When delta is not specified, this is the threshold of p-values for
	defining DML, e.g. Loci with p-values less than this threshold will
	be deemed DMLs. When delta is specified, CpG sites with posterior
	probability greater than 1-p.threshold are deemed DML.
  }
}

\value{A data frame for DMLs. Each row is for a DML. DMLs are sorted by
  statistical significance. The columns are
  \item{chr}{Chromosome number.}
  \item{pos}{Genomic coordinates.}
  \item{mu1, mu2}{Mean methylations of two groups.}
  \item{diff}{Difference of mean methylations of two groups.}
  \item{diff.se}{Standard error of the methylation difference.}
  \item{stat}{Wald statistics.}
  \item{phi1, phi2}{Estimated dispersions in two groups.}
  \item{pval}{P-values. This is obtained from normal distribution.}
  \item{fdr}{False discovery rate.}
  \item{postprob.overThreshold}{The posterior probability of the
	difference in methylation greater than delta. This columns is only
	available when delta>0.}
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
DMLtest, callDMR
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
BSobj <- makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2),
  c("C1","C2", "N1", "N2") )

##  DML test
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1","N2"))

## call DML
dmls <- callDML(dmlTest)
head(dmls)

## call DML with a threshold
dmls2 <- callDML(dmlTest, delta=0.1)
head(dmls2)

## For whole-genome BS-seq data, perform DML test with smoothing
require(bsseqData)
data(BS.cancer.ex)
## takea smallportionof data and test
BSobj <- BS.cancer.ex[10000:15000,]
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"),
   smoothing=TRUE, smoothing.span=500)
dmls <- callDML(dmlTest)
head(dmls)

## from multifactor design
data(RRBS)
DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell+case:cell)
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="cellrN")
dml = callDML(DMLtest.cell, p.threshold=0.05) ## this produce a warning
dml = callDML(DMLtest.cell, delta=0, p.threshold=0.05) ## no warning

head(dml)

}
}
