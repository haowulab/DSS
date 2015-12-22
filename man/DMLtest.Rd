\name{DMLtest}
\alias{DMLtest}
\title{
  Function to perform statistical test of differntially methylated loci
  (DML) for two group comparisons of bisulfite sequencing (BS-seq) data.
}

\description{
  This function takes a BSseq object and two group labels, then perform
  statistical tests for differntial methylation at each CpG site.
}

\usage{
DMLtest(BSobj, group1, group2, equal.disp = FALSE, smoothing = FALSE, smoothing.span = 500) 
}

\arguments{
  \item{BSobj}{An object of BSseq class for the BS-seq data.}
  \item{group1, group2}{Vectors of sample names or indexes for the two
    groups to be tested. See more description in details.}
  \item{equal.disp}{A flag to indicate whether the dispersion in two
    groups are deemed equal. Default is FALSE, and the dispersion
    shrinkages are performed on two conditions independently.}
  \item{smoothing}{A flag to indicate whether to apply smoothing in
    estimating mean methylation levels.}
  \item{smoothing.method}{Method for smoothing. Available options are
    "ma" for moving average, or "BSmooth" for BSmooth smoothing
    method. This will be ignored if smoothing=FALSE.}
  \item{smoothing.span}{The size of smoothing window, in
    basepairs. Default is 500.}
  \item{...}{Other parameters for BSmooth function in "bsseq" package.}
}

\details{
  This is the core function for DML/DMR detection. Tests are performed
  at each CpG site under the null hypothesis that two groups means are
  equal. There is an option for applying smoothing or not in estimating
  mean methylation levels. We recommend to use smoothing=TRUE for
  whole-genome BS-seq data, and smoothing=FALSE for sparser data such
  like from RRBS or hydroxyl-methylation data (TAB-seq). If there is not
  biological replicate, smoothing=TRUE is required. See "Single
  replicate" section for details.  
  
  The BS-seq count data are modeled as Beta-Binomial distribution, where
  the biological variations are captured by the dispersion
  parameter. The dispersion parameters are estimated through a shrinakge
  estimator based on a Bayesian hierarchical model. Then a Wald test is
  performed at each CpG site. 

  Due to the differences in coverages, some CpG sites are not covered in
  both groups, and the test cannot be performed. Those loci will be
  ignored in test and results will be "NA". 
}

\section{Single replicate}{
  When there is no biological replicate (one sample per treatment
  group), smoothing=TRUE is required. With the smoothed means,
  dispersion can still be estimated via the shrinkage estimator. This
  procedure uses data from neighboring CpG sites as "pseudo-replicate"
  for estimating biological variance. Results show that the estimation
  works reasonably well, and the DMR calling results are better than
  simply using the difference of two smoothed means.}

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

\section{Estimating mean methylation levels}{
  When smoothing=FALSE, the mean methylation levels are estimated based
  on the ratios of methylated and total read counts, and the spatial
  correlations among nearby CpG sites are ignored. When smoothing=TRUE,
  smoothing based on moving average or the BSmooth method is used to
  estimate the mean methylaion level at each site. Moving average is
  recommended because it is much faster than BSmooth, and the results
  are reasonable similar in terms of mean estimation, dispersion
  estimation, and DMR calling results.
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
makeBSseqData, callDML, callDMR
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

##  DML test without smoothing 
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
head(dmlTest)

## For whole-genome BS-seq data, perform DML test with smoothing
require(bsseqData)
data(BS.cancer.ex)
## take a small portion of data and test
BSobj <- BS.cancer.ex[10000:15000,]
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"), 
   smoothing=TRUE, smoothing.span=500)
head(dmlTest)
}
}
