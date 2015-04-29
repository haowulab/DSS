\name{callDMR}
\alias{callDMR}

\title{
  Function to detect differntially methylated regions (DMR) for two group
  comparisons of bisulfite sequencing (BS-seq) data.
}
\description{
  This function takes the results from DML testing procedure ('callDML'
  function) and calls DMRs. Regions will CpG sites being statistically
  significant are detected as DMRs. Nearby DMRs are merged into longer
  ones. Some restrictions including the minimum length, minimum number
  of CpG sites, etc. are applied. 
}
\usage{
callDMR(DMLresult, delta=0.1, p.threshold=1e-5,
        minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMLresult}{
    A data frame representing the results for DML detection. This should
    be the result returned from 'DMLtest' function.
  }
  \item{delta}{A threshold for defining DMR. In DML detection procedure, a hypothesis
    test that the two groups means are equal is conducted at each CpG
    site. Here if 'delta' is specified, the function will compute the
    posterior probability that the difference of the means are greater
    than delta, and then construct DMR based on that.}
  \item{p.threshold}{
    A threshold of p-values for calling DMR. Loci with p-values less
    than this threshold will be picked and joint to form the DMRs. See
    'details' for more information.
  }
  \item{minlen}{
    Minimum length (in basepairs) required for DMR. Default is 50 bps. 
  }
  \item{minCG}{
    Minimum number of CpG sites required for DMR. Default is 3. 
  }
  \item{dis.merge}{
    When two DMRs are very close to each other and the distance (in bps)
    is less than this number, they will be merged into one. Default is
    50 bps. 
  }
  \item{pct.sig}{
    In all DMRs, the percentage of CG sites with significant p-values
    (less than p.threshold) must be greater than this threshold. Default
    is 0.5. This is mainly used for correcting the effects of merging of nearby
    DMRs. 
  }
}

\value{A data frame for DMRs. Each row is for a DMR. Rows are sorted by
  "areaStat", which is the sum of test statistics of all CpG sites in
  the region. The columns are:
  \item{chr}{Chromosome number.}
  \item{start, end}{Genomic coordinates.}
  \item{length}{Length of the DMR, in bps.}
  \item{nCG}{Number of CpG sites contained in the DMR.}
  \item{meanMethy1, meanMethy2}{Average methylation levels in two
    conditions.}
  \item{diff.Methy}{The difference in the methylation levels between two
    conditions.}
  \item{areaStat}{The sum of the test statistics of all CpG sites within
    the DMR.}
}

\details{
  The choices of 'delta' and 'p.threshold' are somewhat arbitrary.
  The default value for p-value threshold for calling DMR is 1e-5. The
  statistical test on loci level is less powerful when smoothing is NOT
  applied, so users can consider to use a less stringent criteria, such
  as 0.001, in order to get satisfactory number of DMRs.
  This function is reasonably fast since the computationally
  intesnsive part is in 'DMLtest'. Users can try different 
  p.threshold values to obtain satisfactory results.
  
  When specifying a 'delta' value, the posterior probability (pp) of
  each CpG site being DML is computed. Then the p.threshold is applied on
  1-pp, e.g., sites with 1-pp<p.threshold is deemed significant. In this
  case, the criteria for DMR calling is more stringent and users might
  consider to use a more liberal p.threshold in order to get more
  regions.
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
DMLtest, callDML
}

\examples{
\dontrun{
require(bsseq)
require(bsseqData)
data(BS.cancer.ex)

## takea smallportionof data and test
BSobj <- BS.cancer.ex[140000:150000,]
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"),
   smoothing=TRUE, smoothing.span=500)

## call DMR based on test results
dmrs <- callDMR(dmlTest)
head(dmrs)

## or one can specify a threshold for difference in methylation level
dmrs2 <- callDMR(dmlTest, delta=0.2)
head(dmrs2)

## visualize one DMR
showOneDMR(dmrs[1,], BSobj)
}
}
