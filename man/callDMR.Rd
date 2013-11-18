\name{callDMR}
\alias{callDMR}

\title{
  Function to detect differntially methylated regions (DMR) for two group
  comparisons of bisulfite sequencing (BS-seq) data.
}
\description{
  This function takes the results from DML detection procedure ('callDML'
  function) and calls DMRs. Regions will CpG sites being statistically
  significant are detected as DMRs. Nearby DMRs are merged into longer
  ones. Some restrictions including the minimum length, minimum number
  of CpG sites, etc. are applied. Note that this function doesn't
  consider the spatial correlation among nearby CpG sites, which will be
  in the future development plan. 
}
\usage{
callDMR(DMLresult, p.threshold, minlen = 50, minCG = 3, dis.merge = 50, pct.sig = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMLresult}{
    A data frame representing the results for DML detection. This should
    be the result returned from 'callDML' function.
  }
  \item{p.threshold}{
    A threshold of p-values for calling DMR. Loci with p-values less
    than this threshold will be picked and joint to form the DMRs.
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

\value{A data frame for DMRs. Each row is for a DMR. 
  Rows are sorted by chromosome number and genomic coordinates. The columns
  include:
  \item{chr}{Chromosome number.}
  \item{start, end}{Genomic coordinates.}
  \item{length}{Length of the DMR, in bps.}
  \item{nCG}{Number of CpG sites contained in the DMR.}
  \item{meanMethy1, meanMethy2}{Average methylation levels in two
    conditions.}
  \item{diff.Methy}{The difference in the methylation levels between two
  conditions.} 
}

\details{
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
callDML
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
dmls <- callDML(BS1, BS2)
dmrs <- callDMR(dmls, p.threshold=0.001)
}
}
