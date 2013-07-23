\name{makeBSseqData}
\alias{makeBSseqData}

\title{
  Create an object of BSseq class from several data frames.
}
\description{
  This is an utility function to merge BS-seq data from replicated
  experiment and create an object of BSseq class.
  
  After sequence alignment and proper processing, the BS-seq data can be
  summarized by following information at each C position (mostly CpG
  sites, with a little CH): chromosome number, genomic coordinate, total
  number of reads covering the position, and number of reads showing
  methylation at this position. For replicated samples, the data need to
  be merged based on the chromosome number and genomic coordinates. This
  function provide such functionality. It takes replicated data as a
  list of data frames, merged them, and create a BSseq object.
}
\usage{
makeBSseqData(dat, sampleNames)
}

\arguments{
  \item{dat}{A list of multiple data frames from biological
    replicates. Each element represents data from one replicate. The data
    frame should contain following files: (1) Chromosome number; (2)
    Genomic coordinates; (3) Read coverage of the position from BS-seq
    data; (4) Number of reads showing methylation of the position.}
  \item{sampleNames}{A vector of characters for the sample names. The
    length of the vector should match the length of the input list.}
}

\value{
  An object of 'BSseq' class.
}

\author{
Hao Wu <hao.wu@emory.edu>
}

\seealso{
callDML
}

\examples{
require(bsseq)
## first read in methylation data.
path <- file.path(system.file(package="DSS"), "extdata")
dat1.1 <- read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 <- read.table(file.path(path, "cond1_2.txt"), header=TRUE)

## make BSseq objects
BS1 <- makeBSseqData( list(dat1.1, dat1.2), paste("cond1",1:2,sep=".") )
BS1
}
