\name{showOneDMR}
\alias{showOneDMR}
\title{
  Visualze the count data for one DMR
}
\description{
  Given one DMR and an BSseq object, this function generate a multiple
  panel figure, each for a sample, to visualze the counts. There is a
  bar at each CpG, the gray bar shows the total coverage, and the black
  bar shows the methylated count. 
}
\usage{
showOneDMR(OneDMR, BSobj, ext = 500, ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OneDMR}{A data frame with one row representing one DMR. It must
    have chr, start, and end columns. This is typically a row from the
    result generated from callDMR.}
  \item{BSobj}{An object of class BSseq.}
  \item{ext}{The amount (in bps) the plotting region should be extended
    in both directions.}
  \item{ylim}{Y-axis limit.}
}

\value{
  This function only generates a figure and has no return values.
}

\author{
Hao Wu <hao.wu@emory.edu>
}

\seealso{
  callDMR
}
\examples{
\dontrun{
require(bsseq)
require(bsseqData)
data(BS.cancer.ex)

## takea small portion of data and test
BSobj <- BS.cancer.ex[140000:150000,]
dmlTest <- DMLtest(BSobj, group1=c("C1", "C2", "C3"), group2=c("N1","N2","N3"),
   smoothing=TRUE, smoothing.span=500)

## call DMR based on test results
dmrs <- callDMR(dmlTest)

## visualize one DMR
showOneDMR(dmrs[1,], BSobj)
}

}
