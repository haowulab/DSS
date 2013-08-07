\name{SeqCountSet-class}
\Rdversion{1.1}
\docType{class}
\alias{SeqCountSet-class}
\alias{SeqCountSet}

% Constructor-like function:
\alias{newSeqCountSet}

\title{Class \code{"SeqCountSet"} - container for count data from
  sequencing experiment}

\description{
  This class is the main container for storing data from sequencing
  technology. It is directly inherited fro 'ExpressionSet' class, with two
  more fields 'normalizationFactor' for normalization factors and
  'dispersion' for gene-wise dispersions. 
}

\section{Slots}{
    \describe{
      \item{\code{normalizationFactor}:}{Normalization factor for counts.}
      \item{\code{dispersion}:}{Gene-wise dispersions.}
      \item{\code{experimentData}:}{See 'ExpressionSet'.}
      \item{\code{assayData}:}{See 'ExpressionSet'.}
      \item{\code{phenoData}:}{See 'ExpressionSet'.}
      \item{\code{featureData}:}{See 'ExpressionSet'.}
      \item{\code{annotation}:}{See 'ExpressionSet'.}
      \item{\code{protocolData}:}{See 'ExpressionSet'.}
  }
}
  
\section{Extends}{
  Class \code{"\linkS4class{ExpressionSet}"}, directly.
  Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
  Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{newSeqCountSet(counts,
	designs,normalizationFactor,featureData)}:
      Creates a 'SeqCountSet' object.
      \describe{
	\item{\code{counts}}{A matrix of integers with rows corresponding to
	  genes and columns for samples.}
	\item{\code{designs}}{A vector or data frame representing experimental
	  design. The length of the vector or number of rows of the
	  data frame must match the number of columns of
	  input counts. This field can be accessed using 'pData' function.}
	\item{\code{normalizationFactor}}{A vector or matrix of
	  normalization	factors for the counts. }
	\item{\code{featureData}}{Additional information for genes
	  as an 'AnnotatedDataFrame' object. This field can be
	  access by using 'featureData' function.}
      }
    }
  }
}


\section{Methods}{
  \describe{
    \item{dispersion, dispersion<- :}{Access and set gene-wise dispersions.}
    \item{normalizationFactor, normalizationFactor<- :}{Access and set
      normalization factors.}
  }
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\seealso{
  dispersion, normalizationFactor
}

\note{
  This is similar to 'CountDataSet' in DESeq or 'DGEList' in edgeR. 
}

\examples{
counts=matrix(rpois(600, 10), ncol=6)
designs=c(0,0,0,1,1,1)
seqData=newSeqCountSet(counts, designs)
seqData
pData(seqData)
head(exprs(seqData))

## multiple factor designs
design=data.frame(gender=c(rep("M",2*nreps), rep("F",2*nreps)), strain=rep(c("WT", "Mutant"),2*nreps))
X=model.matrix(~gender+strain, data=design)
counts=matrix(rpois(800, 10), ncol=8)
seqData=newSeqCountSet(counts, as.data.frame(X))
seqData
pData(seqData)

}
\keyword{classes}
