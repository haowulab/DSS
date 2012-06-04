\name{DSS-package}
\docType{package}
\alias{DSS}
\alias{DSS-package}
\title{
  Dispersion shrinkage for sequencing data
}
\description{
  DSS is an R library performing the differential expression analysis for
  RNA-seq count data. Compared with other similar packages (DESeq,
  edgeR), DSS implements a new dispersion shrinkage method to estimate the
  gene-specific biological variance. Extensive simulation results
  showed that DSS performs favorabily compared to DESeq and edgeR when
  the variation of biological variances is large.

  DSS only works for two group comparison at this time. We plan to
  extend the functionalities and make it work for more general
  experimental designs in the near future. 
}

\author{
  Hao Wu <hao.wu@emory.edu>
}

\keyword{package}
