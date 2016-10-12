\name{barcode}
\docType{data}
\alias{barcode}
\alias{barcode.cid}
\title{The barcode dataset}
\description{A flow cytometry data that is used barcode to measure many
  samples simultaneusly}
\usage{data(barcode)}
\format{
  An object (barcode) of data frame with 180912 rows and 3 columns and a
  vector (barcode.cid) for the cluster labels accoring to
  the manual gating. 
}
\source{The data is a random subset of
  the full data set for Figure 3A of the paper (Sugar et al 2010),
   This subset
  was used to do all comparisons in the paper (Ge et al 2012) with other clustering
  algorithms.}
\references{
  Sugar I. P. and Sealfon S. C., Misty Mountain clustering:
  application to fast unsupervised flow cytometry gating,
  \emph{BMC Bioinformatcs, 2010}, 11:502. \cr
  
  Ge Y. et al, flowPeaks: a fast unsupervised clustering for flow
  cytometry data via K-means and density peak finding, 2012,
  Bioinformatics  8(15):2052-8
  }
\keyword{datasets}
