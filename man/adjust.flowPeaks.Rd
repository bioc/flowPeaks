\name{adjust.flowPeaks}
\alias{adjust.flowPeaks}
\title{Adjusting the smoothing and merging behavior of the flowPeaks
  results}
\description{Adjusting the smoothing and merging behavior of the flowPeaks
  results by changing the multiplers of the covariance matrix and the
  tolerance level for joining two peaks}
\usage{
adjust.flowPeaks(object,tol,h0,h,\dots)
}
\arguments{
  \item{object}{The output from the function \code{\link{flowPeaks}}}
  \item{tol}{See \code{\link{flowPeaks}}}
  \item{h0}{See \code{\link{flowPeaks}}}
  \item{h}{See \code{\link{flowPeaks}}}
  \item{\dots}{Optional additional arguments.  At present no additional
    arguments are used.}
}
\value{
  It returns an updated object of class flowPeaks, the detail defintion
  of which
  can be seen in
  \code{\link{flowPeaks}}.
  }
\author{
  Yongchao Ge \email{yongchao.ge@gmail.com}
}
\seealso{\code{\link{flowPeaks}}} 
\keyword{adjust}
