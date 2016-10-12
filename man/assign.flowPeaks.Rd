\name{assign.flowPeaks}
\alias{assign.flowPeaks}
\title{Obtain the flowPeaks cluster lables with the option of identifying
  outliers and applying to a new data set}
\description{The function takes a flowPeaks output and a new data set
  (or could be the same dataset that generated the flowPeaks), and
  compute the cluster label assignment}
\usage{
assign.flowPeaks(fp,A,tol=0.01,fc=0.8)
}
\arguments{
  \item{fp}{an object of class flowPeaks, the output from the function
    \code{\link{flowPeaks}} or  \code{\link{adjust.flowPeaks}}}
  \item{A}{A data matrix with the same number of columns as the data
    that geneterated fp}
  \item{tol}{All points where the probability
    density is less than tol (default is
    1\%) of the peak
    denisty of that cluster are labled as outliers. If tol is set
    0, no outliers according to this rule. The details can be seen in
    the first equation of Section 2.5 in
    the flowPeaks manuscript (Ge et al 2012)}
  \item{fc}{All points  where the classified cluster
    contributes less than fc (default is 80\%)
    of overall denstiy  are labeled as
    outliers. if fc is set to 0\%, no outliers can be found
    according to this rule. The details can be seen in the second
    equation of
    Section 2.5 in
    the flowPeaks manuscript (Ge et al 2012)}
}
\value{
  It returns the class label assignment of each data point, where -1
  indicates outliers. When A is the same data that generated fp,
  If tol is 1 and fc is 0, the returned
  labels are the same as fp$peaks.cluster.
}
\author{
  Yongchao Ge \email{yongchao.ge@gmail.com}
}
\references{
  Ge Y. et al, flowPeaks: a fast unsupervised clustering for flow
  cytometry data via K-means and density peak finding, 2012,
  Bioinformatics,  8(15):2052-8
  }
\keyword{cluster,multivarite}
\seealso{\code{\link{flowPeaks}}}   
