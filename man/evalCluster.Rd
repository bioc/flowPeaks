\name{evalCluster}
\alias{evalCluster}
\title{evaulate the result of a clustering algorihm
 by comparing it with the gold
  standard}
\description{
  This function takes the cluster labels of the two clusterings, one is
  based on the gold standard, the other is a candidate clusterign, and
  compute one of the three metrics to assess the candidate clustering performance.
}
\usage{
evalCluster(gs,cand,method=c("Rand.index","Fmeasure","Vmeasure"),
                      rm.gs.outliers=TRUE)
}
\arguments{
  \item{gs}{A integer-valued vector of length n for the cluster labels of
    the gold standard clustering, where negative numbers
    such as -1 is for the outerliers}
\item{cand}{A integer-valued vector of length n for the cluster label of
a candidate clustering, where -1 is for the outliers}
\item{rm.gs.outliers}{Determining whether the outliers of the gold
standard clustering should be removed in the comparison}
\item{method}{A single character to indicate which one of three metrics
should be used to evaluate the clustering. The details are described in
Ge (2012) and references mentioned in that paper\cr
\describe{
\item{Rand.index}{The adjusted Rand.index}
\item{Fmeasure}{F-measure}
\item{Vmeasure}{V-measure}
}
}
}
\author{
  Yongchao Ge \email{yongchao.ge@gmail.com}
}
\references{
Ge Y. et al, flowPeaks: a fast unsupervised clustering for flow
  cytometry data via K-means and density peak finding, 2012, Bioinformatics, in press. 
  }

\keyword{cluster,utiliteis}
\seealso{\code{\link{flowPeaks}}}   
    
