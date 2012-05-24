\name{flowPeaks}
\alias{flowPeaks}
\title{Doing the flowPeaks analysis}
\description{This is the core function in the flowPeaks package. It
  generates the output of the cluster and information associated with
  each cluster, which can be used by the function \code{plot} for
  visualization}
\usage{
flowPeaks(x,tol=0.1,h0=1,h=1.5)
}
\arguments{
  \item{x}{
    a data matrix for the flow cytometry data, it needs to have
      at least two rows, and the names for each column should be unique.
      For a flowFrame data, its exprssion matrix slot
      should be used as x, where only channles of interest are selected
      (see the example below).
    }
  \item{tol}{The tolerance (between 0 and 1) when neighboring clusters
    should be considered to be merged}
  \item{h0}{The multiplier of the vaiarance matrix S0}
  \item{h}{The multiplier of the variance matrix S}
}
\value{
  It returns an object of class flowPeaks, which is a list of the
  following variables:
  \item{peaks.cluster}{An integer shows the cluster labels (between 1 and K
    for K clusters) for each cell. The clustering is based on the
    flowPeaks algorithm}
  \item{peaks}{A summary of the cluster information. It is a list
    with the following three variables:\itemize{
    \item cid: cluster labels, should always be 1:K;
    \item w: the weights of the K clusters;
    \item mu: The mean of all cells in the K clusters;
    \item S: The variance matrix of the K clusters. Note that each variance
    matrix for each cluster has been stacked as a column vector}
  }
  \item{kmeans.cluster}{An integer shows the cluster labels for the initial
    kmeans clustering}
  \item{kmeans}{A summary of the initial kmeans clustering. The meaning
    of the variables can be seens in the description of peaks above
  }
  \item{info}{The information that can be used for plot, and how the
    initial kmeans clustering and the final flowPeaks clustering are
    connected
  }
  \item{x}{The input data x}
  }
    
\seealso{\code{\link{plot.flowPeaks}}}
\examples{
##demonstrate how to use a flowFrame
\dontrun{
require(flowCore)
samp <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))
##do the clustering based on the asinh transforamtion of
##the first two FL channels
fp<-flowPeaks(asinh(samp@exprs[,3:4]))
plot(fp)
}

data(barcode)
fp<-flowPeaks(barcode[,c(1,3)])
plot(fp)

##to compare it with the gold standard
evalCluster(barcode.cid,fp$peaks.cluster,method="Vmeasure")

#to remove the outliers
fpc<-assign.flowPeaks(fp,fp$x)
plot(fp,classlab=fpc,drawboundary=FALSE,
  drawvor=FALSE,drawkmeans=FALSE,drawlab=TRUE)


#to adjust the cluster by increasing the tol,h0, h, which results
#in a smaller number of clusters
fp2<-adjust.flowPeaks(fp,tol=0.5,h0=2,h=2) 
summary(fp2)
print(fp) #an alternative of using summary(fp) 
}
\author{
  Yongchao Ge \email{yongchao.ge@gmail.com}
}
\references{
  Ge Y. et al, flowPeaks: a fast unsupervised clustering for flow
  cytometry data via K-means and density peak finding, 2012, Bioinformatics, in press.
  }
\keyword{cluster,multivariate,smooth}
  
