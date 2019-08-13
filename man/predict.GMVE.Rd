\name{predict.GMVE}
\alias{predict.GMVE}
\title{Predicts \code{\link{GMVECluster}} clusters}
\description{
	This function predicts the class of a GMVE generated cluster
}
\usage{
	\method{predict}{GMVE}(object,...)
}
\arguments{
	\item{object}{
		An object of class GMVE 
	}
	\item{...}{
		A list with: testdata=testdata. thr=p.value threshold
	}
}
\value{
	a named list with the predicted class of every data sample
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{GMVECluster}}}
\keyword{Cluster_Prediction}
