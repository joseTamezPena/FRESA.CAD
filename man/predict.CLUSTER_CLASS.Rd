\name{predict.CLUSTER_CLASS}
\alias{predict.CLUSTER_CLASS}
\title{Predicts \code{\link{ClustClass}} outcome}
\description{
	This function predicts the outcome from a ClustClass classifier
}
\usage{
	\method{predict}{CLUSTER_CLASS}(object,...)
}
\arguments{
	\item{object}{
		An object of class CLUSTER_CLASS 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the predict of a hierarchical ClustClass classifier
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{ClustClass}}}
\keyword{Model_Prediction}
