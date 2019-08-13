\name{predict.GMVE_BSWiMS}
\alias{predict.GMVE_BSWiMS}
\title{Predicts \code{\link{GMVEBSWiMS}} outcome}
\description{
	This function predicts the outcome from a GMVEBSWiMS classifier
}
\usage{
	\method{predict}{GMVE_BSWiMS}(object,...)
}
\arguments{
	\item{object}{
		An object of class GMVE_BSWiMS 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the predict of a hierarchical GMVE-BSWiMS classifier
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{GMVEBSWiMS}}}
\keyword{Model_Prediction}
