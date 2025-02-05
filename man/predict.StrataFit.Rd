\name{predict.StrataFit}
\alias{predict.StrataFit}
\title{Predicts \code{\link{StrataFit}} targets}
\description{
	This function predicts the outcome of the StrataFit 
}
\usage{
	\method{predict}{StrataFit}(object,...)
}
\arguments{
	\item{object}{
		An object of class StrataFit 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	numeric vector with the target predictions
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{StrataFit}}}
\keyword{Model_fitting}
