\name{predict.FRESA_BESS}
\alias{predict.FRESA_BESS}
\title{Predicts \code{\link{BESS}} models}
\description{
	This function predicts the outcome from a BESS model
}
\usage{
	\method{predict}{FRESA_BESS}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_BESS 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the predict BESS object
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{BESS}}}
\keyword{Model_Prediction}
