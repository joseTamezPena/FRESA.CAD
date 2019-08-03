\name{predict.FRESA_BOOST}
\alias{predict.FRESA_BOOST}
\title{Predicts \code{\link{FRESA_BOOST}} models}
\description{
	This function predicts the outcome from a BOOST_BSWiMS model
}
\usage{
	\method{predict}{FRESA_BOOST}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_BOOST 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the predict of boosted BSWiMS
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{BSWiMS.model}}}
\keyword{Model_Prediction}
