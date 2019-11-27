\name{predict.FRESA_HLCM}
\alias{predict.FRESA_HLCM}
\title{Predicts BOOST_BSWiMS models}
\description{
	This function predicts the outcome from a BOOST_BSWiMS model
}
\usage{
	\method{predict}{FRESA_HLCM}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_HLCM 
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
