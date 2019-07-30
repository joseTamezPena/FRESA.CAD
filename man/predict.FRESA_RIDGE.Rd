\name{predict.FRESA_RIDGE}
\alias{predict.FRESA_RIDGE}
\title{Predicts \code{\link{LM_RIDGE_MIN}} models}
\description{
	This function predicts the outcome from a LM_RIDGE_MIN model
}
\usage{
	\method{predict}{FRESA_RIDGE}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_RIDGE 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{LM_RIDGE_MIN}}}
\keyword{Model_Prediction}
