\name{predict.FRESA_LASSO}
\alias{predict.FRESA_LASSO}
\title{Predicts LASSO fitted objects }
\description{
	This function predicts the outcome from a FRESA_LASSO fitted object
}
\usage{
	\method{predict}{FRESA_LASSO}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_LASSO containing the model to be analyzed
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{LASSO}}}
\keyword{Model_Prediction}
