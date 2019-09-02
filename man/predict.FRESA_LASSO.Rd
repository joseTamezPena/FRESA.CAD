\name{predict.FRESA_GLMNET}
\alias{predict.FRESA_GLMNET}
\title{Predicts GLMNET fitted objects }
\description{
	This function predicts the outcome from a FRESA_GLMNET fitted object
}
\usage{
	\method{predict}{FRESA_GLMNET}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_GLMNET containing the model to be analyzed
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{GLMNET}}}
\keyword{Model_Prediction}
