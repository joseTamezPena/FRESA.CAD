\name{predict.FRESA_NAIVEBAYES}
\alias{predict.FRESA_NAIVEBAYES}
\title{Predicts \code{\link{NAIVE_BAYES}} models}
\description{
	This function predicts the outcome from a FRESA_NAIVEBAYES model
}
\usage{
	\method{predict}{FRESA_NAIVEBAYES}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_NAIVEBAYES 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{NAIVE_BAYES}}}
\keyword{Model_Prediction}
