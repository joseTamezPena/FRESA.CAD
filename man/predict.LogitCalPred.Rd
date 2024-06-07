\name{predict.LogitCalPred}
\alias{predict.LogitCalPred}
\title{Predicts calibrated probabilities}
\description{
	This function predicts the calibrated probability of a binary outcome
}
\usage{
	\method{predict}{LogitCalPred}(object,...)
}
\arguments{
	\item{object}{
		An object of class LogitCalPred 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the calibrated probability
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{calBinProb}}}
\keyword{Model_Prediction}
