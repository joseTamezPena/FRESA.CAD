\name{predict.BAGGS}
\alias{predict.BAGGS}
\title{Predicts \code{\link{baggedModel}} bagged models}
\description{
	This function predicts the class of a BAGGS generated models
}
\usage{
	\method{predict}{BAGGS}(object,...)
}
\arguments{
	\item{object}{
		An object of class BAGGS 
	}
	\item{...}{
		A list with: testdata=testdata. 
	}
}
\value{
	a named list with the predicted class of every data sample
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{baggedModel}}}
\keyword{Bagged_Prediction}
