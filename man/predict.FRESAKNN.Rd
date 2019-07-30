\name{predict.FRESAKNN}
\alias{predict.FRESAKNN}
\title{Predicts \code{class::knn} models}
\description{
	This function predicts the outcome from a FRESAKNN model
}
\usage{
	\method{predict}{FRESAKNN}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESAKNN containing the KNN train set
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{KNN_method}}, \code{class::knn}}
\keyword{Model_Prediction}
