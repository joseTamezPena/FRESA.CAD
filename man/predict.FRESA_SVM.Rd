\name{predict.FRESA_SVM}
\alias{predict.FRESA_SVM}
\title{Predicts \code{\link{TUNED_SVM}} models}
\description{
	This function predicts the outcome from a TUNED_SVM model
}
\usage{
	\method{predict}{FRESA_SVM}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESA_SVM 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	the predict e1071::svm object
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{TUNED_SVM}}}
\keyword{Model_Prediction}
