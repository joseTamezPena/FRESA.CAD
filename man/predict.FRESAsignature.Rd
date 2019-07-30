\name{predict.FRESAsignature}
\alias{predict.FRESAsignature}
\title{Predicts \code{\link{CVsignature}} models}
\description{
	This function predicts the outcome from a FRESAsignature model
}
\usage{
	\method{predict}{FRESAsignature}(object,...)
}
\arguments{
	\item{object}{
		An object of class FRESAsignature 
	}
	\item{...}{
		A list with: testdata=testdata
	}
}
\value{
	A vector of the predicted values
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{CVsignature}},\code{\link{getSignature}},\code{\link{signatureDistance}}}
\keyword{Model_Prediction}
