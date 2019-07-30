\name{predict.fitFRESA}
\alias{predict}
\alias{predict.fitFRESA}
\title{Linear or probabilistic prediction}
\description{
	This function returns the predicted outcome of a specific model.
	The model is used to generate linear predictions.
	The probabilistic values are generated using the logistic transformation on the linear predictors.
}
\usage{
	\method{predict}{fitFRESA}(object,
	                ...)
}
\arguments{
	\item{object}{
		An object of class fitFRESA containing the model to be analyzed
	}
	\item{...}{
		A list with: testdata=testdata;predictType=c("linear","prob") and impute=FALSE.
		If impute is set to TRUE it will use the object model to impute missing data
	}
}
\value{
	A vector with the predicted values
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{nearestNeighborImpute}}}
\keyword{Model_Prediction}
