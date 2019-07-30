\name{residualForFRESA}
\alias{residualForFRESA}
\title{Return residuals from prediction}
\description{
	Given a model and a new data set, this function will return the residuals of the predicted values.
	When dealing with a Cox proportional hazards regression model, the function will return the Martingale residuals.
}
\usage{
	residualForFRESA(object,
	                 testData,
	                 Outcome,
	                 eta = 0.05)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{testData}{
		A data frame  where all variables are stored in different columns, with the data set to be predicted
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{eta}{
		The weight of the contribution of the Martingale residuals, or 1 - the weight of the contribution of the classification residuals (only needed if \code{object} is of class \code{coxph})
	}
}
\value{
	A vector with the residuals (i.e. the differences between the predicted and the real outcome)
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Model_Inspection}
