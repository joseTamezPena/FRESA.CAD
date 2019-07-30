\name{summary.bootstrapValidation_Bin}
\alias{summary.bootstrapValidation_Bin}
\title{Generate a report of the results obtained using the bootstrapValidation_Bin function}
\description{
	This function prints two tables describing the results of the bootstrap-based validation of binary classification models.
	The first table reports the accuracy, sensitivity, specificity and area under the ROC curve (AUC) of the train and test data set, along with their confidence intervals.
	The second table reports the model coefficients and their corresponding integrated discrimination improvement (IDI) and net reclassification improvement (NRI) values.
}
\usage{
	\method{summary}{bootstrapValidation_Bin}(object,
	        ...)
}
\arguments{
	\item{object}{
		An object of class \code{bootstrapValidation_Bin}
	}
	\item{...}{
		Additional parameters for the generic \code{summary} function
	}
}
\value{
	\item{performance}{
		A vector describing the results of the bootstrapping procedure
	}
	\item{summary}{
		An object of class \code{summary.lm}, \code{summary.glm}, or \code{summary.coxph} containing a summary of the analyzed model
	}
	\item{coef}{
		A matrix with the coefficients, IDI, NRI, and the 95\% confidence intervals obtained via bootstrapping
	}
	\item{performance.table}{
		A matrix with the tabulated results of the blind test accuracy, sensitivity, specificities, and area under the ROC curve
	}  
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{summaryReport}}}

\keyword{Model_Inspection}
