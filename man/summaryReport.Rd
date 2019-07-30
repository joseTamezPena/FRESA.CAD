\name{summaryReport}
\alias{summaryReport}
\title{Report the univariate analysis, the cross-validation analysis and the correlation analysis}
\description{
	This function takes the variables of the cross-validation analysis and extracts the results from the univariate and correlation analyses.
	Then, it prints the cross-validation results, the univariate analysis results, and the correlated variables.
	As output, it returns a list of each one of these results.
}
\usage{
	summaryReport(univariateObject,
	              summaryBootstrap,
	              listOfCorrelatedVariables = NULL,
	              digits = 2)
}
\arguments{
	\item{univariateObject}{
		A data frame that contains the results of the \code{univariateRankVariables} function
	}
	\item{summaryBootstrap}{
		A list that contains the results of the \code{summary.bootstrapValidation_Bin} function
	}
	\item{listOfCorrelatedVariables}{
		A matrix that contains the \code{correlated.variables} value from the results obtained with the \code{listTopCorrelatedVariables} function
	}
	\item{digits}{
		The number of significant digits to be used in the print function
	}
}
\value{
	\item{performance.table}{
		A matrix with the tabulated results of the blind test accuracy, sensitivity, specificities, and area under the ROC curve
	}
	\item{coefStats}{
		A data frame that lists all the model features along with its univariate statistics and bootstrapped coefficients
	}
	\item{cor.varibles}{
		A matrix that lists all the features that are correlated to the model variables
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{summary.bootstrapValidation_Bin}}}

\keyword{Model_Inspection}
