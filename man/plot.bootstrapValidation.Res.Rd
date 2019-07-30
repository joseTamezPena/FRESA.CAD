\name{plot.bootstrapValidation_Res}
\alias{plot.bootstrapValidation_Res}
\title{Plot ROC curves of bootstrap results}
\description{
	This function plots ROC curves and a Kaplan-Meier curve (when fitting a Cox proportional hazards regression model) of a bootstrapped model.
}
\usage{
	\method{plot}{bootstrapValidation_Res}(x,
	     xlab = "Years",
	     ylab = "Survival",
	     ...)
}
\arguments{
	\item{x}{
		A \code{bootstrapValidation_Res} object
	}
	\item{xlab}{
		The label of the \emph{x}-axis
	}
	\item{ylab}{
		The label of the \emph{y}-axis
	}
	\item{...}{
		Additional parameters for the plot
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{plot.bootstrapValidation_Bin}}}

\keyword{Model_Inspection}
