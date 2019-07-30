\name{plot.bootstrapValidation_Bin}
\alias{plot}
\alias{plot.bootstrapValidation_Bin}
\title{Plot ROC curves of bootstrap results}
\description{
	This function plots ROC curves and a Kaplan-Meier curve (when fitting a Cox proportional hazards regression model) of a bootstrapped model.
}
\usage{
	\method{plot}{bootstrapValidation_Bin}(x,
	     xlab = "Years",
	     ylab = "Survival",
		 strata.levels=c(0),
	     main = "ROC",
		 cex=1.0,
		 ...)
}
\arguments{
	\item{x}{
		A \code{bootstrapValidation_Bin} object
	}
	\item{xlab}{
		The label of the \emph{x}-axis
	}
	\item{ylab}{
		The label of the \emph{y}-axis
	}
	\item{strata.levels}{
		stratification level for the Kaplan-Meier plots
	}
	\item{main}{
		Main Plot title
	}
	\item{cex}{
		The text cex
	}
	\item{...}{
		Additional parameters for the generic \code{plot} function
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{plot.bootstrapValidation_Res}}}

\keyword{Model_Inspection}
