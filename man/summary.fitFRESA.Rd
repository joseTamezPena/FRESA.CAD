\name{summary.fitFRESA}
\alias{summary}
\alias{summary.fitFRESA}
\title{Returns the summary of the fit}
\description{
	Returns a summary of fitted model created by the modelFitting function with the fitFRESA parameter set to TRUE
}
\usage{
	\method{summary}{fitFRESA}(object,
	type=c("Improvement","Residual"),
	ci=c(0.025,0.975),
	data=NULL,
	...)
}
\arguments{
	\item{object}{
		fitted model with the \code{modelFitting} function
	}
	\item{type}{
		the type of coefficient estimation
	}
	\item{ci}{
		lower and upper limit of the ci estimation
	}
	\item{data}{
		the data to be used for 95%CI estimation
	}
	\item{...}{
		parameters of the boostrap method
	}
}
\value{
	a list with the analysis results.
 }
\author{Jose G. Tamez-Pena}%%\note{
\seealso{\code{\link{modelFitting}},\code{\link{bootstrapValidation_Bin}},\code{\link{bootstrapValidation_Res}}}

\keyword{Model_Inspection}