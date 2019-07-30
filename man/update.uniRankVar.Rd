\name{update.uniRankVar}
\alias{update}
\alias{update.uniRankVar}
\title{Update the univariate analysis using new data}
\description{
	This function updates the results from an univariate analysis using a new data set
}
\usage{
	\method{update}{uniRankVar}(object,
	           ...)
}
\arguments{
	\item{object}{
		A list with the results from the \code{uniRankVar} function
	}
	\item{...}{
		Additional parameters  to be passed to the \code{uniRankVar} function, used to update the univariate analysis
	}
}
\value{
	A list with the same format as the one yielded by the \code{uniRankVar} function
 }
\author{Jose G. Tamez-Pena}%%\note{
\seealso{\code{\link{uniRankVar}}}

\keyword{Data_Inspection}