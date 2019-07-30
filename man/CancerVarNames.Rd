\name{cancerVarNames}
\alias{cancerVarNames}
\docType{data}
\title{Data frame used in several examples of this package}
\description{
	This data frame contains two columns, one with names of variables, and the other with descriptions of such variables.
	It is used in several examples of this package.
	Specifically, it is used in examples working with the stage C prostate cancer data from the \code{rpart} package
}
\usage{
	data(cancerVarNames)
}
\format{
	A data frame with names and descriptions of the variables used in several examples
	\describe{
		\item{\code{Var}}{
			A column with the names of the variables
		}
		\item{\code{Description}}{
			A column with a short description of the variables
		}
  }
}
\examples{
	data(cancerVarNames)
}
\keyword{Datasets}