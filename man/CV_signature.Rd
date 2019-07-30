\name{CVsignature}
\alias{CVsignature}
\title{Cross-validated Signature}
\description{
	A  formula based wrapper of the \code{\link{getSignature}} function
}
\usage{
	CVsignature(formula = formula,data=NULL,...)
}
\arguments{
	\item{formula}{
	The base formula
	}
	\item{data}{
	The data to be used for training the signature method
	}
	\item{...}{
	Parameters for the \code{\link{getSignature}} function 
	}
}

\value{
	\item{fit}{
	A \code{\link{getSignature}} object.
	}
	\item{method}{
	The distance method 
	}
	\item{variable.importance}{
	The named vector of relevant features
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{getSignature}},\code{\link{signatureDistance}} }

\keyword{Model_Generation}
