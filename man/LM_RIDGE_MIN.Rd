\name{LM_RIDGE_MIN}
\alias{LM_RIDGE_MIN}
\title{Ridge Linear Models}
\description{
	FRESA wrapper to fit \code{MASS::lm.ridge} object to the data and returning the coef with minimum GCV
}
\usage{
	LM_RIDGE_MIN(formula = formula,data=NULL,...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the method
	}
	\item{...}{
	Parameters to be passed to the MASS::lm.ridge function
	}
}

\value{
	\item{fit}{
	The \code{MASS::lm.ridge} fitted object
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{MASS::lm.ridge}}
\keyword{Model_Generation}
