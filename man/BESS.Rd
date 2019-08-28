\name{BESS}
\alias{BESS}
\title{CV BeSS fit}
\description{
	Fits a \code{BeSS::bess} object to the data, and return the selected features  
}
\usage{
	BESS(formula = formula,data=NULL,...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the bess model
	}
	\item{...}{
	Parameters to be passed to the bess function
	}
}

\value{
	\item{fit}{
	The \code{BsSS::bess} fitted object
	}
	\item{formula}{
	The formula
	}
	\item{usedFeatures}{
	The list of features used by fit
	}
	\item{selectedfeatures}{
	The character vector of the model features acording to BeSS type
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{BeSS::bess}}
\keyword{Model_Generation}
