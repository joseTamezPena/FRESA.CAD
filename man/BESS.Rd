\name{BESS}
\alias{BESS}
\alias{BESS_GSECTION}
\title{CV BeSS fit}
\description{
	Fits a \code{BeSS::bess} object to the data, and return the selected features  
}
\usage{
	BESS(formula = formula, data=NULL, method="sequential", ic.type="BIC",...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the bess model
	}
	\item{method}{
	BeSS: Methods to be used to select the optimal model size
	}
	\item{ic.type}{
	BeSS: Types of best model returned.
	}
	\item{...}{
	Parameters to be passed to the \code{BeSS::bess} function
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
	The character vector of the model features according to BeSS type
	}
}
\author{Jorge Orozco}
\seealso{\code{BeSS::bess}}
\keyword{Model_Generation}
