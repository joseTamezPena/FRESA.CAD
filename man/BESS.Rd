\name{BESS}
\alias{BESS}
\title{CV BeSS fit}
\description{
	Fits a \code{BeSS::bess} object to the data, and sets the prediction to use the features that created the minimum CV error or one SE. 
}
\usage{
	BESS(formula = formula,data=NULL,...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the KNN method
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
	The list of features to be used
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{glmnet::cv.glmnet}}
\keyword{Model_Generation}
