\name{NAIVE_BAYES}
\alias{NAIVE_BAYES}
\title{Naive Bayes Modeling}
\description{
	FRESA wrapper to fit \code{naivebayes::naive_bayes} object to the data 
}
\usage{
	NAIVE_BAYES(formula = formula,data=NULL,...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the method
	}
	\item{...}{
	Parameters to be passed to the naivebayes::naive_bayes function
	}
}

\value{
	\item{fit}{
	The \code{naivebayes::naive_bayes} fitted object
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{naivebayes::naive_bayes}}
\keyword{Model_Generation}
