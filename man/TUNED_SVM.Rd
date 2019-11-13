\name{TUNED_SVM}
\alias{TUNED_SVM}
\title{Naive Bayes Modeling}
\description{
	FRESA wrapper to fit grid-tuned \code{ e1071::svm} object 
}
\usage{
	TUNED_SVM(formula = formula,
	          data=NULL,
	          gamma = 10^(-5:-1),
	          cost = 10^(-3:1),
	          ...
	          )
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the method
	}
	\item{gamma}{
	The vector of possible gamma values
	}
	\item{cost}{
	The vector of possible cost values 
	}
	\item{...}{
	Parameters to be passed to the e1071::svm function
	}
}

\value{
	\item{fit}{
	The \code{e1071::svm} fitted object
	}
	\item{tuneSVM}{
	The \code{e1071::tune.svm} object
	}
}

\author{Jose G. Tamez-Pena}
\seealso{\code{e1071::svm}}
\keyword{Model_Generation}
