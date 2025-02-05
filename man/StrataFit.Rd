\name{StrataFit}
\alias{StrataFit}
\alias{StrataFit_EM}
\title{Stratified fit of datasets}
\description{
	strata-wise fitting of models. The function will model the strata to predict non-labeled data sets}
\usage{
StrataFit(formula = formula, 
			data=NULL, 
			strataFitMethod=KNN_method, 
			targetFitMethod=LASSO_1SE,
			strataFit.control=NULL,
			...
		  )
}

\arguments{
	\item{formula}{
	the base formula to extract the outcome
	}
	\item{data}{
	the data to be used for training the method
	}
	\item{strataFitMethod}{
	the classification modeling function for the strata
	}
	\item{targetFitMethod}{
	The classification modeling for the target
	}
	\item{strataFit.control}{
	the parameters to the strata fitting.
	}
	\item{...}{
	parameters for the classification function
	}
}

\value{
	\item{fit_strata}{
	The strata modeling 
	}
	\item{fit_target}{
	A list of fitting models per strata
	}
	\item{theclasses}{
	The labels of the strata
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Model_Generation}
