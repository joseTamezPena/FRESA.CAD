\name{baggedModel}
\alias{baggedModel}
\title{Get the bagged model from a list of models}
\description{
	This function will take the frequency-ranked of variables and the list of models to create a single bagged model
}
\usage{
	baggedModel(modelFormulas,
	            data,
	            type=c("LM","LOGIT","COX"),
	            Outcome=NULL,
	            timeOutcome=NULL,
	            frequencyThreshold=0.025,
	            univariate=NULL,
				useFreq=TRUE,
				n_bootstrap=1,
				equifreqCorrection=0
	            )
}
\arguments{
	\item{modelFormulas}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the time to outcome
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{frequencyThreshold}{
		set the frequency the threshold of the frequency of features to be included in the model)
	}
	\item{univariate}{
		The FFRESA.CAD univariate analysis matrix
	}
	\item{useFreq}{
		Use the feature frequency to order the formula terms. If set to a positive value is the number of estimation loops
	}
	\item{n_bootstrap}{
		if greater than 1, defines the number of bootstraps samples to be used
	}
	\item{equifreqCorrection}{
		Indicates the average size of repeated features in an equivalent model
	}
}
\value{
  \item{bagged.model}{
		the bagged model
	}
	\item{formula}{
		the formula of the model
	}
	\item{frequencyTable}{
		the table of variables ranked by their model frequency
	}
	\item{faverageSize}{
		the average size of the models
	}
	\item{formulaNetwork}{
        The matrix of interaction between formulas
	}
	\item{Jaccard.SM}{
	   The Jaccard Stability Measure of the formulas
	}
	\item{coefEvolution}{
	   The evolution of the coefficients
	}
	\item{avgZvalues}{
	   The average Z value of each coefficient
	}	
	\item{featureLocation}{
	   The average location of the feature in the formulas
	}	
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{ensemblePredict}}}

\keyword{Model_Generation}