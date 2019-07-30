\name{reportEquivalentVariables}
\alias{reportEquivalentVariables}
\title{Report the set of variables that will perform an equivalent IDI discriminant function}
\description{
	Given a model, this function will report a data frame with all the variables that may be interchanged in the model without affecting its classification performance.
	For each variable in the model, this function will loop all candidate variables and report all of which result in an equivalent or better zIDI than the original model.
}
\usage{
	reportEquivalentVariables(object,
	                          pvalue = 0.05,
	                          data,
	                          variableList,
	                          Outcome = "Class",
	                          timeOutcome=NULL,
	                          type = c("LOGIT", "LM", "COX"),
	                          description = ".",
	                          method="BH",
	                          osize=0,
	                          fitFRESA=TRUE)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the IDI , allowed for a pair of variables to be considered equivalent
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{description}{
		The name of the column in \code{variableList} that stores the variable description
	}
	\item{method}{
		The method used by the p-value adjustment algorithm
	}
	\item{osize}{
		The number of features used for p-value adjustment
	}
	\item{fitFRESA}{
		if TRUE it will use the cpp based fitting method 
	}
}
\value{
	\item{pvalueList}{
		A list with all the unadjusted p-values of the equivalent features per model variable
	}
	\item{equivalentMatrix}{
	A data frame with three columns. The first column is the original variable of the model. 
	The second column lists all variables that, if interchanged, will not statistically affect the performance of the model.
	The third column lists the corresponding \emph{z}-scores of the IDI for each equivalent variable.
	}
	\item{formulaList}{
	 a character vector with all the equivalent formulas 
	}
	\item{equivalentModel}{
	 a bagged model that used all the equivalent formulas. The model size is limited by the number of observations
	}
	
}
\author{Jose G. Tamez-Pena}

\keyword{Model_Inspection}
