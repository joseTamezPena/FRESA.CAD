\name{updateModel.Res}
\alias{updateModel.Res}
\title{Update the NeRI-based model using new data or new threshold values}
\description{
	This function will take the frequency-ranked set of variables and will generate a new model with terms that meet the net residual improvement (NeRI) threshold criteria.
}
\usage{
	updateModel.Res(Outcome, 
	                covariates = "1", 
	                pvalue = c(0.025, 0.05),
	                VarFrequencyTable, 
	                variableList, 
	                data, 
	                type = c("LM", "LOGIT", "COX"),
	                testType=c("Binomial", "Wilcox", "tStudent"), 
	                lastTopVariable = 0, 
	                timeOutcome = "Time",
	                maxTrainModelSize = -1,
	                p.thresholds = NULL
	                )
}
\arguments{
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{covariates}{
		A string of the type "1 + var1 + var2" that defines which variables will always be included in the models (as covariates)
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model
	}
	\item{VarFrequencyTable}{
		An array with the ranked frequencies of the features, (e.g. the \code{ranked.var} value returned by the \code{ForwardSelection.Model.Res} function)
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{lastTopVariable}{
		The maximum number of variables to be tested
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
	\item{p.thresholds}{
		The p.value thresholds estimated in forward selection
	}
}
\value{
  \item{final.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the final model
	}
	\item{var.names}{
		A vector with the names of the features that were included in the final model
	}
	\item{formula}{
		An object of class \code{formula} with the formula used to fit the final model
	}
	\item{z.NeRI}{
		A vector in which each element represents the \emph{z}-score of the NeRI, associated to the \code{testType}, for each feature found in the final model
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{updateModel.Bin}}}

\keyword{Model_Generation}