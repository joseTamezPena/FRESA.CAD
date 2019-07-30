\name{modelFitting}
\alias{modelFitting}
\title{Fit a model to the data}
\description{
	This function fits a linear, logistic, or Cox proportional hazards regression model to given data
}
\usage{
	modelFitting(model.formula,
	             data,
	             type = c("LOGIT", "LM", "COX","SVM"),
	             fitFRESA=TRUE,
	              ...)
}
\arguments{
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), Cox proportional hazards ("COX") or "SVM"
	}
	\item{fitFRESA}{
		if true it will perform use the FRESA cpp code for fitting
	}
	\item{...}{
		Additional parameters for fitting a default \code{glm} object
	}
}
\value{
	A fitted model of the type defined in \code{type}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Model_Inspection}
