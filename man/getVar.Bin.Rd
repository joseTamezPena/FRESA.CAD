\name{getVar.Bin}
\alias{getVar.Bin}
\title{Analysis of the effect of each term of a binary classification model by analysing its reclassification performance}
\description{
	This function provides an analysis of the effect of each model term by comparing the binary classification performance between the Full model and the model without each term.
	The model is fitted using the train data set, but probabilities are predicted for the train and test data sets.
	Reclassification improvement is evaluated using the \code{improveProb} function (\code{Hmisc} package).
	Additionally, the integrated discrimination improvement (IDI) and the net reclassification improvement (NRI) of each model term are reported.
}
\usage{
	getVar.Bin(object,
	                       data,
	                       Outcome = "Class", 
	                       type = c("LOGIT", "LM", "COX"),
	                       testData = NULL,
	                       callCpp=TRUE)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analysed
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testData}{
		A data frame similar to \code{data}, but with a data set to be independently tested. If \code{NULL}, \code{data} will be used.
	}
	\item{callCpp}{
		is set to true it will use the c++ implementation of improvement.
	}
}
\value{
	\item{z.IDIs}{
		A vector in which each term represents the \emph{z}-score of the IDI obtained with the Full model and the model without one term
	}
	\item{z.NRIs}{
		A vector in which each term represents the \emph{z}-score of the NRI obtained with the Full model and the model without one term
	}
	\item{IDIs}{
		A vector in which each term represents the IDI obtained with the Full model and the model without one term
	}
	\item{NRIs}{
		A vector in which each term represents the NRI obtained with the Full model and the model without one term
	}
	\item{testData.z.IDIs}{
		A vector similar to \code{z.IDIs}, where values were estimated in \code{testdata}
	}
	\item{testData.z.NRIs}{
		A vector similar to \code{z.NRIs}, where values were estimated in \code{testdata}
	}
	\item{testData.IDIs}{
		A vector similar to \code{IDIs}, where values were estimated in \code{testdata}
	}
	\item{testData.NRIs}{
		A vector similar to \code{NRIs}, where values were estimated in \code{testdata}
	}
	\item{uniTrainAccuracy}{
		A vector with the univariate train accuracy of each model variable
	}
	\item{uniTestAccuracy}{
		A vector with the univariate test accuracy of each model variable
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{getVar.Res}}}

\keyword{Model_Inspection}