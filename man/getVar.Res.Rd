\name{getVar.Res}
\alias{getVar.Res}
\title{Analysis of the effect of each term of a linear regression model by analysing its residuals}
\description{
	This function provides an analysis of the effect of each model term by comparing the residuals of the Full model and the model without each term.
	The model is fitted using the train data set, but analysis of residual improvement is done on the train and test data sets.
	Residuals are compared by a paired \emph{t}-test, a paired Wilcoxon rank-sum test, a binomial sign test and the \emph{F}-test on residual variance.
	Additionally, the net residual improvement (NeRI) of each model term is reported.
}
\usage{
	getVar.Res(object,
	           data,
	           Outcome = "Class",
	           type = c("LM", "LOGIT", "COX"),
	           testData = NULL,
	           callCpp=TRUE)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
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
		is set to true it will use the c++ implementation of residual improvement.
	}
}
\value{
	\item{tP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the paired \emph{t}-test comparing the absolute values of the residuals obtained with the Full model and the model without one term
	}
	\item{BinP.value}{
		A vector in which each element represents the \emph{p}-value associated with a significant improvement in residuals according to the binomial sign test
	}
	\item{WilcoxP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the Wilcoxon rank-sum test comparing the absolute values of the residuals obtained with the Full model and the model without one term
	}
	\item{FP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the \emph{F}-test comparing the residual variances of the residuals obtained with the Full model and the model without one term
	}
	\item{NeRIs}{
		A vector in which each element represents the net residual improvement between the Full model and the model without one term
	}
	\item{testData.tP.value}{
		A vector similar to \code{tP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.BinP.value}{
		A vector similar to \code{BinP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.WilcoxP.value}{
		A vector similar to \code{WilcoxP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.FP.value}{
		A vector similar to \code{FP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.NeRIs}{
		A vector similar to \code{NeRIs}, where values were estimated in \code{testdata}
	}
	\item{unitestMSE}{
		A vector with the univariate residual mean sum of squares of each model variable on the test data
	}
	\item{unitrainMSE}{
		A vector with the univariate residual mean sum of squares of each model variable on the train data
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{getVar.Bin}}}

\keyword{Model_Inspection}
