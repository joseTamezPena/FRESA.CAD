\name{uniRankVar}
\alias{uniRankVar}
\title{Univariate analysis of features (additional values returned)}
\description{
	This function reports the mean and standard deviation for each feature in a model, and ranks them according to a user-specified score.
	Additionally, it does a Kolmogorov-Smirnov (KS) test on the raw and \emph{z}-standardized data.
	It also reports the raw and \emph{z}-standardized \emph{t}-test score, the \emph{p}-value of the Wilcoxon rank-sum test, the integrated discrimination improvement (IDI), the net reclassification improvement (NRI), the net residual improvement (NeRI), and the area under the ROC curve (AUC).
	Furthermore, it reports the \emph{z}-value of the variable significance	on the fitted model.
	Besides reporting an ordered data frame, this function returns all arguments as values, so that the results can be updates with the \code{update.uniRankVar} if needed.
}
\usage{
	uniRankVar(variableList,
	           formula,
	           Outcome,
	           data,
	           categorizationType = c("Raw",
	                                   "Categorical",
	                                   "ZCategorical",
	                                   "RawZCategorical",
	                                   "RawTail",
	                                   "RawZTail",
	                                   "Tail",
	                                   "RawRaw"),
	           type = c("LOGIT", "LM", "COX"),
	           rankingTest = c("zIDI",
	                           "zNRI",
	                           "IDI",
	                           "NRI",
	                           "NeRI",
	                           "Ztest",
	                           "AUC",
	                           "CStat",
	                           "Kendall"),
	            cateGroups = c(0.1, 0.9),
	            raw.dataFrame = NULL,
	            description = ".",
	            uniType = c("Binary", "Regression"),
	            FullAnalysis=TRUE,
	            acovariates = NULL, 
	            timeOutcome = NULL)
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores an optional binary outcome that may be used to show the stratified analysis
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{categorizationType}{
		How variables will be analyzed : As given in \code{data} ("Raw"); broken into the \emph{p}-value categories given by \code{cateGroups} ("Categorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, and weighted by the \emph{z}-score ("ZCategorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, weighted by the \emph{z}-score, plus the raw values ("RawZCategorical"); raw values, plus the tails ("RawTail"); or raw values, weighted by the \emph{z}-score, plus the tails ("RawZTail")
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{rankingTest}{
		Variables will be ranked based on: The \emph{z}-score of the IDI ("zIDI"), the \emph{z}-score of the NRI ("zNRI"), the IDI ("IDI"), the NRI ("NRI"), the NeRI ("NeRI"), the \emph{z}-score of the model fit ("Ztest"), the AUC ("AUC"), the Somers' rank correlation ("Cstat"), or the Kendall rank correlation ("Kendall")
	}
	\item{cateGroups}{
		A vector of percentiles to be used for the categorization procedure
	}
	\item{raw.dataFrame}{
		A data frame similar to \code{data}, but with unadjusted data, used to get the means and variances of the unadjusted data
	}
	\item{description}{
		The name of the column in \code{variableList} that stores the variable description
	}
	\item{uniType}{
		Type of univariate analysis: Binary classification ("Binary") or regression ("Regression")
	}
	\item{FullAnalysis}{
		If FALSE it will only order the features according to its z-statistics of the linear model
	}
	\item{acovariates}{
		the list of covariates
	}
	\item{timeOutcome}{
		the name of the Time to event feature
	}
}
\details{
	This function will create valid dummy categorical variables if, and only if, \code{data} has been \emph{z}-standardized.
	The \emph{p}-values provided in \code{cateGroups} will be converted to its corresponding \emph{z}-score, which will then be used to create the categories.
	If non \emph{z}-standardized data were to be used, the categorization analysis would return wrong results.
}
\value{
	\item{orderframe}{
		A sorted list of model variables stored in a data frame 
	}
	\item{variableList}{
		The argument \code{variableList}
	}
	\item{formula}{
		The argument \code{formula}
	}
	\item{Outcome}{
		The argument \code{Outcome}
	}
	\item{data}{
		The argument \code{data}
	}
	\item{categorizationType}{
		The argument \code{categorizationType}
	}
	\item{type}{
		The argument \code{type}
	}
	\item{rankingTest}{
		The argument \code{rankingTest}
	}
	\item{cateGroups}{
		The argument \code{cateGroups}
	}
	\item{raw.dataFrame}{
		The argument \code{raw.dataFrame}
	}
	\item{description}{
		The argument \code{description}
	}
	\item{uniType}{
		The argument \code{uniType}
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{update.uniRankVar}, 
				 \link{univariateRankVariables}}}

\keyword{Data_Inspection}