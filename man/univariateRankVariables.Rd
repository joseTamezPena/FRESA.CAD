\name{univariateRankVariables}
\alias{univariateRankVariables}
\title{Univariate analysis of features}
\description{
	This function reports the mean and standard deviation for each feature in a model, and ranks them according to a user-specified score.
	Additionally, it does a Kolmogorov-Smirnov (KS) test on the raw and \emph{z}-standardized data.
	It also reports the raw and \emph{z}-standardized \emph{t}-test score, the \emph{p}-value of the Wilcoxon rank-sum test, the integrated discrimination improvement (IDI), the net reclassification improvement (NRI), the net residual improvement (NeRI), and the area under the ROC curve (AUC).
	Furthermore, it reports the \emph{z}-value of the variable significance	on the fitted model.
}
\usage{
	univariateRankVariables(variableList,
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
	                        uniType = c("Binary","Regression"),
	                        FullAnalysis=TRUE,
	                        acovariates = NULL,
	                        timeOutcome = NULL
)
}
\arguments{
	\item{variableList}{
		A data frame with the candidate variables to be ranked
	}
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{categorizationType}{
		How variables will be analyzed: As given in \code{data} ("Raw"); broken into the \emph{p}-value categories given by \code{cateGroups} ("Categorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, and weighted by the \emph{z}-score ("ZCategorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, weighted by the \emph{z}-score, plus the raw values ("RawZCategorical"); raw values, plus the tails ("RawTail"); or raw values, weighted by the \emph{z}-score, plus the tails ("RawZTail")
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
	A sorted data frame.
	In the case of a binary classification analysis, the data frame will have the following columns:
	\item{Name}{
		Name of the raw variable or of the dummy variable if the data has been categorized
	}
	\item{parent}{
		Name of the raw variable from which the dummy variable was created
	}
	\item{descrip}{
		Description of the parent variable, as defined in \code{description}
	}
	\item{cohortMean}{
		Mean value of the variable
	}
	\item{cohortStd}{
		Standard deviation of the variable
	}
	\item{cohortKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the variable
	}
	\item{cohortKSP}{
		Associated \emph{p}-value to the \code{cohortKSD}
	}
	\item{caseMean}{
		Mean value of cases (subjects with \code{Outcome} equal to 1)
	}
	\item{caseStd}{
		Standard deviation of cases
	}
	\item{caseKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the variable only for cases
	}
	\item{caseKSP}{
		Associated \emph{p}-value to the \code{caseKSD}
	}
	\item{caseZKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the \emph{z}-standardized variable only for cases
	}
	\item{caseZKSP}{
		Associated \emph{p}-value to the \code{caseZKSD}
	}
	\item{controlMean}{
		Mean value of controls (subjects with \code{Outcome} equal to 0)
	}
	\item{controlStd}{
		Standard deviation of controls
	}
	\item{controlKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the variable only for controls
	}
	\item{controlKSP}{
		Associated \emph{p}-value to the \code{controlsKSD}
	}
	\item{controlZKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the \emph{z}-standardized variable only for controls
	}
	\item{controlZKSP}{
		Associated \emph{p}-value to the \code{controlsZKSD}
	}
	\item{t.Rawvalue}{
		Normal inverse \emph{p}-value (\emph{z}-value) of the \emph{t}-test performed on \code{raw.dataFrame}
	}
	\item{t.Zvalue}{
		\emph{z}-value of the \emph{t}-test performed on \code{data}
	}
	\item{wilcox.Zvalue}{
		\emph{z}-value of the Wilcoxon rank-sum test performed on \code{data}
	}
	\item{ZGLM}{
		\emph{z}-value returned by the \code{lm}, \code{glm}, or \code{coxph} functions for the \code{z}-standardized variable
	}
	\item{zNRI}{
		\emph{z}-value returned by the \code{improveProb} function (\code{Hmisc} package) when evaluating the NRI
	}
	\item{zIDI}{
		\emph{z}-value returned by the \code{improveProb} function (\code{Hmisc} package) when evaluating the IDI
	}
	\item{zNeRI}{
		\emph{z}-value returned by the \code{improvedResiduals} function when evaluating the NeRI
	}
	\item{ROCAUC}{
		Area under the ROC curve returned by the \code{roc} function (\code{pROC} package)
	}
	\item{cStatCorr}{
		\emph{c} index of Somers' rank correlation returned by the \code{rcorr.cens} function (\code{Hmisc} package)
	}
	\item{NRI}{
		NRI returned by the \code{improveProb} function (\code{Hmisc} package)
	}
	\item{IDI}{
		IDI returned by the \code{improveProb} function (\code{Hmisc} package)
	}
	\item{NeRI}{
		NeRI returned by the \code{improvedResiduals} function
	}
	\item{kendall.r}{
		Kendall \eqn{\tau} rank correlation coefficient between the variable and the binary outcome
	}
	\item{kendall.p}{
		Associated \emph{p}-value to the \code{kendall.r}
	}
	\item{TstudentRes.p}{
		\emph{p}-value of the improvement in residuals, as evaluated by the paired \emph{t}-test
	}
	\item{WilcoxRes.p}{
		\emph{p}-value of the improvement in residuals, as evaluated by the paired Wilcoxon rank-sum test
	}
	\item{FRes.p}{
		\emph{p}-value of the improvement in residual variance, as evaluated by the \emph{F}-test
	}
	\item{caseN_Z_Low_Tail}{
		Number of cases in the low tail
	}
	\item{caseN_Z_Hi_Tail}{
		Number of cases in the top tail
	}
	\item{controlN_Z_Low_Tail}{
		Number of controls in the low tail
	}
	\item{controlN_Z_Hi_Tail}{
		Number of controls in the top tail
	}
	In the case of regression analysis, the data frame will have the following columns:
	\item{Name}{
		Name of the raw variable or of the dummy variable if the data has been categorized
	}
	\item{parent}{
		Name of the raw variable from which the dummy variable was created
	}
	\item{descrip}{
		Description of the parent variable, as defined in \code{description}
	}
	\item{cohortMean}{
		Mean value of the variable
	}
	\item{cohortStd}{
		Standard deviation of the variable
	}
	\item{cohortKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the variable
	}
	\item{cohortKSP}{
		Associated \emph{p}-value to the \code{cohortKSP}
	}
	\item{cohortZKSD}{
		D statistic of the KS test when comparing a normal distribution and the distribution of the \emph{z}-standardized variable
	}
	\item{cohortZKSP}{
		Associated \emph{p}-value to the \code{cohortZKSD}
	}
	\item{ZGLM}{
		\emph{z}-value returned by the glm or Cox procedure for the z-standardized variable
	}
	\item{zNRI}{
		\emph{z}-value returned by the \code{improveProb} function (\code{Hmisc} package) when evaluating the NRI
	}
	\item{NeRI}{
		NeRI returned by the \code{improvedResiduals} function
	}
	\item{cStatCorr}{
		\emph{c} index of Somers' rank correlation returned by the \code{rcorr.cens} function (\code{Hmisc} package)
	}
	\item{spearman.r}{
		Spearman \eqn{\rho} rank correlation coefficient between the variable and the outcome
	}
	\item{pearson.r}{
		Pearson \emph{r} product-moment correlation coefficient between the variable and the outcome
	}
	\item{kendall.r}{
		Kendall \eqn{\tau} rank correlation coefficient between the variable and the outcome
	}
	\item{kendall.p}{
		Associated \emph{p}-value to the \code{kendall.r}
	}
	\item{TstudentRes.p}{
		\emph{p}-value of the improvement in residuals, as evaluated by the paired \emph{t}-test
	}
	\item{WilcoxRes.p}{
		\emph{p}-value of the improvement in residuals, as evaluated by the paired Wilcoxon rank-sum test
	}
	\item{FRes.p}{
		\emph{p}-value of the improvement in residual variance, as evaluated by the \emph{F}-test
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena}

\keyword{Data_Inspection}