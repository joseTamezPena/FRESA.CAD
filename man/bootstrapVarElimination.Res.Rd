\name{bootstrapVarElimination_Res}
\alias{bootstrapVarElimination_Res}
\title{NeRI-based backwards variable elimination with bootstrapping}
\description{
	This function removes model terms that do not improve the bootstrapped net residual improvement (NeRI) significantly.
}
\usage{
	bootstrapVarElimination_Res(object,
	                            pvalue = 0.05,
	                            Outcome = "Class",
	                            data,
	                            startOffset = 0, 
	                            type = c("LOGIT", "LM", "COX"),
	                            testType = c("Binomial",
	                                         "Wilcox",
	                                         "tStudent",
	                                         "Ftest"),
	                            loops = 64,
	                            setIntersect = 1,
	                            print=TRUE,
	                            plots=TRUE
                                )

}	
\arguments{
	\item{object}{ 
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analysed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{startOffset}{
		Only terms whose position in the model is larger than the \code{startOffset} are candidates to be removed
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{setIntersect}{
		The intersect of the model (To force a zero intersect, set this value to 0)
	}
	\item{print}{
		Logical. If \code{TRUE}, information will be displayed
	}
	\item{plots}{
		Logical. If \code{TRUE},  plots are displayed
	}
}
\details{
	For each model term \eqn{x_i}, the residuals are computed for the Full model and the reduced model( where the term \eqn{x_i} removed).
	The term whose removal results in the smallest drop in bootstrapped test residuals improvement is selected. The hypothesis: the 
	term improves residuals is tested by checking the p-value of average improvement. If \eqn{p(residuals better than reduced residuals)>pvalue}, then the term is removed. 
	In other words, only model terms that significantly aid in improving residuals are kept.
	The procedure is repeated until no term fulfils the removal criterion.
	The p-values of improvement can be computed via a sign-test (Binomial) a paired Wilcoxon test, paired t-test or f-test. The first three tests compare the absolute values of
	the residuals, while the f-test test if the variance of the residuals is improved significantly.
}
\value{
	\item{back.model}{
		An object of the same class as \code{object} containing the reduced model
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{reclas.info}{
		A list with the NeRI statistics of the reduced model, as given by the \code{getVar.Res} function
	}
	\item{bootCV}{
		An object of class \code{bootstrapValidation_Res} containing the results of the bootstrap validation in the reduced model
	}
	\item{back.formula}{
		An object of class \code{formula} with the formula used to fit the reduced model
	}
	\item{lastRemoved}{
		The name of the last term that was removed (-1 if all terms were removed)
	}
	\item{at.opt.model}{
		The model with close to minimum bootstrapped RMSE 
	}
	\item{beforeFSC.formula}{
		The formula of the model before the FSC stage
	}
	\item{at.RMSE.formula}{
		the string formula of the model that had the minimum or close to minimum RMSE
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{bootstrapVarElimination_Bin},
				 \link{backVarElimination_Res},
				 \link{bootstrapValidation_Res}}}

\keyword{Model_Generation}