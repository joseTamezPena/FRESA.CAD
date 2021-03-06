\name{backVarElimination_Bin}
\alias{backVarElimination_Bin}
\title{IDI/NRI-based backwards variable elimination}
\description{
	This function removes model terms that do not significantly affect the integrated discrimination improvement (IDI) or the net reclassification improvement (NRI) of the model.
}
\usage{
	backVarElimination_Bin(object,
	                   pvalue = 0.05,
	                   Outcome = "Class",
	                   data,
	                   startOffset = 0,
	                   type = c("LOGIT", "LM", "COX"),
	                   selectionType = c("zIDI", "zNRI")
					   )
}
\arguments{
	\item{object}{ 
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to either IDI or NRI, allowed for a term in the model
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
	\item{selectionType}{
		The type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-score of IDI or of NRI
	}
}
\details{
	For each model term \eqn{x_i}, the IDI or NRI is computed for the Full model and the reduced model( where the term \eqn{x_i} removed).
	The term whose removal results in the smallest drop in improvement is selected. The hypothesis: the 
	term adds classification improvement is tested by checking the pvalue of improvement. If \eqn{p(IDI or NRI)>pvalue}, then the term is removed. 
	In other words, only model terms that significantly aid in subject classification are kept.
	The procedure is repeated until no term fulfils the removal criterion.
}
\value{
	\item{back.model}{
		An object of the same class as \code{object} containing the reduced model
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{reclas.info}{
		A list with the NRI and IDI statistics of the reduced model, as given by the \code{getVar.Bin} function
	}
	\item{back.formula}{
		An object of class \code{formula} with the formula used to fit the reduced model
	}
	\item{lastRemoved}{
		The name of the last term that was removed (-1 if all terms were removed)
	}
	\item{at.opt.model}{
		the model before the BH procedure
	}
	\item{beforeFSC.formula}{
		the string formula of the model before the BH procedure
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{backVarElimination_Res}, 
				 \link{bootstrapVarElimination_Bin},
				 \link{bootstrapVarElimination_Res}}}

\keyword{Model_Generation}