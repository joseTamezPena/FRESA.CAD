\name{bootstrapValidation_Res}
\alias{bootstrapValidation_Res}
\title{Bootstrap validation of regression models}
\description{
	This function bootstraps the model \emph{n} times to estimate for each variable the empirical bootstrapped distribution of model coefficients, and net residual improvement (NeRI).
	At each bootstrap the non-observed data is predicted by the trained model, and statistics of the test prediction are stores and reported.
}
\usage{
	bootstrapValidation_Res(fraction = 1,
	                        loops = 200,
	                        model.formula,
	                        Outcome,
	                        data,
	                        type = c("LM", "LOGIT", "COX"),
	                        plots = FALSE,
							bestmodel.formula=NULL)
}
\arguments{
	\item{fraction}{
		The fraction of data (sampled with replacement) to be used as train
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{plots}{
		Logical. If \code{TRUE}, density distribution plots are displayed
	}
	\item{bestmodel.formula}{
		An object of class \code{formula} with the best formula to be compared
	}
}
\details{
	The bootstrap validation will estimate the confidence interval of the model coefficients and the NeRI.
	It will also compute the train and blind test root-mean-square error (RMSE), as well as the distribution of the NeRI \emph{p}-values.
}
\value{
	\item{data}{
		The data frame used to bootstrap and validate the model
	}
	\item{outcome}{
		A vector with the predictions made by the model
	}
	\item{boot.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing a model whose coefficients are the median of the coefficients of the bootstrapped models
	}
	\item{NeRIs}{
		A matrix with the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{tStudent.pvalues}{
		A matrix with the \emph{t}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{wilcox.pvalues}{
		A matrix with the Wilcoxon rank-sum test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{bin.pvalues}{
		A matrix with the binomial test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{F.pvalues}{
		A matrix with the \emph{F}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{test.tStudent.pvalues}{
		A matrix with the \emph{t}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.wilcox.pvalues}{
		A matrix with the Wilcoxon rank-sum test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.bin.pvalues}{
		A matrix with the binomial test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.F.pvalues}{
		A matrix with the \emph{F}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{testPrediction}{
		A vector that contains all the individual predictions used to validate the model in the bootstrap test sets
	}
	\item{testOutcome}{
		A vector that contains all the individual outcomes used to validate the model in the bootstrap test sets
	}
	\item{testResiduals}{
		A vector that contains all the residuals used to validate the model in the bootstrap test sets
	}
	\item{trainPrediction}{
		A vector that contains all the individual predictions used to validate the model in the bootstrap train sets
	}
	\item{trainOutcome}{
		A vector that contains all the individual outcomes used to validate the model in the bootstrap train sets
	}
	\item{trainResiduals}{
		A vector that contains all the residuals used to validate the model in the bootstrap train sets
	}
	\item{testRMSE}{
		The global RMSE, estimated using the bootstrap test sets
	}
	\item{trainRMSE}{
		The global RMSE, estimated using the bootstrap train sets
	}
	\item{trainSampleRMSE}{
		A vector with the RMSEs in the bootstrap train sets
	}
	\item{testSampledRMSE}{
		A vector with the RMSEs in the bootstrap test sets
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{bootstrapValidation_Bin},
				 \link{plot.bootstrapValidation_Res}}}

\keyword{Model_Diagnosis}