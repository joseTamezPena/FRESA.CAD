\name{bootstrapValidation_Bin}
\alias{bootstrapValidation_Bin}
\title{Bootstrap validation of binary classification models}
\description{
	This function bootstraps the model \emph{n} times to estimate for each variable the empirical distribution of model coefficients, area under ROC curve (AUC), integrated discrimination improvement (IDI) and net reclassification improvement (NRI).
	At each bootstrap the non-observed data is predicted by the trained model, and statistics of the test prediction are stored and reported.
	The method keeps track of predictions and plots the bootstrap-validated ROC.
	It may  plots the blind test accuracy, sensitivity, and specificity, contrasted with the bootstrapped trained distributions.
}
\usage{
	bootstrapValidation_Bin(fraction = 1,
	                    loops = 200,
	                    model.formula,
	                    Outcome,
	                    data,
	                    type = c("LM", "LOGIT", "COX"),
	                    plots = FALSE,
						best.model.formula=NULL)
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
	\item{best.model.formula}{
		An object of class \code{formula} with the formula to be used for the best model
	}
}
\details{
	The bootstrap validation will estimate the confidence interval of the model coefficients and the NRI and IDI.
	The non-sampled values will be used to estimate the blind accuracy, sensitivity, and specificity.  
	A plot to monitor the evolution of the bootstrap procedure will be displayed if \code{plots} is set to TRUE.
	The plot shows the train and blind test ROC.
	The density distribution of the train accuracy, sensitivity, and specificity are also shown, with the blind test results drawn along the y-axis.
}
\value{
	\item{data}{
		The data frame used to bootstrap and validate the model
	}
	\item{outcome}{
		A vector with the predictions made by the model
	}
	\item{blind.accuracy}{
		The accuracy of the model in the blind test set
	}
	\item{blind.sensitivity}{
		The sensitivity of the model in the blind test set
	}
	\item{blind.specificity}{
		The specificity of the model in the blind test set
	}
	\item{train.ROCAUC}{
		A vector with the AUC in the bootstrap train sets
	}
	\item{blind.ROCAUC}{
		An object of class \code{roc} containing the AUC in the bootstrap blind test set
	}
	\item{boot.ROCAUC}{
		An object of class \code{roc} containing the AUC using the mean of the bootstrapped coefficients
	}
	\item{fraction}{
		The fraction of data that was sampled with replacement
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{base.Accuracy}{
		The accuracy of the original model
	}
	\item{base.sensitivity}{
		The sensitivity of the original model
	}
	\item{base.specificity}{
		The specificity of the original model
	}
	\item{accuracy}{
		A vector with the accuracies in the bootstrap test sets
	}
	\item{sensitivities}{
		A vector with the sensitivities in the bootstrap test sets
	}
	\item{specificities}{
		A vector with the specificities in the bootstrap test sets
	}
	\item{train.accuracy}{
		A vector with the accuracies in the bootstrap train sets
	}
	\item{train.sensitivity}{
		A vector with the sensitivities in the bootstrap train sets
	}
	\item{train.specificity}{
		A vector with the specificities in the bootstrap train sets
	}
	\item{s.coef}{
		A matrix with the coefficients in the bootstrap train sets
	}
	\item{boot.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing a model whose coefficients are the median of the coefficients of the bootstrapped models
	}
	\item{boot.accuracy}{
		The accuracy of the \code{mboot.model} model
	}
	\item{boot.sensitivity}{
		The sensitivity of the \code{mboot.model} model
	}
	\item{boot.specificity}{
		The specificity of the \code{mboot.model} model
	}
	\item{z.NRIs}{
		A matrix with the \emph{z}-score of the NRI for each model term, estimated using the bootstrap train sets
	}
	\item{z.IDIs}{
		A matrix with the \emph{z}-score of the IDI for each model term, estimated using the bootstrap train sets
	}
	\item{test.z.NRIs}{
		A matrix with the \emph{z}-score of the NRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.z.IDIs}{
		A matrix with the \emph{z}-score of the IDI for each model term, estimated using the bootstrap test sets
	}
	\item{NRIs}{
		A matrix with the NRI for each model term, estimated using the bootstrap test sets
	}
	\item{IDIs}{
		A matrix with the IDI for each model term, estimated using the bootstrap test sets
	}
	\item{testOutcome}{
		A vector that contains all the individual outcomes used to validate the model in the bootstrap test sets
	}
	\item{testPrediction}{
		A vector that contains all the individual predictions used to validate the model in the bootstrap test sets
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{bootstrapValidation_Res},
				 \link{plot.bootstrapValidation_Bin},
				 \link{summary.bootstrapValidation_Bin}}}

\keyword{Model_Diagnosis}