\name{predictionStats}
\alias{predictionStats_binary}
\alias{predictionStats_regression}
\alias{predictionStats_ordinal}
\alias{predictionStats_survival}
\title{Prediction Evaluation}
\description{
	This function returns the statistical metrics describing the association between model predictions and the ground truth outcome 
}
\usage{

      predictionStats_binary(predictions, plotname="", center=FALSE,...)
      predictionStats_regression(predictions, plotname="",...)
      predictionStats_ordinal(predictions,plotname="",...)
      predictionStats_survival(predictions,plotname="",...)

}
\arguments{
	\item{predictions}{
		A matrix whose first column is the ground truth, and the second is the model prediction
	}
	\item{plotname}{
		The main title to be used by the plot function. If empty, no plot will be provided 
	}
	\item{center}{
		For binary predictions indicates if the prediction is around zero 
	}
	\item{...}{
		Extra parameters to be passed to the plot function.
	}
}
\details{
	These functions will analyze the prediction outputs and will compare to the ground truth.
	The output will depend on the prediction task: Binary classification, Linear Regression, Ordinal regression or Cox regression.
}

\value{
  \item{accc}{
		The classification accuracy with its95\% confidence intervals (95/%CI)
	}
  \item{berror}{
		The balanced error rate with its 95\%CI
	}
  \item{aucs}{
		The ROC area under the curve (ROC AUC) of the binary classifier with its 95\%CI
	}
  \item{specificity}{
		The specificity with its 95\%CI
	}
  \item{sensitivity}{
		The sensitivity with its 95\%CI
	}
  \item{ROC.analysis}{
		The output of the ROC function
	}
  \item{CM.analysis}{
		The output of the \code{epiR::epi.tests} function
	}
  \item{corci}{
		the Pearson correlation with its 95\%CI
    }
  \item{biasci}{
		the regression bias and its 95\%CI
	}
  \item{RMSEci}{
		the root mean square error (RMSE) and its 95\%CI
	}
  \item{spearmanci}{
		the Spearman correlation and its 95\%CI
	}
  \item{MAEci}{
		the mean absolute difference(MAE) and its 95\%CI
	}
  \item{pearson}{
		the output of the \code{cor.test} function
	}
  \item{Kendall}{
		the Kendall correlation and its 95\%CI
	}
  \item{Bias}{
		the ordinal regression bias and its 95\%CI
	}
  \item{BMAE}{
		the balanced mean absolute difference for ordinal regression
	}
  \item{class95ci}{
		the output of the bootstrapped estimation of accuracy, sensitivity, and ROC AUC
	}
  \item{KendallTauB}{
		the output of the \code{DescTools::KendallTauB} function
	}
  \item{Kappa.analysis}{
		the output of the \code{irr::kappa2} function
	}
	\item{CIFollowUp}{
	The follow-up concordance index with its95\% confidence intervals (95/%CI)
	}
  \item{CIRisk}{
	The risks concordance index with its95\% confidence intervals (95/%CI)
	}
  \item{LogRank}{
	The LogRank test with its95\% confidence intervals (95/%CI) 
	}
}

\author{Jose G. Tamez-Pena}
\seealso{\code{\link{randomCV}}}

\keyword{Model_Inspection}
