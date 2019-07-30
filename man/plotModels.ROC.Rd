\name{plotModels.ROC}
\alias{plotModels.ROC}
\title{Plot test ROC curves of each cross-validation model}
\description{
	This function plots test ROC curves of each model found in the cross validation process.
	It will also aggregate the models into a single prediction performance, plotting the resulting ROC curve (models coherence).
	Furthermore, it will plot the mean sensitivity for a given set of specificities.
}
\usage{
   plotModels.ROC(modelPredictions,
    number.of.models=0,
    specificities=c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),
    theCVfolds=1,
    predictor="Prediction",
	cex=1.0,
    ...)
}
\arguments{
	\item{modelPredictions}{
		A data frame returned by the \code{crossValidationFeatureSelection_Bin} function, either the \code{Models.testPrediction}, the \code{FullBSWiMS.testPrediction},\cr the \code{Models.CVtestPredictions}, the \code{TestRetrained.blindPredictions},\cr the \code{KNN.testPrediction}, or the \code{LASSO.testPredictions} value
	}
	\item{number.of.models}{
		The maximum number of models to plot
	}
	\item{specificities}{
		Vector containing the specificities at which the ROC sensitivities will be calculated
	}
	\item{theCVfolds}{
		The number of folds performed in a Cross-validation experiment
	}
	\item{predictor}{
		The name of the column to be plotted
	}
	\item{cex}{
		Controlling the font size of the text inside the plots
	}
	\item{...}{
		Additional parameters for the \code{roc} function (\code{pROC} package)
	}
}
\value{
	\item{ROC.AUCs}{
		A vector with the AUC of each ROC
	}
	\item{mean.sensitivities}{
		A vector with the mean sensitivity at the specificities given by \code{specificities}
	}
	\item{model.sensitivities}{
		A matrix where each row represents the sensitivity at the specificity given by \code{specificities} for a different ROC
	}
	\item{specificities}{
		The specificities used to calculate the sensitivities
	}
	\item{senAUC}{
		The AUC of the ROC curve that resulted from using \code{mean.sensitivities}
	}
	\item{predictionTable}{
		The confusion matrix between the outcome and the ensemble prediction
	}
	\item{ensemblePrediction}{
		The ensemble (median prediction) of the repeated predictions
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Model_Inspection}
