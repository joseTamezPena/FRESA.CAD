\name{ensemblePredict}
\alias{ensemblePredict}
\title{The median prediction from a list of models}
\description{
	Given a list of model formulas, this function will train such models and return the a single(ensemble) prediction from the list of formulas on a test data set.
	It may also provides a \emph{k}-nearest neighbors (KNN) prediction using the features listed in such models.
}
\usage{
	ensemblePredict(formulaList,
	              trainData,
	              testData = NULL, 
	              predictType = c("prob", "linear"),
	              type = c("LOGIT", "LM", "COX","SVM"),
	              Outcome = NULL,
	              nk = 0
				  )
}
\arguments{
	\item{formulaList}{
		A list made of objects of class \code{formula}, each representing a model formula to be fitted and predicted with
	}
	\item{trainData}{
		A data frame with the data to train the model, where all variables are stored in different columns
	}
	\item{testData}{
		A data frame similar to \code{trainData}, but with the data set to be predicted. If \code{NULL}, \code{trainData} will be used
	}
	\item{predictType}{
		Prediction type: Probability ("prob") or linear predictor ("linear")
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{nk}{
		The number of neighbors used to generate the KNN classification. If zero, \emph{k} is set to the square root of the number of cases. If less than zero, it will not perform the KNN classification
	}
}
\value{
	\item{ensemblePredict}{
		A vector with the median prediction for the \code{testData} data set, using the models from \code{formulaList}
	}
	\item{medianKNNPredict}{
		A vector with the median prediction for the \code{testData} data set, using the KNN models
	}
	\item{predictions}{
		A matrix, where each column represents the predictions made with each model from \code{formulaList}
	}
	\item{KNNpredictions}{
		A matrix, where each column represents the predictions made with a different KNN model
	}
	\item{wPredict}{
		A vector with the weighted mean ensemble
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Model_Inspection}
