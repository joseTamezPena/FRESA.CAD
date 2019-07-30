\name{getKNNpredictionFromFormula}
\alias{getKNNpredictionFromFormula}
\title{Predict classification using KNN}
\description{
	This function will return the classification of the samples of a test set using a \emph{k}-nearest neighbors (KNN) algorithm with euclidean distances, given a formula and a train set.
}
\usage{
	getKNNpredictionFromFormula(model.formula,
	                            trainData,
	                            testData,
	                            Outcome = "CLASS",
	                            nk = 3)
}
\arguments{
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{trainData}{
		A data frame with the data to train the model, where all variables are stored in different columns
	}
	\item{testData}{
		A data frame similar to \code{trainData}, but with the data set to be predicted
	}
	\item{Outcome}{
		The name of the column in \code{trainData} that stores the variable to be predicted by the model
	}
	\item{nk}{
		The number of neighbors used to generate the KNN classification
	}
}
\value{
	\item{prediction}{
		A vector with the predicted outcome for the \code{testData} data set
	}
	\item{prob}{
		The proportion of \emph{k} neighbors that predicted the class to be the one being reported in \code{prediction}
	}
	\item{binProb}{
		The proportion of \emph{k} neighbors that predicted the class of the outcome to be equal to 1
	}
	\item{featureList}{
		A vector with the names of the features used by the KNN procedure
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{predict.fitFRESA},
				 \link{knn}}}

\keyword{Model_Inspection}
