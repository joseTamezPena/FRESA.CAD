\name{KNN_method}
\alias{KNN_method}
\title{KNN Setup for KNN prediction}
\description{
	Prepares the KNN function to be used to predict the class of a new set
}
\usage{
	KNN_method(formula = formula,data=NULL,...)
}
\arguments{
	\item{formula}{
	the base formula to extract the outcome
	}
	\item{data}{
	the data to be used for training the KNN method
	}
	\item{...}{
	parameters for the KNN function and the data scaling method
	}
}

\value{
	\item{trainData}{
	The data frame to be used to train the KNN prediction
	}
	\item{scaledData}{
	The scaled training set
	}
	\item{classData}{
	A vector with the outcome to be used by the KNN function
	}
	\item{outcome}{
	The name of the outcome
	}
	\item{usedFeatures}{
	The list of features to be used by the KNN method
	}
	\item{mean_col}{
	A vector with the mean of each training feature
	}
	\item{disp_col}{
	A vector with the dispesion of each training feature
	}
	\item{kn}{
	The number of neigbors to be used by the predict function
	}
	\item{scaleMethod}{
	The scaling method to be used by FRESAScale() function
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{class::knn},\code{\link{FRESAScale}} }

\keyword{Model_Generation}
