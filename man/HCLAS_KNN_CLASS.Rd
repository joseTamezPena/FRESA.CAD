\name{HCLAS_KNN_CLASS}
\alias{HCLAS_KNN_CLASS}
\title{Binary Hierarchical Classification: Classifiers-KNN-Classifiers}
\description{
	Hierarchical use of Binary Classifier. Wrongly classified subjects are classified using a KNN method followed by the final classifier
}
\usage{
	HCLAS_KNN_CLASS(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.0,...)
}
\arguments{
	\item{formula}{
	the base formula to extract the outcome
	}
	\item{data}{
	the data to be used for training the KNN method
	}
	\item{method}{
	the binary classification function
	}
	\item{hysteresis}{
	the hysteresis shift for detecting wrongly classified subjects
	}
	\item{...}{
	parameters for the classification function
	}
}

original = orgModel,alternativeModel = alternativeModel,classModel = classModel,accuracy=accuracy,selectedfeatures = selectedfeatures

\value{
	\item{original}{
	The original model trained with all the dataset
	}
	\item{alternativeModel}{
	The model used to classify the wrongly classified samples
	}
	\item{classModel}{
	The KNN model that identifies correct vs wrongly classified samples
	}
	\item{accuracy}{
	The original accuracy
	}
	\item{selectedfeatures}{
	The character vector of selected features
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{class::knn}} 

\keyword{Model_Generation}
