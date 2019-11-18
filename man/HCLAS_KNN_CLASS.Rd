\name{HCLAS_CLUSTER}
\alias{HCLAS_CLUSTER}
\title{Binary Hierarchical Classification: Classifiers-CLASS-Classifiers}
\description{
	Hierarchical use of Binary Classifier. Wrongly classified subjects are classified using a class method 
}
\usage{
	HCLAS_CLUSTER(formula = formula, 
	                data=NULL,
	                method=BSWiMS.model,
	                hysteresis = 0.0,
					classMethod=KNN_method,
					classModel.Control=NULL,
	                ...
					)
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
	\item{classMethod}{
	the function name of the class classifiers
	}
	\item{classModel.Control}{
	the parameters to be passed to the class-classifiers
	}
	\item{...}{
	parameters for the classification function
	}
}

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
	\item{hysteresis}{
	The used hysteresis
	}
	\item{classSet}{
	The discoverd class label of each sample
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{class::knn}} 

\keyword{Model_Generation}
