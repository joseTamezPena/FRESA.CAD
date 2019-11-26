\name{HCLAS_CLUSTER}
\alias{HCLAS_CLUSTER}
\alias{HCLAS_EM_CLUSTER}
\title{Binary Hierarchical Classification: Classifiers-CLASS-Classifiers}
\description{
	Hierarchical use of Binary Classifier. The method find clusters of subjects that can be classified by binary classifiers
}
\usage{
	HCLAS_CLUSTER(formula = formula, 
	                data=NULL,
	                method=BSWiMS.model,
	                hysteresis = 0.0,
					classMethod=KNN_method,
					classModel.Control=NULL,
					minsize=10,
	                ...
					)
}

\arguments{
	\item{formula}{
	the base formula to extract the outcome
	}
	\item{data}{
	the data to be used for training 
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
	the parameters to be passed to the class classifiers
	}
	\item{minsize}{
	the minimum size of the discovered clusters
	}
	\item{...}{
	parameters for the classification function
	}
}

\value{
	\item{original}{
	The original model trained with all the dataset or the model of the first set.
	}
	\item{alternativeModel}{
	The list of models used to classify the discovered data clusters
	}
	\item{classModel}{
	The class model that identifies the data clusters
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
