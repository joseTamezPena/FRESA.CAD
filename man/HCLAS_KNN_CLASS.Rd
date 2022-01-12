\name{HLCM}
\alias{HLCM}
\alias{HLCM_EM}
\title{Latent class based modeling of binary outcomes}
\description{
	Modeling a binary outcome via the the discovery of latent clusters. Each discovered latent cluster is modeled by the user provided fit function. Discovered clusters will be modeled by KNN or SVM.
}
\usage{
	HLCM(formula = formula, 
	                data=NULL,
	                method=BSWiMS.model,
	                hysteresis = 0.25,
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
	the data to be used for training the method
	}
	\item{method}{
	the binary classification function
	}
	\item{hysteresis}{
	the hysteresis shift for detecting wrongly classified subjects
	}
	\item{classMethod}{
	the function name for modeling the discovered latent clusters
	}
	\item{classModel.Control}{
	the parameters to be passed to the latent-class fitting function
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
	The original model trained with all the dataset
	}
	\item{alternativeModel}{
	The model used to classify the wrongly classified samples
	}
	\item{classModel}{
	The method that models the latent class
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
	The discovered class label of each sample
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{class::knn}} 

\keyword{Model_Generation}
