\name{filteredFit}
\alias{filteredFit}
\title{A generic fit method with a filtered step for feature selection}
\description{
	Sequential application of feature selection method followed by optional scaling, PCA then fitting algorithm
}
\usage{	
    filteredFit(formula = formula, 
	            data=NULL,
	            filtermethod=univariate_Wilcoxon,
	            fitmethod=e1071::svm,
	            filtermethod.control=list(pvalue=0.10,limit=0.1),
	            Scale="none",
	            PCA=FALSE,
	            DECOR=FALSE,
	            DECOR.control=list(thr=0.80,method="fast",type="NZLM"),
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
	\item{filtermethod}{
	the method for feature selection
	}
	\item{fitmethod}{
	the fit function to be used
	}
	\item{filtermethod.control}{
	the set of parameters required by the feature selection function
	}
	\item{Scale}{
	Scale the data using the provided method
	}
	\item{PCA}{
	Estimate the PCA
	}
	\item{DECOR}{
	Decorrelate the input data using FCA and the GDSTM
	}
	\item{DECOR.control}{
	Parameters to the decorrelation function
	}
	\item{...}{
	parameters for the fitting function
	}
}


\value{
	\item{fit}{
	The fitted model
	}
	\item{filter}{
	The output of the feature selection function
	}
	\item{selectedfeatures}{
	The character vector with all the selected features
	}
	\item{usedFeatures}{
	The set of features used for training
	}
	\item{parameters}{
	The parameters passed to the fitting method
	}
	\item{asFactor}{
	Indicates if the fitting was to a factor
	}
	\item{classLen}{
	The number of possible outcomes
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Model_Generation}
