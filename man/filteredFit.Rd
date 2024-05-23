\name{filteredFit}
\alias{filteredFit}
\title{A generic pipeline of Feature Selection, Transformation, Scale and fit}
\description{
	Sequential application of feature selection, linear transformation, data scaling then fit
}
\usage{	
    filteredFit(formula = formula, 
	            data=NULL,
	            filtermethod=univariate_KS,
	            filtermethod.control=list(limit=0),
				Transf=c("none","PCA","CCA","ILAA"),
				Transf.control=list(thr=0.8),
	            Scale="none",
				Scale.control=list(strata=NA),
	            fitmethod=e1071::svm,
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
	\item{filtermethod.control}{
	the set of parameters required by the feature selection function
	}
	\item{Scale}{
	Scale the data using the provided method
	}
	\item{Scale.control}{
	Scale parameters
	}
	\item{Trans}{
	Data transformations: "none","PCA","CCA" or "ILAA",
	}
	\item{Trans.control}{
	Parameters to the transformation function
	}
	\item{fitmethod}{
	The fit function to be used
	}
	\item{...}{
	Parameters for the fitting function
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
