\name{LASSO}
\alias{LASSO}
\alias{LASSO_MIN}
\alias{LASSO_1SE}
\title{CV LASSO fit with s="lambda.min" or s="lambda.1se"}
\description{
	Fits a \code{glmnet::cv.glmnet} object to the data, and sets the prediction to use the features that created the minimum CV error or one SE. 
}
\usage{
	LASSO_MIN(formula = formula,data=NULL,...)
	LASSO_1SE(formula = formula,data=NULL,...)
}

\arguments{
	\item{formula}{
	The base formula to extract the outcome
	}
	\item{data}{
	The data to be used for training the KNN method
	}
	\item{...}{
	Parameters to be passed to the cv.glmnet function
	}
}

\value{
	\item{fit}{
	The \code{glmnet::cv.glmnet} fitted object
	}
	\item{s}{
	The s. Set to "lambda.min" or  "lambda.1se" for prediction
	}
	\item{formula}{
	The formula
	}
	\item{outcome}{
	The name of the outcome
	}
	\item{usedFeatures}{
	The list of features to be used
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{glmnet::cv.glmnet}}
\keyword{Model_Generation}
