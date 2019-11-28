\name{GLMNET}
\alias{GLMNET}
\alias{LASSO_MIN}
\alias{LASSO_1SE}
\alias{GLMNET_ELASTICNET_MIN}
\alias{GLMNET_RIDGE_MIN}
\title{GLMNET fit with feature selection"}
\description{
	Fits a \code{glmnet::cv.glmnet} object to the data, and sets the prediction to use the features that created the minimum CV error or one SE. 
}
\usage{
	GLMNET(formula = formula,data=NULL,coef.thr=0.001,s="lambda.min",...)
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
	\item{coef.thr}{
	The threshold for feature selection when alpha < 1.
	}
	\item{s}{
	The lambda threshold to be use at prediction and feature selection
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
