\name{getSignature}
\alias{getSignature}
\title{Returns a CV signature template}
\description{
	This function returns the matrix template [mean,sd,IQR] that maximizes the ROC AUC between cases of controls. 
}
\usage{
	getSignature(
	              data,
	              varlist=NULL,
	              Outcome=NULL,
	              target=c("All","Control","Case"),
	              CVFolds=3,
	              repeats=9,
	              distanceFunction=signatureDistance,
	              ...
	)
}
\arguments{
	\item{data}{
		A data frame whose rows contains the sampled "subject" data, and each column is a feature.
	}
	\item{varlist}{
		The varlist is a character vector that list all the features to be searched by the Backward elimination forward selection procedure.
	}
	\item{Outcome}{
		The name of the column that has the binary outcome. 1 for cases, 0 for controls
	}
	\item{target}{
		The target template that will be used to maximize the AUC.
	}
	\item{CVFolds}{
		The number of folds to be used
	}
	\item{repeats}{
		how many times the CV procedure will be repeated
	}
	\item{distanceFunction}{
		The function to be used to compute the distance between the template and each sample
	}
	\item{...}{
		the parameters to be passed to the distance function
	}
}
\details{
   The function repeats full cycles of a Cross Validation (RCV) procedure. At each CV cycle the algorithm estimate the mean template and the distance between the template and the test samples.
The ROC AUC is computed after the RCV is completed. A forward selection scheme. The set of features that maximize the AUC during the Forward loop is returned.      
}
\value{
	\item{controlTemplate}{
		the control matrix  with quantile probs[0.025,0.25,0.5,0.75,0.975] that maximized the AUC (template of controls subjects)
	}
	\item{caseTamplate}{
		the case matrix  with quantile probs[0.025,0.25,0.5,0.75,0.975] that maximized the AUC (template of case subjects)
	}
	\item{AUCevolution}{
		The AUC value at each cycle
	}
	\item{featureSizeEvolution}{
		The number of features at each cycle
	}
	\item{featureList}{
		The final list of features
	}
	\item{CVOutput}{
		A data frame with four columns: ID, Outcome, Case Distances, Control Distances. Each row contains the CV test results 
	}
	\item{MaxAUC}{
		The maximum ROC AUC
	}
}
\author{Jose G. Tamez-Pena}
\keyword{Feature_Selection}