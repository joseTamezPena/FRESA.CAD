\name{featureDecorrelation}
\alias{featureDecorrelation}
\alias{predictDecorrelate}
\title{Supervised decorrelation of dataframe features}
\description{
	All features significatly associated to the outcome will be decorrelated
}
\usage{
	featureDecorrelation(data=NULL,refdata=NULL,Outcome=NULL,loops=20,thr=0.80,unipvalue=0.05,useWhite=TRUE,...)
	predictDecorrelate(decorrelatedobject,testData)
}
\arguments{
	\item{data}{
	The dataframe whose features will de decorrelated
	}
	\item{Outcome}{
	The target outcome
	}
	\item{refdata}{
	The frame used to get the correlation and lm formula 
	}
	\item{loops}{
	the maxumum number of loops
	}
	\item{thr}{
	correlation threshold for feature selection
	}
	\item{unipvalue}{
	association p-value
	}
	\item{useWhite}{
	 if TRUE, the returned matrix will use the estimated whitening Matrix
	}
	\item{...}{
	parameters passed to the adjuting method.
	}
	\item{decorrelatedobject}{
	The dataframe created by the featureDecorrelation function
	}
	\item{testData}{
	Test data to decorrelate
	}
}

\details{
 The data-frame will be analyzed and correlated features will be decorrelated. Feature selection may be based on Outcome or by default in Top correlation
}

\value{
	\item{adjustedFrame}{
	The decorrelated data frame with the follwing attributes 
	}
	\item{topFeatures}{
		The list of features that were decorrelated
	}
	\item{TotalAdjustments}{
		The count of how many times was decorrelted
	}
	\item{DeCorrmatrix}{
		The Decorrelation matrix for linear models
	}
	\item{varincluded}{
		The list of variables in the whitening matrix
	}
	\item{baseFeatures}{
		The list of feetures used as base features
	}
	\item{useDeCorr}{
		If TRUE the estimated DeCorrmatrix matrix was used for decorrelation
	}
}

\author{Jose G. Tamez-Pena}
\seealso{\code{featureAdjustment}}
\keyword{Data Conditioning}
