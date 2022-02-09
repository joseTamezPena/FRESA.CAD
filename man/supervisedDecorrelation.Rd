\name{GDSMTDecorrelation}
\alias{GDSMTDecorrelation}
\alias{predictDecorrelate}
\title{Supervised decorrelation of dataframe features}
\description{
	All features significatly associated to the outcome will be decorrelated
}
\usage{
	GDSMTDecorrelation(data=NULL,thr=0.80,
	                     refdata=NULL,Outcome=NULL,
	                     baseFeatures=NULL,unipvalue=0.05,
	                     useDeCorr=TRUE,maxLoops=100,
	                     verbose=FALSE,
	                     method=c("fast","pearson","spearman","kendall"),
						 skipRelaxed=TRUE,
	                     ...)

	
	predictDecorrelate(decorrelatedobject,testData)
}
\arguments{
	\item{data}{
	The dataframe whose features will de decorrelated
	}
	\item{thr}{
	correlation threshold for feature selection
	}
	\item{refdata}{
	The frame used to get the correlation and lm formula 
	}
	\item{Outcome}{
	The target outcome
	}
	\item{baseFeatures}{
	A vector of features to be used as the base vectors. 
	}
	\item{unipvalue}{
	association p-value
	}
	\item{useDeCorr}{
	 if TRUE, the returned matrix will use the estimated decorrelation Matrix
	}
	\item{maxLoops}{
	the maxumum number of loops
	}
	\item{verbose}{
	if TRUE, it will display internal evolution of algorithm.
	}
	\item{method}{
	method pased to correlation function.
	}
	\item{skipRelaxed}{
	is set to FALSE it will use relaxed convergence
	}
	\item{...}{
	parameters passed to the adjuting method.
	}
	\item{decorrelatedobject}{
	The returned datafrmae of the GDSMTDecorrelation function
	}
	\item{testData}{
	the new dataframe to be decorrelated
	}
}

\details{
 The data-frame will be analyzed and correlated features will be decorrelated. Feature selection may be based on Outcome or by default in Top correlation
}

\value{
	\item{adjustedFrame}{
	The decorrelated data frame with the follwing attributes 
	}
	\item{attr:topFeatures)}{
		Attribute of adjustedFrame: The list of features that were decorrelated
	}
	\item{attr:TotalAdjustments}{
		Attribute of adjustedFrame: The count of how many times was decorrelted
	}
	\item{attr:DeCorrmatrix}{
		Attribute of adjustedFrame: The Decorrelation matrix for linear models
	}
	\item{attr:varincluded}{
		Attribute of adjustedFrame: The list of variables in the whitening matrix
	}
	\item{attr:baseFeatures}{
		Attribute of adjustedFrame: The list of feetures used as base features
	}
	\item{attr:useDeCorr}{
		Attribute of adjustedFrame: If TRUE the estimated DeCorrmatrix matrix was used for decorrelation
	}
	\item{attr:correlatedToBase}{
		Attribute of adjustedFrame: List of correlated featuers to baseFeatures 
	}
	\item{attr:AbaseFeatures}{
		Attribute of adjustedFrame: List of independent base features suplementing baseFeatures
	}
	\item{attr:fsocre}{
		Attribute of adjustedFrame: List of the correlation score of each feature.
	}
}



\author{Jose G. Tamez-Pena}
\seealso{\code{featureAdjustment}}
\keyword{Data Conditioning}
