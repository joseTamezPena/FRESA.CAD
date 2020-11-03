\name{steerDecorrelation}
\alias{steerDecorrelation}
\title{Supervised decorrelation of dataframe features}
\description{
	All features significatly associated to the outcome will be decorrelated
}

\usage{
	steerDecorrelation(data=NULL, Outcome=NULL,refdata=NULL,loops=10,unipvalue=0.1,method=NULL,...)
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
	\item{unipvalue}{
	the maxumum p-value for linear association
	}
	\item{method}{
	fitting type of featureAdjustment
	}	
	\item{...}{
	parameters passed to the univariate filter.
	}
}

\details{
 The data-frame will be analyzed: 1) Features associated to the outcome will be extracted. 2) Features associeted to the significant features will be decorrelated
}

	
\value{
	\item{adjustedFrame}{
	The scaled data set
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{featureAdjustment}}
\keyword{Data Conditioning}
