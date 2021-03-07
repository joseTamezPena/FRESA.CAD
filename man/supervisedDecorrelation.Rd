\name{featureDecorrelation}
\alias{featureDecorrelation}
\title{Supervised decorrelation of dataframe features}
\description{
	All features significatly associated to the outcome will be decorrelated
}

\usage{
	featureDecorrelation(data=NULL,refdata=NULL,Outcome=NULL,loops=c(10,10),thr=0.75,unipvalue=0.05,...)
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
	\item{...}{
	parameters passed to the adjuting method.
	}
}

\details{
 The data-frame will be analyzed and correlated features will be decorrelated. Feature selection may be based on Outcome or by default in Top correlation
}

	
\value{
	\item{adjustedFrame}{
	The scaled data set
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{featureAdjustment}}
\keyword{Data Conditioning}
