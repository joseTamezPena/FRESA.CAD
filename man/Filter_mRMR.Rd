\name{mRMR.classic_FRESA}
\alias{mRMR.classic_FRESA}
\title{FRESA.CAD wrapper of mRMRe::mRMR.classic}
\description{
	Returns the positive MI-scored set of maximum relevance minimum redundancy (mRMR) features returned by the mMRM.classic function
}
\usage{

mRMR.classic_FRESA(data=NULL, Outcome=NULL,feature_count=0,...)

}
\arguments{
	\item{data}{
	The data frame
	}
	\item{Outcome}{
	The outcome feature
	}
	\item{feature_count}{
	The number of features to return
	}
	\item{...}{
	Extra parameters to be passed to the \code{mRMRe::mMRM.classic} function
	}
}

\value{
Named vector with the MI-score of the selected features
}
\author{Jose G. Tamez-Pena}
\seealso{\code{mRMRe::mRMR.classic}}

\keyword{Feature_Filtering}
