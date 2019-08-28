\name{trajectoriesPolyFeatures}
\alias{trajectoriesPolyFeatures}
\title{Extract the per patient polynomial Coefficients of a feature trayectory}
\description{
	Given a longituinal data set, it will extract the associated polynomial coefficients for each sample.
}

\usage{
	trajectoriesPolyFeatures(data,
                              feature="v1", 
                              degree=2, 
                              time="t", 
                              group="ID",
                              timeOffset=0,
                              strata=NULL,
                              plot=TRUE,
                              ...)
}
\arguments{
	\item{data}{
	The dataframe 
	}
	\item{feature}{
	The name of the outcome
	}
	\item{degree}{
	The fitting function used to model the data
	}
	\item{time}{
	The percentage of the data to be used for training
	}
	\item{group}{
	The number of times that the CV process will be repeated
	}
	\item{timeOffset}{
	The time offset
	}
	\item{strata}{
	Data strafication
	}
	\item{plot}{
	if TRUE it will plot the data
	}
	\item{...}{
	parameters passed to plot
	}
}

\value{
	\item{coef}{
	The trayaectory coefficient matrix
	}
}
\author{Jose G. Tamez-Pena}


\keyword{Trajectory}
