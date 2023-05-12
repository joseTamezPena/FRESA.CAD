\name{ppoisGzero}
\alias{adjustProb}
\alias{ppoisGzero}
\alias{meanTimeToEvent}

\title{Probability of more than zero events}
\description{
 	Returns the probability of having 1 or more Poisson events
}


\usage{
	ppoisGzero(index,h0)
	adjustProb(probGZero,gain)
	meanTimeToEvent(probGZero,timeInterval)

}
\arguments{
	\item{index}{
	The hazard index 
	}	
	\item{h0}{
	Baseline hazard
	}
	\item{probGZero}{
	The uncalibrated probability of nozero events
	}
	\item{gain}{
	The calibration gain
	}	
	\item{probGZero}{
	The probability of the event
	}
	\item{timeInterval}{
	The time interval 
	}
	
}

\details{
 Auxiliary functions for the estimation of the probability of having at least one Poisson event. Or the mean time to event. 
 }

\value{
	The probability of nozero events.
	Or the expected time to event (meanTimeToEvent)
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{RRPlot}}
\keyword{Model Calibration}


