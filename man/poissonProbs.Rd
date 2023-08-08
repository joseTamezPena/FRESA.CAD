\name{ppoisGzero}
\alias{adjustProb}
\alias{ppoisGzero}
\alias{meanTimeToEvent}
\alias{expectedEventsPerInterval}

\title{Probability of more than zero events}
\description{
 	Returns the probability of having 1 or more Poisson events
	the adjusted probability (adjustProb)
	the exptected time to event (meanTimeToEvent)
	or the exected number of events per interval (expectedEventsPerInterval)
}


\usage{
	ppoisGzero(index,h0)
	adjustProb(probGZero,gain)
	meanTimeToEvent(probGZero,timeInterval)
	expectedEventsPerInterval(probGZero)

}
\arguments{
	\item{index}{
	The hazard index 
	}	
	\item{h0}{
	Baseline hazard
	}
	\item{probGZero}{
	The probability of having any event
	}
	\item{gain}{
	The calibration gain
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
	Or the expected number of events per interval (expectedEventsPerInterval)
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{RRPlot}}
\keyword{Model Calibration}


