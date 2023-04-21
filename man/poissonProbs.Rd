\name{ppoisGzero}
\alias{adjustProb}
\alias{ppoisGzero}
\title{Probability of more than zero events}
\description{
 	Returns the probability of having 1 or more Poisson events
}


\usage{
	ppoisGzero(index,h0)
	adjustProb(probGZero,gain)

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
}

\details{
 Auxiliary functions for the estimation of the probability of having at least one Poisson event
 }

\value{
	The probability of nozero events
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{RRPlot}}
\keyword{Model Calibration}


