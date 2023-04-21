\name{ppoisGzero}
\alias{ppoisGzero}
\alias{adjustProb}
\title{Poisson Progability of no Events}
\description{
 	Returs the Poisson probability of any Events 
}

\usage{
	ppoisGzero(index,h0)
	adjustProb(probGZero,gain)

}
\arguments{
	\item{index}{
	linear predictor as estimated using Cox fit.
	}	
	\item{h0}{
	The baseline hazard
	}
	\item{probGZero}{
	The probability of having any event
	}
	\item{gain}{
	The calibration gain
	}
	
}

\details{
 The ppoisGzero()function is used to get the probability of getting any number of events given the linear predictor and the baseline hazard. The adjustProb() function will return the adjusted proability of any event based on the calibration gain.
 }

\value{
	The probability of no having an event.
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{RRPlot}}
\keyword{Model Calibration}


