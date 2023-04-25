\name{CalibrationProbPoissonRisk}
\alias{CoxRiskCalibration}
\alias{CalibrationProbPoissonRisk}
\title{Baseline hazard and interval time Estimations}
\description{
 	It will estimate the baseline hazard (ho) and the time interval that best describes a estimations of the probabilities of time-to-event Poisson events
}

\usage{
	CalibrationProbPoissonRisk(Riskdata,trim=0.10)
	CoxRiskCalibration(ml,data,outcome,time,trim=0.10)

}
\arguments{
	\item{Riskdata}{
	The data frame with thre columns:  Event, Probability of event, time to event
	}	
	\item{trim}{
	The percentge of tails of data not to be used to estimate the time interval
	}	
	\item{ml}{
	A Cox model of the events
	}
	\item{data}{
	the new dataframe to predict the model
	}
	\item{outcome}{
	The name of the columnt that has the event: 1 uncensored, 2; Censored
	}
	\item{time}{
	The time to event, or time to last observation.
	}
	
}

\details{
 The function will estimate the baseline hazard of Poisson events and its corresponding time interval from a list of predicted probability that the event will occur for censored (Outome=0) of the actual event happened (Outcome=1). The time to event will be used to get the best time interval that models the rate of events.
 }

\value{
	\item{index}{
		A vector with the prognistic index based on the provided probabilities
	}
	\item{probGZero}{
		The vector with the calibrated probabilites of the event happening
	}
	\item{hazard}{
		The predicted hazard of each event
	}
	\item{h0}{
		The estimated bsaeline hazard
	}
	\item{hazardGain}{
		The calibration gain
	}
	\item{timeInterval}{
		The estimated time interval of the Poisson event
	}
	\item{meaninterval}{
		The mean observed interval of events
	}
	\item{Ahazard}{
		The cumulated hazzard after calibration
	}
	\item{delta}{
		The relative difference between observed and estimated number of events.
	}
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{RRPlot}}
\keyword{Model Calibration}


