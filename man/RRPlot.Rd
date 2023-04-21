\name{RRPlot}
\alias{RRPlot}
\title{Plot and Analysis of Indices of Risk}
\description{
	Plots of calibration and performance of risk probabilites 
}

\usage{
	RRPlot(riskData=NULL,
	atProb=c(0.90,0.80),
	atThr=NULL,
	title="",
	timetoEvent=NULL,
	ysurvlim=c(0,1.0),
	riskTimeInterval=NULL,
	ExpectedPrevalence=NULL
	)	
}
\arguments{
	\item{riskData}{
	The data frame with two columns: First: Event label (event=1, censored=0). Second: Probability of any future event within the riskTimeInterval
	}
	\item{atProb}{
	The desired specificity or sensitivity of the computed risk at threshold
	}
	\item{atThr}{
	The risk threshold
	}
	\item{title}{
	The title to be appended on generated each plot
	}
	\item{timetoEvent}{
	The time to event vector
	}
	\item{ysurvlim}{
	The y limits of the survival plot
	}
	\item{riskTimeInterval}{
	The time interval of the probability estimations
	}
	\item{ExpectedPrevalence}{
	For Case-Control Studies: The expected prevalence of events.
	}
	
}

\details{
The RRPlot function will analyze the provided probabilities of risk and its associated events to generate calibration plots and plots of Relative Risk (RR) vs all the sensitivity values. Furthermore, it will compute and analyze the RR of the computed threshold that contains the prescribed percentage of negative cases. If the user provides the time to event data, the function will also plot the Kaplan-Meier curve and return the logrank probability of differences between risk categories. For the calibration plot it will use the user provided riskTimeInterval to get the expected number of events. If the user does not provide the riskTimeInterval the function will use the maximum time of observations with events.
 }

\value{
	\item{CumulativeOvs}{
		Matrix with the Cumulative and Observed Events
	}
	\item{OEData}{
		Matrix with the Estimated and Observed Events
	}
	\item{DCA}{
		Decision Curve Analysis data matrix
	}
	\item{RRData}{
		The risk ratios data matrix for the ploted observations
	}
	\item{OERatio}{
		The Observed/Expected ratio with 95CI
	}
	\item{fit}{
		The loess fit of the Risk Ratios
	}
	\item{prevalence}{
		The prevalence of events
	}
	\item{thr_atP}{
		The p-value that contains atProb of the negative subjects
	}
	\item{SEN_atP}{
		The observed sensitivity at threshold
	}
	\item{LowEventsFrac_atP}{
		The fraction of subjects at low risk
	}
	\item{HighEventsFrac_atP}{
		The fraction of subjects at high risk
	}
	\item{RR_atP}{
		The fraction of subjects at high risk
	}
	\item{c.index}{
		The c-index
	}
	\item{surfit}{
		The survival fit object
	}
	\item{sufdif}{
		The logrank test analysis
	}
	\item{LogRankE}{
		The bootstreped p-value of the logrank test
	}
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{EmpiricalSurvDiff}}
\keyword{Data Visualization}


