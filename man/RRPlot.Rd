\name{RRPlot}
\alias{RRPlot}
\title{Plot and Analysis of Indices of Risk}
\description{
	Given a list of the probabilities of risks per subject, the function will plot the risk ratios as a function of sensitivity
}

riskData=NULL,atProb=c(0.90,0.80),atThr=NULL,title="",timetoEvent=NULL,ysurvlim=c(0,1.0),riskTimeInterval=NULL

\usage{
	RRPlot(riskData=NULL,
	atProb=c(0.90,0.80),
	atThr=NULL,
	title="Relative Risk",
	timetoEvent=NULL,
	ysurvlim=c(0,1.0),
	riskTimeInterval=NULL,
	ExpectedPrevalence=NULL
	)	
}
\arguments{
	\item{riskData}{
	The data frame with two columns: Subjects with event, Probability of Future Event
	}
	\item{atProb}{
	The desired specificity or sensitivity of the computed risk at threshold
	}
	\item{atThr}{
	The risk threshold
	}
	\item{title}{
	The title of the plot
	}
	\item{timetoEvent}{
	The time to event vector
	}
	\item{ysurvlim}{
	The y limits of the survival plot
	}
	\item{riskTimeInterval}{
	The time interval of the risk estimations
	}
	\item{ExpectedPrevalence}{
	For Case-Control Studies: The expected prevalence of events.
	}
	
}

\details{
 The RRPlot function will analyze and plot the Relative Risk (RR) vs all the sensitiviy values. Furthermore, it will compute and analyze the RR of the computed threshold that contains the prescribed percentage of negative cases. If the user provides the time to event data, the function will also plot the Kaplan-Meier curve and the logrank probability of differences. between subjects at risk.
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


