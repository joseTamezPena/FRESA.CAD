\name{RRPlot}
\alias{RRPlot}
\title{Plot and Analysis of Indices of Risk}
\description{
	Given a list of the probabilities of risks per subject, the function will plot the risk ratios as a function of Positive predictive values (PPV)
}

RRPlot(riskData=NULL,atProb=0.90,title="Relative Risk",timetoEvent=NULL,titleS="Kaplan-Meier")

\usage{
	RRPlot(riskData=NULL,
	atProb=0.90,
	title="Relative Risk",
	timetoEvent=NULL,
	titleS="Kaplan-Meier")	
}
\arguments{
	\item{riskData}{
	The data frame with two columns: Subjects with event, Probability of Future Event
	}
	\item{atProb}{
	The desired negative predictive value of the computed risk at threshold
	}
	\item{title}{
	The title of the plot
	}
	\item{timetoEvent}{
	The time to event vector
	}
	\item{titleS}{
	The title of the Kaplan-Meier plot
	}
}

\details{
 The RRPlot function will plot the Relative Risk (RR) values for every possible probability threshold of the positive cases. Furthermore, it will compute and analyze the RR of the computed threshold that contains the prescribed percentage of negative cases. If the user provides the time to event data, the function will also plot the Kaplan-Meier curve and the logrank probability of differences. between subjects at risk.
 }

\value{
	\item{EventsThr}{
	The explored thresholds
	}
	\item{isEvent}{
		The list of events per threshold values
	}
	\item{RR}{
		The RR for every possible threshold
	}
	\item{PPV}{
		The positive predictive values for every threshold
	}
	\item{NPV}{
		The negative predictive value for every threshold
	}
	\item{Sensitivity}{
		The sensitivity for every threshold
	}
	\item{fit}{
		The loess fit 
	}
	\item{thr_atP}{
		The p-value that contains atProb of the negative subjects
	}
	\item{PPV_atP}{
		The precession at the threshold
	}
	\item{LowEventsFrac_atP}{
		The fraction of subjects at low risk
	}
	\item{HighEventsFrac_atP}{
		The fraction of subjects at high risk
	}
	\item{RR_atP}{
		The estimated RR based on the loess fit
	}
	\item{surfit}{
		The survival fit object
	}
	\item{LogRankE}{
		The p-value of the logrank test
	}
}

\author{Jose G. Tamez-Pena}

\examples{
  #TBD
}

\seealso{\code{EmpiricalSurvDiff}}
\keyword{Data Visualization}


