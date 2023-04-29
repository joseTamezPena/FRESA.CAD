\name{RRPlot}
\alias{RRPlot}
\title{Plot and Analysis of Indices of Risk}
\description{
	Plots of calibration and performance of risk probabilites 
}

\usage{
	RRPlot(riskData=NULL,
	timetoEvent=NULL,
	riskTimeInterval=NULL,
	ExpectedPrevalence=NULL,
	atProb=c(0.90,0.80),
	atThr=NULL,
    plotRR=TRUE,
	title="",
	ysurvlim=c(0,1.0)
	)	
}
\arguments{
	\item{riskData}{
	The data frame with two columns: First: Event label (event=1, censored=0). Second: Probability of any future event within the riskTimeInterval
	}
	\item{timetoEvent}{
	The time to event vector
	}
	\item{riskTimeInterval}{
	The time interval of the probability estimations
	}
	\item{ExpectedPrevalence}{
	For Case-Control Studies: The expected prevalence of events.
	}
	\item{atProb}{
	The desired specificity or sensitivity of the computed risk at threshold
	}
	\item{atThr}{
	The risk threshold
	}
	\item{plotRR}{
	If set to FALSE it will not generate the plots
	}
	\item{title}{
	The title postfix to be appended on each one of the generated plot titles
	}
	\item{ysurvlim}{
	The y limits of the survival plot
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
	\item{surdif}{
		The logrank test analysis
	}
	\item{LogRankE}{
		The bootstreped p-value of the logrank test
	}
}

\author{Jose G. Tamez-Pena}

\examples{
	\dontrun{

        ### RR Plot Example ####
        # Start the graphics device driver to save all plots in a pdf format
        pdf(file = "RRPlot.pdf",width = 8, height = 6)
		library(survival)
		library(FRESA.CAD)

		op <- par(no.readonly = TRUE)

		## The lung data from the survival package
		data(cancer, package="survival")
		lungD <- lung
		lungD$inst <- NULL
		lungD$status <- lungD$status - 1
		lungD <- lungD[complete.cases(lungD),]


		## Exploring Raw Features with RRPlot

		convar <- colnames(lungD)[lapply(apply(lungD,2,unique),length) > 10]
		convar <- convar[convar != "time"]
		topvar <- univariate_BinEnsemble(lungD[,c("status",convar)],"status")
		print(names(topvar))
		topv <- min(5,length(topvar))
		topFive <- names(topvar)[1:topv]
		RRanalysis <- list();
		idx <- 1
		for (topf in topFive)
		{
		  RRanalysis[[idx]] <- RRPlot(cbind(lungD$status,lungD[,topf]),
									  atProb=c(0.90),
									  timetoEvent=lungD$time,
									  title=topf
									  #plotRR=FALSE
		  )
		  idx <- idx + 1
		}
		names(RRanalysis) <- topFive

		## Reporting the Metrics

		ROCAUC <- NULL
		CstatCI <- NULL
		RRatios <- NULL
		LogRangp <- NULL
		Sensitivity <- NULL
		Specificity <- NULL

		for (topf in topFive)
		{
		  CstatCI <- rbind(CstatCI,RRanalysis[[topf]]$c.index$cstatCI)
		  RRatios <- rbind(RRatios,RRanalysis[[topf]]$RR_atP)
		  LogRangp <- rbind(LogRangp,RRanalysis[[topf]]$surdif$pvalue)
		  Sensitivity <- rbind(Sensitivity,RRanalysis[[topf]]$ROCAnalysis$sensitivity)
		  Specificity <- rbind(Specificity,RRanalysis[[topf]]$ROCAnalysis$specificity)
		  ROCAUC <- rbind(ROCAUC,RRanalysis[[topf]]$ROCAnalysis$aucs)
		}
		rownames(CstatCI) <- topFive
		rownames(RRatios) <- topFive
		rownames(LogRangp) <- topFive
		rownames(Sensitivity) <- topFive
		rownames(Specificity) <- topFive
		rownames(ROCAUC) <- topFive

		print(ROCAUC)
		print(CstatCI)
		print(RRatios)
		print(LogRangp)
		print(Sensitivity)
		print(Specificity)

		meanMatrix <- cbind(ROCAUC[,1],
   		                    CstatCI[,1],
   		                    Sensitivity[,1],
   		                    Specificity[,1],
   		                    RRatios[,1])
		colnames(meanMatrix) <- c("ROCAUC","C-Stat","Sen","Spe","RR")
		print(meanMatrix)

		## COX Modeling
		ml <- BSWiMS.model(Surv(time,status)~1,data=lungD,NumberofRepeats = 10)
		sm <- summary(ml)
		print(sm$coefficients)

		### Cox Model Performance


		timeinterval <- 2*mean(subset(lungD,status==1)$time)

		h0 <- sum(lungD$status & lungD$time <= timeinterval)
		h0 <- h0/sum((lungD$time > timeinterval) | (lungD$status==1))
		print(t(c(h0=h0,timeinterval=timeinterval)),caption="Initial Parameters")

		index <- predict(ml,lungD)

		rdata <- cbind(lungD$status,ppoisGzero(index,h0))

		rrAnalysisTrain <- RRPlot(rdata,atProb=c(0.90),
								  timetoEvent=lungD$time,
								  title="Raw Train: lungD Cancer",
								  ysurvlim=c(0.00,1.0),
								  riskTimeInterval=timeinterval)


		### Reporting Performance 


		print(t(rrAnalysisTrain$OERatio),caption="O/E Ratio")
		print(t(rrAnalysisTrain$OE95ci),caption="O/E Ratio")
		print(t(rrAnalysisTrain$OAcum95ci),caption="O/Acum Ratio")
		print(rrAnalysisTrain$c.index$cstatCI,caption="C. Index")
		print(t(rrAnalysisTrain$ROCAnalysis$aucs),caption="ROC AUC")
		print((rrAnalysisTrain$ROCAnalysis$sensitivity),caption="Sensitivity")
		print((rrAnalysisTrain$ROCAnalysis$specificity),caption="Specificity")
		print(t(rrAnalysisTrain$thr_atP),caption="Probability Thresholds")
		print(t(rrAnalysisTrain$RR_atP),caption="Risk Ratio")
		print(rrAnalysisTrain$surdif,caption="Logrank test")

  
        dev.off()
	}
}

\seealso{\code{EmpiricalSurvDiff}}
\keyword{Data Visualization}


