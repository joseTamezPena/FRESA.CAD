\name{calBinProb}
\alias{calBinProb}
\title{Calibrates Predicted Binary Probabilities}
\description{
	The predicted binary probabilities are calibrated to match the observed event rate.
	A logistic model is used to calibrate the predicted probability to the actual event rate.
}
\usage{	
    calBinProb(BinaryOutcome=NULL, 
	            OutcomeProbability=NULL
				)
}


\arguments{
	\item{BinaryOutcome}{
	The observed binary outcome
	}
	\item{OutcomeProbability}{
	The predicted probability
	}

}


\value{
	The logistic model calibrated to the observed outcome rate
	}
	
\author{Jose G. Tamez-Pena}

\keyword{Model_Generation}
