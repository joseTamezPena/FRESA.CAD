\name{getMedianSurvCalibratedPrediction}
\alias{getMedianSurvCalibratedPrediction}
\alias{getMedianLogisticCalibratedPrediction}
\title{Binary Predictions Calibration of Random CV}
\description{
	Remove the bias from the test predictions generated via RandomCV
}
\usage{

      getMedianSurvCalibratedPrediction(testPredictions)
      getMedianLogisticCalibratedPrediction(testPredictions)

}
\arguments{
	\item{testPredictions}{
		A matrix with the test predictions from the randomCV() function
	}
}
\details{
	There is one function for binary predictions and one for survival predictions.
	For each trained-test prediction partition. The funciton will subtract the bias. Then it will compute the median prediction.
	Warning: This procedure is not blinded to the outcome hence it has infromation leakage.
}

\value{
  The median estimation of each calibrated predictions
}

\author{Jose G. Tamez-Pena}
\seealso{\code{\link{randomCV}}}

\keyword{Model_Calibration}
