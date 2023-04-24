\name{metric95ci}
\alias{metric95ci}
\alias{concordance95ci}
\alias{sperman95ci}
\alias{MAE95ci}
\alias{ClassMetric95ci}
\title{Estimators and 95CI}
\description{
	Bootstraped estimation of mean and 95CI 
}
\usage{
      metric95ci(metric,nss=1000,ssize=0)
      concordance95ci(datatest,nss=1000)
      sperman95ci(datatest,nss=4000)
      MAE95ci(datatest,nss=4000)
      ClassMetric95ci(datatest,nss=4000)
}
\arguments{
	\item{datatest}{
		A matrix whose first column is the model predictionground truth, and the second the prediction
	}
	\item{nss}{
		The number of bootstrap samples
	}
	\item{metric}{
		A vector with metric estimations
	}
	\item{ssize}{
		The maximim number of samples to use
	}
}
\details{
	A set of auxiliary samples to bootstrap estimations of the 95CI
}

\value{
  the mean estimation of the metrics with its corresponding  95CI
}

\author{Jose G. Tamez-Pena}
\seealso{\code{\link{randomCV}}}

\keyword{Model_Inspection}
