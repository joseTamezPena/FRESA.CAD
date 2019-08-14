\name{nearestCentroid}
\alias{nearestCentroid}
\title{Class Label Based on the Minimum Mahalanobis Distance}
\description{
	The function will return the set of labels of a data set
}
\usage{
	nearestCentroidfunction((dataset,clustermean=NULL,clustercov=NULL, p.threshold=1.0e-6)

}
\arguments{
	\item{dataset}{
		The data set to be labeled
	}
	\item{clustermean}{
		The list of cluster centers.
	}
	\item{clustercov}{
		The list of cluster covariances
	}
	\item{p.threshold}{
		The minimum aceptance p.value
	}
	}
}
\details{ The data set will be labeled based on the nearest cluster label. Points distance with membershiop probability
lower tha the aceptance threshold will have the label "zero" 
}
\value{
	\item{ClusterLabels}{
		The labels of each point
	}
}
\references{}

\examples{
	\dontrun{
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Data Labeling}