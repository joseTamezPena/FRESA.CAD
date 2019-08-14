\name{clusterISODATA}
\alias{clusterISODATA}
\title{Cluster Clustering using the Isodata Approach}
\description{
   Returns the set of Gaussian Ellipsoids that best model the data
}
\usage{
   clusterISODATA(dataset,
                 clusteringMethod=GMVECluster,
                 trainFraction=0.99,
                 randomTests=10,
                 jaccardThreshold=0.45,
                 isoDataThreshold=0.75,
                 plot=TRUE,
                 ...)

}
\arguments{
   \item{dataset}{
      The data set to be clustered
   }
   \item{clusteringMethod}{
      The clustering method.
   }
   \item{trainFraction}{
      The fraction of the data used to train the clusters
   }
   \item{randomTests}{
      The number of clustering sets that will be generated
   }
   \item{jaccardThreshold}{
      The minimum Jaccard index to be considered for data clustering
   }
   \item{isoDataThreshold}{
      The minimum distance (as p.value) between gaussian clusters
   }
   \item{plot}{
      If true it will plot the clustered points
   }
   \item{...}{
      Parameter list to be passed to the clustering method
   }
}
\details{
The data will be clustered N times as defined by number of randomTests. Once clustered the Jaccard Index map
will generated and ordered from high to low. The mean clusters parameters (Covariance sets) associated with the point with the highest Jaccard index
will define the first cluster. A cluster will be added if the mahalanobis distance between clusters is greater than the given aceptance p.value (isoDataThreshold)
Only clusters associated with points with a Jaccard index greater than jaccardThreshold will be considered.
}
\value{
   \item{cluster}{
      The numeric vector with the cluster label of each point
   }
   \item{classification}{
      The numeric vector with the cluster label of each point
   }
   \item{robustCovariance}{
      The list of robust covariances per cluster
   }
   \item{pointjaccard}{
      The mean of jaccard index per data point
   }
	\item{centers}{
		The list of cluster centers
	}
	\item{covariances}{
		The list of cluster covariance
	}
	\item{features}{
		The characer vector with the names of the features used
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Cluster_Generation}