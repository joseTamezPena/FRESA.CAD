\name{GMVECluster}
\alias{GMVECluster}
\title{Set Clustering using the Generalized Minimum Volume Ellipsoid (GMVE)}
\description{
	The Function will return the set of Gaussian Ellipsoids that best model the data
}
\usage{
	GMVECluster(dataset, 
	            p.threshold=0.975,
	            samples=10000,
	            p.samplingthreshold=0.50,
	            sampling.rate = 3,
	            jitter=TRUE,
	            tryouts=25,
	            pca=TRUE,
	            verbose=FALSE)

}
\arguments{
	\item{dataset}{
		The data set to be clustered
	}
	\item{p.threshold}{
		The p-value threshold of point acceptance into a set.
	}
	\item{samples}{
		If the set is large, The number of random samples
	}
	\item{p.samplingthreshold}{
		Defines the maximum distance between set candidate points
	}
	\item{sampling.rate}{
		Uniform sampling rate for candidate clusters
	}
	\item{jitter}{
		If true, will jitter the data set
	}
	\item{tryouts}{
		The number of cluster candidates that will be analyed per sampled point
	}
	\item{pca}{
		If TRUE, it will use the PCA transform for dimension reduction
	}
	\item{verbose}{
		If true it will print the clustering evolution
	}
}
\details{
Implementation of the GMVE clustering algorithm as proposed by Jolion et al. (1991).
}
\value{
	\item{cluster}{
		The numeric vector with the cluster label of each point
	}
	\item{classification}{
		The numeric vector with the cluster label of each point
	}
	\item{centers}{
		The list of cluster centers
	}
	\item{covariances}{
		The list of cluster covariance
	}
	\item{robCov}{
		The list of robust covariances per cluster
	}
	\item{k}{
		The number of discovered clusters
	}
	\item{features}{
		The characer vector with the names of the features used
	}
	\item{jitteredData}{
		The jittered dataset
	}
}
\references{Jolion, Jean-Michel, Peter Meer, and Samira Bataouche. "Robust clustering with applications in computer vision." IEEE Transactions on Pattern Analysis & Machine Intelligence 8 (1991): 791-802.}

\author{Jose G. Tamez-Pena}

\keyword{Cluster_Generation}