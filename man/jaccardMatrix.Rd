\name{jaccardMatrix}
\alias{jaccardMatrix}
\title{Jaccard Index of two labeled sets}
\description{
	The Jaccard Index analysis of two labeled sets
}
\usage{
	jaccardMatrix(clustersA=NULL,clustersB=NULL)

}
\arguments{
	\item{clustersA}{
		The first labeled point set
	}
	\item{clustersB}{
		The second labeled point set
	}
}
\details{
This function will compute the Jaccard Index Matrix: \eqn{[(A=i) \cap (B=j)]/[(A=i) \cup (B=j)]} for all \eqn{(i,j)} possible label pairs presenet in A and B
}
\value{
	\item{jaccardMat}{
		The numeric matrix of Jaccard Indexes of all possible paired sets 
	}
	\item{elementJaccard}{
		The corresponding Jaccard index for each data point
	}
	\item{balancedMeanJaccard}{
		The average of all marginal Jaccards
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Cluster_Evaluation}