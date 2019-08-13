\name{jaccardMatrix}
\alias{jaccardMatrix}
\title{Jaccard Index of two labeled sets}
\description{
	The Jaccard Index analysis of two labeled point sets
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
This function will compute the Jaccard Index Matrix: ((A==i) and (B==j))/((A==i) or (B==j)) for all {i,j} possible label pairs presenet in A and B
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
\references{}

\examples{
	\dontrun{
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Cluster_Evaluation}