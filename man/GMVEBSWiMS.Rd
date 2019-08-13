\name{GMVEBSWiMS}
\alias{GMVEBSWiMS}
\title{Hybrid Hierarchical Modeling with GMVE and BSWiMS}
\description{
	This function returns the BSWiMS supervised-classifier present at each one of the GMVE unsupervised Gaussian data clusters
}
\usage{
	GMVEBSWiMS(formula = formula,
	            data=NULL,
				GMVE.control = list(p.threshold = 0.85,p.samplingthreshold = 0.5),
				...,
	            )

}
\arguments{
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{GMVE.control}{
		Control parameters of the GMVECluster function 
	}
	\item{...}{
		Parameters to be passed to the BSWiMS.model function
	}
}
\details{
First, the funtion calls the BSWiMS function that return the relevant features associated with the outcome. 
Then, it calls the GMVE clustering algorithm (GMVECluster) that returns a relevant data partition based on Gaussian clusters.
Finally, the function will execute the BSWiMS.model classification function on each cluster returned by GMVECluster.
}
\value{
	\item{features}{
		The character vector with the releavant BSWiMS features.
	}
	\item{cluster}{
		The GMVECluster object
	}
	\item{models}{
		The list of BSWiMS.model models per cluster
	}
}
\references{}

\examples{
	\dontrun{
	# Get the Sonar data set
      library(mlbench)
      data(Sonar)
      Sonar$Class <- 1*(Sonar$Class == "M")
	  #Train hierachical classifier
      mc <- GMVEBSWiMS(Class~.,Sonar)
	  #report the classification
      pb <- predict(mc,Sonar)
      print(table(1*(pb>0.0),Sonar$Class))
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Model_Generation}