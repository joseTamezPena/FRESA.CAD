\name{nearestNeighborImpute}
\alias{nearestNeighborImpute}
\title{nearest neighbor NA imputation}
\description{
	The function will replace any NA present in the data-frame with the median values of the nearest neighbours.
}
\usage{
	nearestNeighborImpute(tobeimputed,
	                      referenceSet=NULL,
						  catgoricCol=NULL,
	                      distol=1.05,
						  useorder=TRUE
	                     )

}
\arguments{
	\item{tobeimputed}{
		a data frame with missing values (NA values)
	}
	\item{referenceSet}{
		An optional data frame with a set of complete observations. This data frame will be added to the search set
	}
	\item{catgoricCol}{
		An optional list of columns names that should be consider categorical
	}
	\item{distol}{
		The tolerance used to define if a particular set of row observations is similar to the minimum distance 
		
	}
	\item{useorder}{
		Impute using the last observation on startified by categorical data	
	}
}
\details{
 This function will find any NA present in the data set and it will search for the row set of complete observations that have the closest IQR normalized Manhattan distance to the row with missing values. 
 If a set of rows have similar minimum distances (toldis*(minimum distance) > row set distance) the median value will be used.
}
\value{
	A data frame, where each NA has been replaced with the value of the nearest neighbors
}

\examples{
	\dontrun{
	# Get the stage C prostate cancer data from the rpart package
	library(rpart)
	data(stagec)
	# Set the options to keep the na
	options(na.action='na.pass')
	# create a model matrix with all the NA values imputed
	stagecImputed <- nearestNeighborImpute(model.matrix(~.,stagec)[,-1])
	}
}


\author{Jose G. Tamez-Pena}

\keyword{Data_Conditioning}