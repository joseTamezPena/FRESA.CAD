\name{rankInverseNormalDataFrame}
\alias{rankInverseNormalDataFrame}
\title{rank-based inverse normal transformation of the data}
\description{
	This function takes a data frame and a reference control population to return a \emph{z}-transformed data set conditioned to the reference population.
	Each sample data for each feature column in the data frame is conditionally \emph{z}-transformed using a rank-based inverse normal transformation, based on the rank of the sample in the reference frame.
}
\usage{
	rankInverseNormalDataFrame(variableList,
	                           data,
	                           referenceframe,
	                           strata=NA)
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{referenceframe}{
		A data frame similar to \code{data}, but with only the control population
	}
	\item{strata}{
		The name of the column in \code{data} that stores the variable that will be used to stratify the model
	}
}
\value{
	A data frame where each observation has been conditionally \emph{z}-transformed, given control data
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\examples{
	\dontrun{
	# Start the graphics device driver to save all plots in a pdf format
	pdf(file = "Example.pdf")
	# Get the stage C prostate cancer data from the rpart package
	library(rpart)
	data(stagec)
	# Split the stages into several columns
	dataCancer <- cbind(stagec[,c(1:3,5:6)],
	                    gleason4 = 1*(stagec[,7] == 4),
	                    gleason5 = 1*(stagec[,7] == 5),
	                    gleason6 = 1*(stagec[,7] == 6),
	                    gleason7 = 1*(stagec[,7] == 7),
	                    gleason8 = 1*(stagec[,7] == 8),
	                    gleason910 = 1*(stagec[,7] >= 9),
	                    eet = 1*(stagec[,4] == 2),
	                    diploid = 1*(stagec[,8] == "diploid"),
	                    tetraploid = 1*(stagec[,8] == "tetraploid"),
	                    notAneuploid = 1-1*(stagec[,8] == "aneuploid"))
	# Remove the incomplete cases
	dataCancer <- dataCancer[complete.cases(dataCancer),]
	# Load a pre-established data frame with the names and descriptions of all variables
	data(cancerVarNames)
	# Set the group of no progression
	noProgress <- subset(dataCancer,pgstat==0)
	# z-transform g2 values using the no-progression group as reference
	dataCancerZTransform <- rankInverseNormalDataFrame(variableList = cancerVarNames[2,],
	                                                   data = dataCancer,
	                                                   referenceframe = noProgress)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Data_Conditioning}
