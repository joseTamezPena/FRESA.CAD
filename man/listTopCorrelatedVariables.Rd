\name{listTopCorrelatedVariables}
\alias{listTopCorrelatedVariables}
\title{List the variables that are highly correlated with each other}
\description{
	This function computes the Pearson, Spearman, or Kendall correlation for each specified variable in the data set and returns a list of the variables that are correlated to them.
	It also provides a short variable list without the highly correlated variables.
}
\usage{
	listTopCorrelatedVariables(variableList,
	                           data,
	                           pvalue = 0.001,
	                           corthreshold = 0.9,
	                           method = c("pearson", "kendall", "spearman"))
}
\arguments{
	\item{variableList}{ 
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to \code{method}, allowed for a pair of variables to be defined as significantly correlated
	}
	\item{corthreshold}{
		The minimum correlation score, associated to \code{method}, allowed for a pair of variables to be defined as significantly correlated
	}
	\item{method}{
		Correlation method: Pearson product-moment ("pearson"), Spearman's rank ("spearman"), or Kendall rank ("kendall")
	}
}
\value{
	\item{correlated.variables}{
		A data frame with two columns:\enumerate{\item cor.var.names: The variables that are correlated	\item cor.var.value: The correlation value}
	}
	\item{short.list}{
		A vector with a list of variables that are not correlated to each other. For every correlated pair, only the variable that first entered the correlation analysis was kept
	}
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
	# Load a pre-stablished data frame with the names and descriptions of all variables
	data(cancerVarNames)
	# Get the variables that have a correlation coefficient larger 
	# than 0.65 at a p-value of 0.05
	cor <- listTopCorrelatedVariables(variableList = cancerVarNames,
	                                  data = dataCancer,
	                                  pvalue = 0.05,
	                                  corthreshold = 0.65,
	                                  method = "pearson")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Data_Inspection}
