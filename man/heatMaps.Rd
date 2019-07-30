\name{heatMaps}
\alias{heatMaps}
\title{Plot a heat map of selected variables}
\description{
	This function creates a heat map for a data set based on a univariate or frequency ranking
}
\usage{
	heatMaps(variableList=NULL,
	         varRank = NULL,
	         Outcome,
	         data,
	         title = "Heat Map",
	         hCluster = FALSE,
	         prediction = NULL,
	         Scale = FALSE,
	         theFiveColors=c("blue","cyan","black","yellow","red"),
	         outcomeColors = c("blue","lightgreen","yellow","orangered","red"),
	         transpose=FALSE,
	         ...)
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{varRank}{
		A data frame with the name of the variables in \code{variableList}, ranked according to a certain metric
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{title}{
		The title of the plot
	}
	\item{hCluster}{
		Logical. If \code{TRUE}, variables will be clustered
	}
	\item{prediction}{
		A vector with a prediction for each subject, which will be used to rank the heat map
	}
	\item{Scale}{
		An optional value to force the data normalization \code{outcome}
	}
	\item{theFiveColors}{
		the colors of the heatmap 
	}
	\item{outcomeColors}{
		the colors of the outcome bar 
	}
	\item{transpose}{
		transpose the heatmap 
	}
	\item{...}{
		additional parameters for the heatmap.2 function
	}
}
\value{
	\item{dataMatrix}{
		A matrix with all the terms in \code{data} described by \code{variableList}
	}
	\item{orderMatrix}{
		A matrix similar to \code{dataMatrix}, where rows are ordered according to the \code{outcome}
	}
	\item{heatMap}{
		A list with the values returned by the \code{heatmap.2} function (\code{gplots} package)
	}
}

\examples{
	\dontrun{

		library(rpart)
		data(stagec)

		# Set the options to keep the na
		options(na.action='na.pass')
		# create a model matrix with all the NA values imputed
		stagecImputed <- as.data.frame(nearestNeighborImpute(model.matrix(~.,stagec)[,-1]))

		# the simple heat map
		hm <- heatMaps(Outcome="pgstat",data=stagecImputed,title="Heat Map",Scale=TRUE) 

		# transposing the heat-map with clustered colums
		hm <- heatMaps(Outcome="pgstat",data=stagecImputed,title="Heat Map",Scale=TRUE,
					   transpose= TRUE,hCluster = TRUE,
					   cexRow=0.80,cexCol=0.50,srtCol=35) 

		# transposing the heat-map with reds and time to event as outcome
		hm <- heatMaps(Outcome="pgtime",data=stagecImputed,title="Heat Map",Scale=TRUE,
					   theFiveColors=c("black","red","orange","yellow","white"),
					   cexRow=0.50,cexCol=0.80,srtCol=35) 
	}
}


\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Data_Inspection}
