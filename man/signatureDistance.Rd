\name{signatureDistance}
\alias{signatureDistance}
\title{Distance to the signature template}
\description{
	This function returns a normalized distance to the signature template 
}
\usage{
	signatureDistance(
	              template,
	              data=NULL,
	              method = c("pearson","spearman","kendall","RSS","MAN","NB"),
	              fwts=NULL
	)
}
\arguments{
	\item{template}{
		A list with a template matrix of the signature described with 
		quantiles = [0.025,0.100,0.159,0.250,0.500,0.750,0.841,0.900,0.975]
	}
	\item{data}{
		A data frame that will be used to compute the distance
	}
	\item{method}{
		The distance method.
	}
	\item{fwts}{
		A numeric vector defining the weight of each feature
	}
}
\details{
   The distance to the template:
   "pearson","spearman" and "kendall" distances are computed using the correlation function i.e. 1-r. 
   "RSS" distance is the normalized root sum square distance
   "MAN" Manhattan. The standardized L^1 distance  
   "NB" Weighted Naive-Bayes distance  
}
\value{
	\item{result}{
		the distance to the template
	}
}
\author{Jose G. Tamez-Pena}
\keyword{Feature_Selection}