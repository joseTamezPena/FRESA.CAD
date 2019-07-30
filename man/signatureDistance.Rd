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
	              method = c("pearson","spearman","kendall","RSS","MAN")
	)
}
\arguments{
	\item{template}{
		A list with a template matrix of the signature described with 
		quantiles = [0.025,0.25,0.5,0.75,0.975]
	}
	\item{data}{
		A data frame that will be used to compute the distance
	}
	\item{method}{
		The distance method.
	}
}
\details{
   The distance to the template:
   "pearson","spearman" and "kendall" distances are computed using the correlation function i.e. 1-r. 
   "RSS" distance is the root sum square distance
   "MAN" Manhattan. The standardized L^1 distance  
}
\value{
	\item{result}{
		the distance to the template
	}
}
\author{Jose G. Tamez-Pena}
\keyword{Feature_Selection}