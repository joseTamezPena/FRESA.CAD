\name{getLatentCoefficients}
\alias{getLatentCoefficients}
\alias{getObservedCoef}
\title{Derived Features of the UPLTM transform }
\description{
	Returs the list latent features, and their corresponding coeficients, from the UPLTM transform
}
\usage{
	getLatentCoefficients(decorrelatedobject)
    getObservedCoef(decorrelatedobject,latentModel)
}
\arguments{
	\item{decorrelatedobject}{
	The returned dataframe of the \code{IDeA} function
	}
	\item{latentModel}{
	A linear model with coefficients
	}
}

\details{
 The UPLTM transformation extracted by the \code{IDeA} function is analyzed and a named list of latent features will be returned with their required formula used to compute the latent varible.
 Given a coeficient vector of latent variables. The getObservedCoef will return a vector of coefficients associated with the observed variables.
}

\value{
	The list of derived coefficients of each one of latent feature or vector of coefficients	
}

\author{Jose G. Tamez-Pena}

\examples{

	# load FRESA.CAD library
#	library("FRESA.CAD")

# iris data set
	data('iris')


	#Decorrelating with usupervised basis and correlation goal set to 0.25
	system.time(irisDecor <- IDeA(iris,thr=0.25))
	print(getLatentCoefficients(irisDecor));
}

\seealso{\code{IDeA}}
\keyword{Data Transformations}


