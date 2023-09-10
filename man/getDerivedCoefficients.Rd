\name{getLatentCoefficients}
\alias{getLatentCoefficients}
\alias{getObservedCoef}
\title{Derived Features of the UPSTM transform }
\description{
	Returs the list latent features, and their corresponding coeficients, from the UPSTM transform
}
\usage{
	getLatentCoefficients(decorrelatedobject)
    getObservedCoef(decorrelatedobject,latentCoef)
}
\arguments{
	\item{decorrelatedobject}{
	The returned dataframe of the \code{IDeA} function
	}
	\item{latentCoef}{
	A named vector of linear cofficients
	}
}

\details{
 The UPSTM transformation extracted by the \code{IDeA} function is analyzed and a named list of latent features will be returned with their required formula used to compute the latent varible.
}

\value{
	The list of derived features 
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


