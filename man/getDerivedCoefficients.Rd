\name{getLatentCoefficients}
\alias{getLatentCoefficients}
\title{Derived Features of the GDSTM transform }
\description{
	Returs the list latent features, and their corresponding coeficients, from the GDSTM transform
}
\usage{
	getLatentCoefficients(decorrelatedobject)
}
\arguments{
	\item{decorrelatedobject}{
	The returned dataframe of the \code{GDSTMDecorrelation} function
	}
}

\details{
 The GDSTM transformation extracted by the \code{GDSTMDecorrelation} function is analyzed and a named list of latent features will be returned with their required formula used to compute the latent varible.
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
	system.time(irisDecor <- GDSTMDecorrelation(iris,thr=0.25))
	print(getLatentCoefficients(irisDecor));
}

\seealso{\code{GDSTMDecorrelation}}
\keyword{Data Transformations}


