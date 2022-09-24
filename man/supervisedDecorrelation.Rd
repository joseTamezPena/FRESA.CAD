\name{getDerivedCoefficients}
\alias{getDerivedCoefficients}
\title{Derived Features of the GDSTM transform }
\description{
	Returs the coefficient list of derived features of the GDSTM transform
}
\usage{
	getDerivedCoefficients(decorrelatedobject)
}
\arguments{
	\item{decorrelatedobject}{
	The returned dataframe of the \code{GDSTMDecorrelation} function
	}
}

\details{
 The GDSTM transformation extracted by the \code{GDSTMDecorrelation} function is analyzed and a named list of derived features will be returned with their required formula used to compute the derived feature.
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

	colors <- c("red","green","blue")
	names(colors) <- names(table(iris$Species))
	classcolor <- colors[iris$Species]

	#Decorrelating with usupervised basis and correlation goal set to 0.25
	system.time(irisDecor <- GDSTMDecorrelation(iris,thr=0.25))
	print(getDerivedCoefficients(irisDecor));
}

\seealso{\code{GDSTMDecorrelation}}
\keyword{Data Transformations}
