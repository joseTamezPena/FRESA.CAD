\name{GDSTMDecorrelation}
\alias{GDSTMDecorrelation}
\alias{predictDecorrelate}
\title{Decorrelation of data frames}
\description{
	All continous features that with significant correlation will be decorrelated
}
\usage{
	GDSTMDecorrelation(data=NULL,thr=0.80,
	                     refdata=NULL,Outcome=NULL,
	                     baseFeatures=NULL,unipvalue=0.05,
	                     useDeCorr=TRUE,maxLoops=100,
	                     verbose=FALSE,
	                     method=c("fast","pearson","spearman","kendall"),
						 skipRelaxed=FALSE,
						 corRank=FALSE,
	                     ...)

	
	predictDecorrelate(decorrelatedobject,testData)
}
\arguments{
	\item{data}{
	The dataframe whose features will de decorrelated
	}
	\item{thr}{
	The maximum allowed correlation.
	}
	\item{refdata}{
	Option: A data frame that may be used to decorrelate the target dataframe 
	}
	\item{Outcome}{
	The target outcome for supervised basis
	}
	\item{baseFeatures}{
	A vector of features to be used as basis vectors. 
	}
	\item{unipvalue}{
	Maximum p-value for correlation significance
	}
	\item{useDeCorr}{
	 if TRUE, the transformation matrix (GDSTM) will be computed
	}
	\item{maxLoops}{
	the maxumum number of iteration loops
	}
	\item{verbose}{
	if TRUE, it will display internal evolution of algorithm.
	}
	\item{method}{
	if not set to "fast" the method will be pased to the \code{cor()} function.
	}
	\item{skipRelaxed}{
	is set to FALSE it will use relaxed convergence
	}
	\item{corRank}{
	is set to TRUE it will correlation matrix to break ties.
	}
	\item{...}{
	parameters passed to the \code{featureAdjustment} function.
	}
	\item{decorrelatedobject}{
	The returned dataframe of the \code{GDSTMDecorrelation} function
	}
	\item{testData}{
	the new dataframe to be decorrelated
	}
}

\details{
 The dataframe will be analyzed and significantly correlated features whose correlation 
 is larger than the user supplied threshold will be decorrelated. 
 Basis feature selection may be based on Outcome association or by an unsupervised method. 
 The default options will run the decorrelation using fast matrix operations using \code{Rfast};
 hence, Pearson correlation will be used to estimate the GDSTM. 
}

\value{
	\item{decorrelatedDataframe}{
	The decorrelated data frame with the follwing attributes 
	}
	\item{attr:topFeatures)}{
		Attribute of adjustedFrame: The list of features that were decorrelated
	}
	\item{attr:TotalAdjustments}{
		Attribute of adjustedFrame: The count of how many iteration were required
	}
	\item{attr:GDSTM}{
		Attribute of adjustedFrame: The Decorrelation matrix with the beta coefficients
	}
	\item{attr:varincluded}{
		Attribute of adjustedFrame: The list of input variables used in GDSTM
	}
	\item{attr:baseFeatures}{
		Attribute of adjustedFrame: The list of features used as base features for supervised basis
	}
	\item{attr:useDeCorr}{
		Attribute of adjustedFrame: If TRUE the estimated GDSTM is used for decorrelation
	}
	\item{attr:correlatedToBase}{
		Attribute of adjustedFrame: List of correlated features to the base features 
	}
	\item{attr:AbaseFeatures}{
		Attribute of adjustedFrame: List of unsupervised basis features
	}
	\item{attr:fscore}{
		Attribute of adjustedFrame: The score of each feature.
	}
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
	
	## The transformation matrix is stored at "GDSTM" attribute
	GDSTM <- attr(irisDecor,"GDSTM")
	print(GDSTM)

	#Decorrelating with supervised basis and correlation goal set to 0.25
	system.time(irisDecorOutcome <- GDSTMDecorrelation(iris,Outcome="Species",thr=0.25))
	## The transformation matrix is stored at "GDSTM" attribute
	GDSTM <- attr(irisDecorOutcome,"GDSTM")
	print(GDSTM)

	## Compute PCA 
	features <- colnames(iris[,sapply(iris,is,"numeric")])
	irisPCA <- prcomp(iris[,features]);
	## The PCA transformation
	print(irisPCA$rotation)

	## Plot the transformed sets
	plot(iris[,features],col=classcolor,main="Raw IRIS")

	plot(as.data.frame(irisPCA$x),col=classcolor,main="PCA IRIS")

	featuresDecor <- colnames(irisDecor[,sapply(irisDecor,is,"numeric")])
	plot(irisDecor[,featuresDecor],col=classcolor,main="Unsupervised FCA IRIS")


	featuresDecor <- colnames(irisDecorOutcome[,sapply(irisDecorOutcome,is,"numeric")])
	plot(irisDecorOutcome[,featuresDecor],col=classcolor,main="Supervised FCA IRIS")
}

\seealso{\code{featureAdjustment}}
\keyword{Data Transformations}


