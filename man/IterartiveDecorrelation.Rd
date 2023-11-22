\name{IDeA}
\alias{IDeA}
\alias{ILAA}
\alias{predictDecorrelate}
\title{Decorrelation of data frames}
\description{
  All continous features that with significant correlation will be decorrelated
}

\usage{
  ILAA(data=NULL,
                thr=0.80,
                method=c("pearson","spearman"),
                Outcome=NULL,
                drivingFeatures=NULL,
                maxLoops=100,
                verbose=FALSE,
                bootstrap=0
                )
                
  IDeA(data=NULL,thr=0.80,
                       method=c("fast","pearson","spearman","kendall"),
                       Outcome=NULL,
                       refdata=NULL,
                       drivingFeatures=NULL,
                       useDeCorr=TRUE,
                       relaxed=TRUE,
                       corRank=TRUE,
                       maxLoops=100,
                       unipvalue=0.05,
                       verbose=FALSE,
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
  \item{drivingFeatures}{
  A vector of features to be used as basis vectors. 
  }
  \item{unipvalue}{
  Maximum p-value for correlation significance
  }
  \item{useDeCorr}{
   if TRUE, the transformation matrix (UPLTM) will be computed
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
  \item{relaxed}{
  is set to TRUE it will use relaxed convergence
  }
  \item{corRank}{
  is set to TRUE it will correlation matrix to break ties.
  }
  \item{...}{
  parameters passed to the \code{featureAdjustment} function.
  }
  \item{decorrelatedobject}{
  The returned dataframe of the \code{IDeA} function
  }
  \item{testData}{
  The new dataframe to be decorrelated
  }
  \item{bootstrap}{
   If greater than 1 the number of boostrapping loops
  }
}

\details{
 The dataframe will be analyzed and significantly correlated features whose correlation 
 is larger than the user supplied threshold will be decorrelated. 
 Basis feature selection may be based on Outcome association or by an unsupervised method. 
 The default options will run the decorrelation using fast matrix operations using \code{Rfast};
 hence, Pearson correlation will be used to estimate the unit-preserving spatial transformation matrix (UPLTM).
 ILAA is a wrapper of the more comprensive IDeA method. It estimates linear transforms and allows for boosted transform estimations 
}

\value{
  \item{decorrelatedDataframe}{
  The decorrelated data frame with the follwing attributes 
  }
  \item{attr:UPLTM}{
  Attribute of decorrelatedDataframe: The Decorrelation matrix with the beta coefficients
  }
  \item{attr:fscore}{
  Attribute of decorrelatedDataframe: The score of each feature.
  }
  \item{attr:drivingFeatures}{
  Attribute of decorrelatedDataframe: The list of features used as base features for supervised basis
  }
  \item{attr:unipvalue}{
  Attribute of decorrelatedDataframe: The p-value used to check for fit significance
  }
  \item{attr:R.critical}{
  Attribute of decorrelatedDataframe: The pearson correlation critical value
  }
  \item{attr:IDeAEvolution}{
  Attribute of decorrelatedDataframe: The R measure history and the sparcity
  }
  \item{attr:VarRatio}{
  Attribute of decorrelatedDataframe: The variance ratio between the output latent variable and the observed 
  }
}

\author{Jose G. Tamez-Pena}

\examples{
  \dontrun{
  # load FRESA.CAD library
  #  library("FRESA.CAD")

  # iris data set
  data('iris')

  colors <- c("red","green","blue")
  names(colors) <- names(table(iris$Species))
  classcolor <- colors[iris$Species]

  #Decorrelating with usupervised basis and correlation goal set to 0.25
  system.time(irisDecor <- IDeA(iris,thr=0.25))
  
  ## The transformation matrix is stored at "UPLTM" attribute
  UPLTM <- attr(irisDecor,"UPLTM")
  print(UPLTM)

  #Decorrelating with supervised basis and correlation goal set to 0.25
  system.time(irisDecorOutcome <- IDeA(iris,Outcome="Species",thr=0.25))
  ## The transformation matrix is stored at "UPLTM" attribute
  UPLTM <- attr(irisDecorOutcome,"UPLTM")
  print(UPLTM)

  ## Compute PCA 
  features <- colnames(iris[,sapply(iris,is,"numeric")])
  irisPCA <- prcomp(iris[,features]);
  ## The PCA transformation
  print(irisPCA$rotation)

  ## Plot the transformed sets
  plot(iris[,features],col=classcolor,main="Raw IRIS")

  plot(as.data.frame(irisPCA$x),col=classcolor,main="PCA IRIS")

  featuresDecor <- colnames(irisDecor[,sapply(irisDecor,is,"numeric")])
  plot(irisDecor[,featuresDecor],col=classcolor,main="Outcome-Blind IDeA IRIS")


  featuresDecor <- colnames(irisDecorOutcome[,sapply(irisDecorOutcome,is,"numeric")])
  plot(irisDecorOutcome[,featuresDecor],col=classcolor,main="Outcome-Driven IDeA IRIS")
  }
}
\seealso{\code{featureAdjustment}}
\keyword{Data Transformations}


