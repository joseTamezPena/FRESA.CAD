\name{FRESAScale}
\alias{FRESAScale}
\title{Data frame normalization}
\description{
	All features from the \code{data} will be normalized based on the distribution of the reference data-frame
}

\usage{
	FRESAScale(data,refFrame=NULL,method=c("Norm","Order","OrderLogit","RankInv"),
    refMean=NULL,refDisp=NULL,strata=NA)
}
\arguments{
	\item{data}{
	The dataframe to be normalized
	}
	\item{refFrame}{
	The reference frame that will be used to extract the feature distribution
	}
	\item{method}{
	The normalization method. Norm: Mean and Std, Order: Median and IQR,OrderLogit order plus logit, RankInv: \code{\link{rankInverseNormalDataFrame}} 
	}
	\item{refMean}{
	The mean vector of the reference frame
	}
	\item{refDisp}{
	the data dispersion method of the reference frame
	}
	\item{strata}{
	the data stratification variable for the RankInv method
	}
}

\details{
 The data-frame will be normalized according to the distribution of the reference frame or the mean vector(\code{refMean}) scaled by the reference dispersion vector(\code{refDisp}).
}

	
\value{
	\item{scaledData}{
	The scaled data set
	}
	\item{refMean}{
	The mean or median vector of the reference frame
	}
	\item{refDisp}{
	The data dispersion (standard deviation or IQR)
	}
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{rankInverseNormalDataFrame}}}
\keyword{Model_CV}
