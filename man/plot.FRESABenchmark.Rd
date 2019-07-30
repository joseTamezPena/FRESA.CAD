\name{plot.FRESA_benchmark}
\alias{plot.FRESA_benchmark}
\title{Plot the results of the model selection benchmark}
\description{
	The different output metrics of the benchmark (BinaryBenchmark,RegresionBenchmark or OrdinalBenchmark) are plotted.
	It returns data matrices that describe the different plots.
}
\usage{
	\method{plot}{FRESA_benchmark}(x,...) 
}
\arguments{
	\item{x}{
		A \code{FRESA_benchmark} object
	}
	\item{...}{
		Additional parameters for the generic \code{plot} function
	}
}
\value{
  \item{metrics}{
		The model test performance based on the \code{predictionStats_binary}, \code{predictionStats_regression} or \code{predictionStats_ordinal} functions.
	}
  \item{barPlotsCI}{
		The \code{barPlotCiError} outputs for each metric.
	}
  \item{metrics_filter}{
		The model test performance for each filter method based on the \code{predictionStats_binary} function.
	}
  \item{barPlotsCI_filter}{
		The \code{barPlotCiError} outputs for each metric on the filter methods
	}
  \item{minMaxMetrics}{
		Reports the min and maximum value for each reported metric.
	}
}

\author{Jose G. Tamez-Pena}
\seealso{\code{\link{BinaryBenchmark}}, \code{\link{predictionStats_binary}}}

\keyword{Model_Inspection}
