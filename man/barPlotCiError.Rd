\name{barPlotCiError}
\alias{barPlotCiError}
\title{Bar plot with error bars}
\description{
	Ranked Plot a set of measurements with error bars or confidence intervals (CI)
}
\usage{
barPlotCiError(ciTable, 
               metricname, 
			   thesets, 
			   themethod, 
			   main, 
			   angle = 0, 
			   offsets = c(0.1,0.1),
               scoreDirection = ">",
               ho=NULL,		   
			   ...)
}

\arguments{
	\item{ciTable}{
	A matrix with three columns: the value, the low CI value and the high CI value 
	}
	\item{metricname}{
	The name of the plotted values
	}
	\item{thesets}{
	A character vector with the names of the sets
	}
	\item{themethod}{
	A character vector with the names of the methods
	}
	\item{main}{
	The plot title
	}
	\item{angle}{
	The angle of the x labels
	}
	\item{offsets}{
	The offset of the x-labels
	}
	\item{scoreDirection}{
	Indicates how to aggregate the supMethod score and the ingMethod score.
	}
	\item{ho}{
	the null hypothesis
	}
	\item{...}{
	Extra parametrs pased to the barplot function
	}
}

\value{

  \item{barplot}{
		the x-location of the bars
	}
  \item{ciTable}{
		the ordered matrix with the 95 CI
	}
  \item{barMatrix}{
		the mean values of the bars
	}
  \item{supMethod}{
		A superiority score equal to the numbers of methods that were inferior
	}
  \item{infMethod}{
		A inferiority score equal to the number of methods that were superior
	}
  \item{interMethodScore}{
		the sum of supMethod and infMethod defined by the score direction.
	}

}

\author{Jose G. Tamez-Pena}

\keyword{Benchmarking}