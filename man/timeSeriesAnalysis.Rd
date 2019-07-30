\name{timeSerieAnalysis}
\alias{timeSerieAnalysis}
\title{Fit the listed time series variables to a given model}
\description{
	This function plots the time evolution and does a longitudinal analysis of time dependent features. 
	Features listed are fitted to the provided time model (mixed effect model) with a generalized least squares (GLS) procedure.
	As output, it returns the coefficients, standard errors, \emph{t}-values, and corresponding \emph{p}-values.
}
\usage{
	timeSerieAnalysis(variableList,
	                  baseModel,
	                  data,
	                  timevar = "time",
	                  contime = ".",
	                  Outcome = ".",
	                  ...,
	                  description = ".",
	                  Ptoshow = c(1),
	                  plegend = c("p"),
	                  timesign = "-",
	                  catgo.names = c("Control", "Case")
	                  )
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{baseModel}{
		A string of the type "1 + var1 + var2" that defines the model to which variables will be fitted
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{timevar}{
		The name of the column in \code{data} that stores the visit ID
	}
	\item{contime}{
		The name of the column in \code{data} that stores the continuous time (e.g. days or months) that has elapsed since the baseline visit
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores an optional binary outcome that may be used to show the stratified analysis
	}
	\item{description}{
		The name of the column in \code{variableList} that stores the variable description
	}
	\item{Ptoshow}{
		Index of the \emph{p}-values to be shown in the plot
	}
	\item{plegend}{
		Legend of the \emph{p}-values to be shown in the plot
	}
	\item{timesign}{
		The direction of the arrow of time
	}
	\item{catgo.names}{
		The legends of the binary categories
	}
	\item{...}{
		Additional parameters to be passed to the \code{gls} function
	}
}
\details{
	This function will plot the evolution of the mean value of the listed variables with its corresponding error bars.
	Then, it will fit the data to the provided time model with a GLS procedure and it will plot the fitted values. 
	If a binary variable was provided, the plots will contain the case and control data.
	As output, the function will return the model coefficients and their corresponding \emph{t}-values, and the standard errors and their associated \emph{p}-values.
}
\value{
	\item{coef}{
		A matrix with the coefficients of the GLS fitting
	}
	\item{std.Errors}{
		A matrix with the standardized error of each coefficient
	}
	\item{t.values}{
		A matrix with the \emph{t}-value of each coefficient
	}
	\item{p.values}{
		A matrix with the \emph{p}-value of each coefficient
	}
	\item{sigmas}{
		The root-mean-square error of the fitting
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\keyword{Data_Inspection}