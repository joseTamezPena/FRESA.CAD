\name{improvedResiduals}
\alias{improvedResiduals}
\title{Estimate the significance of the reduction of predicted residuals}
\description{
	This function will test the hypothesis that, given a set of two residuals (new vs. old), the new ones are better than the old ones as measured with non-parametric tests.
	Four \emph{p}-values are provided: one for the binomial sign test, one for the paired Wilcoxon rank-sum test, one for the paired \emph{t}-test, and one for the \code{F}-test.
	The proportion of subjects that improved their residuals, the proportion that worsen their residuals, and the net residual improvement (NeRI) will be returned.
}
\usage{
	improvedResiduals(oldResiduals,
	                  newResiduals,
	                  testType = c("Binomial", "Wilcox", "tStudent", "Ftest"))
}
\arguments{
	\item{oldResiduals}{
		A vector with the residuals of the original model
	}
	\item{newResiduals}{
		A vector with the residuals of the new model
	}
	\item{testType}{
		Type of non-parametric test to be evaluated: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
}
\details{
	This function will test the hypothesis that the new residuals are "better" than the old residuals.
	To test this hypothesis, four types of tests are performed:
	\enumerate{
		\item The paired \emph{t}-test, which compares the absolute value of the residuals
		\item The paired Wilcoxon rank-sum test, which compares the absolute value of residuals
		\item The binomial sign test, which evaluates whether the number of subjects with improved residuals is greater than the number of subjects with worsened residuals
		\item The \emph{F}-test, which is the standard test for evaluating whether the residual variance is "better" in the new residuals.
	}
	The proportions of subjects that improved and worsen their residuals are returned, and so is the NeRI.
}
\value{
	\item{p1}{
		Proportion of subjects that improved their residuals to the total number of subjects
	}
	\item{p2}{
		Proportion of subjects that worsen their residuals to the total number of subjects
	}
	\item{NeRI}{
		The net residual improvement (\code{p1}-\code{p2})
	}
	\item{p.value}{
		The one tail \emph{p}-value of the test specified in \emph{testType}
	}
	\item{BinP.value}{
		The \emph{p}-value associated with a significant improvement in residuals
	}
	\item{WilcoxP.value}{
		The single sided \emph{p}-value of the Wilcoxon rank-sum test comparing the absolute values of the new and old residuals
	}
	\item{tP.value}{
		The single sided \emph{p}-value of the paired t-test comparing the absolute values of the new and old residuals
	}
	\item{FP.value}{
		The single sided \emph{p}-value of the F-test comparing the residual variances of the new and old residuals
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}

\keyword{Model_Inspection}
