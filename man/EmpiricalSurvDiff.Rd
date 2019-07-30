\name{EmpiricalSurvDiff}
\alias{EmpiricalSurvDiff}
\title{Estimate the LR value and its associated p-values}
\description{
	Permutations or Bootstrapping computation of the standardized log-rank (SLR) or the Chi=SLR^2 p-values for differences in survival times 
}
\usage{
	EmpiricalSurvDiff(times=times,
	                  status=status,
	                  groups=groups,
	                  samples=1000,
	                  type=c("SLR","Chi"),
	                  plots=FALSE,
	                  minAproxSamples=100,
	                  computeDist=FALSE,
	                  ...
	                  )

}
\arguments{
	\item{times}{
      A numeric vector with he observed times to event
	}
	\item{status}{
      A numeric vector indicating if the time to event is censored
	}
	\item{groups}{
      A numeric vector indicating the label of the two survival groups 
	}
	\item{samples}{
      The number of bootstrap samples
	}
	\item{type}{
      The type of log-rank statistics. SLR or Chi
	}
	\item{plots}{
      If TRUE, the Kaplan-Meier plot will be plotted 
	}
	\item{minAproxSamples}{
      The number of tail samples used for the normal-distribution approximation
	}
	\item{computeDist}{
      If TRUE, it will compute the bootstrapped distribution of the SLR
	}
	\item{...}{
      Additional parameters for the plot
	}
	
}
\details{
 It will compute the null distribution of the SRL or the square SLR (Chi) via permutations, and it will return the p-value of differences between survival times between two groups.
 It may also be used to compute the empirical distribution of the difference in SLR using bootstrapping. (computeDist=TRUE)
 The p-values will be estimated based on the sampled distribution, or normal-approximated along the tails.
}
\value{
	\item{pvalue}{
      the minimum one-tailed p-value : min[p(SRL < 0),p(SRL > 0)] for type="SLR" or the two tailed p-value: 1-p(|SRL| > 0) for type="Chi"
	}
	\item{LR}{
      A list of LR statistics: LR=Expected, VR=Variance, SLR=Standardized LR.
	}
	\item{p.equal}{
      The two tailed p-value: 1-p(|SRL| > 0)
	}
	\item{p.sup}{
      The one tailed p-value: p(SRL < 0), return NA for type="Chi"
	}
	\item{p.inf}{
      The one tailed p-value: p(SRL > 0), return NA for type="Chi"
	}
	\item{nullDist}{
      permutation derived probability density function of the null distribution
	}
	\item{LRDist}{
      bootstrapped derived probability density function of the SLR (computeDist=TRUE)
	}
}
\examples{
	\dontrun{

      library(rpart)
      data(stagec)

      # The Log-Rank Analysis using survdiff

      lrsurvdiff <- survdiff(Surv(pgtime,pgstat)~grade>2,data=stagec)
      print(lrsurvdiff)

      # The Log-Rank Analysis: permutations of the null Chi distribution
      lrp <- EmpiricalSurvDiff(stagec$pgtime,stagec$pgstat,stagec$grade>2,
                         type="Chi",plots=TRUE,samples=10000,
                         main="Chi Null Distribution")
      print(list(unlist(c(lrp$LR,lrp$pvalue))))

      # The Log-Rank Analysis: permutations of the null SLR distribution
      lrp <- EmpiricalSurvDiff(stagec$pgtime,stagec$pgstat,stagec$grade>2,
                         type="SLR",plots=TRUE,samples=10000,
                         main="SLR Null Distribution")
      print(list(unlist(c(lrp$LR,lrp$pvalue))))

      # The Log-Rank Analysis: Bootstraping the SLR distribution
      lrp <- EmpiricalSurvDiff(stagec$pgtime,stagec$pgstat,stagec$grade>2,
                         computeDist=TRUE,plots=TRUE,samples=100000,
                         main="SLR Null and SLR bootrapped")
      print(list(unlist(c(lrp$LR,lrp$pvalue))))
	
	}
}

\author{Jose G. Tamez-Pena}

\keyword{Hypothesis_Testing}
