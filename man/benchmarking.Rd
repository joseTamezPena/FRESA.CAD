\name{benchmarking}
\alias{BinaryBenchmark}
\alias{RegresionBenchmark}
\alias{OrdinalBenchmark}
\alias{CoxBenchmark}
\title{Compare performance of different model fitting/filtering algorithms}
\description{
	Evaluates a data set with a set of fitting/filtering methods and returns the observed cross-validation performance
}
\usage{

BinaryBenchmark(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,
                referenceCV = NULL,referenceName = "Reference"
                ,referenceFilterName="Reference")
RegresionBenchmark(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,
                referenceCV = NULL,referenceName = "Reference"
                ,referenceFilterName="Reference")
OrdinalBenchmark(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,
                referenceCV = NULL,referenceName = "Reference"
                ,referenceFilterName="Reference")

CoxBenchmark(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,
                referenceCV = NULL,referenceName = "Reference"
                ,referenceFilterName="COX.BSWiMS")
}

\arguments{
	\item{theData}{
	The data frame
	}
	\item{theOutcome}{
	The outcome feature
	}
	\item{reps}{
	The number of times that the random cross-validation will be performed
	}
	\item{trainFraction}{
	The fraction of the data used for training.
	}
	\item{referenceCV}{
	A  single random cross-validation object to be benchmarked or a list of CVObjects to be compared
	}
	\item{referenceName}{
	The name of the reference classifier to be used in the reporting tables
	}
	\item{referenceFilterName}{
	The name of the reference filter to be used in the reporting tables
	}
}
\details{
	The benchmark functions provide the performance of different classification algorithms (BinaryBenchmark), registration algorithms (RegresionBenchmark) or ordinal regression algorithms (OrdinalBenchmark)
	The evaluation method is based on applying the random cross-validation method (\code{\link{randomCV}}) that randomly splits the data into train and test sets.
	The user can provide a Cross validated object that will define the train-test partitions.

	The BinaryBenchmark compares: BSWiMS,Random Forest ,RPART,LASSO,SVM/mRMR,KNN and the ensemble of them in their ability to correctly classify the test data.
	Furthermore, it evaluates the ability of the following feature selection algorithms:
		BSWiMS or ReferenceCV,
		LASSO,
		RPART,
		RF/BSWiMS,
		IDI,
		NRI,
		t-test,
		Wilcoxon,
		Kendall,
		and mRMR
	in their ability to select the best set of features for the following classification methods:
		SVM,
		KNN,
		Naive Bayes,
		Random Forest
		Nearest Centroid (NC) with root sum square (RSS) ,
		and NC with Spearman correlation
	
	The RegresionBenchmark compares: BSWiMS,Random Forest ,RPART,LASSO,SVM/mRMR and the ensemble of them in their ability to correctly predict the test data.
	Furthermore, it evaluates the ability of the following feature selection algorithms:
		BSWiMS or referenceCV,
		LASSO,
		RPART,
		RF/BSWiMS,
		F-Test,
		W-Test,
		Pearson
		Kendall,
		and mRMR
	in their ability to select the best set of features for the following regression methods:
		Linear Regression,
		Robust Regression,
		Ridge Regression,
		LASSO,
		SVM,
		and Random Forest.

	The OrdinalBenchmark  compares: BSWiMS,Random Forest ,RPART,LASSO,KNN ,SVM and the ensemble of them in their ability to correctly predict the test data.
	Furthermore, it evaluates the ability of the following feature selection algorithms:
		BSWiMS or referenceCV,
		LASSO,
		RPART,
		RF/BSWiMS,
		F-Test,
		Kendall,
		and mRMR
	in their ability to select the best set of features for the following regression methods:
		Ordinal,
		KNN,
		SVM,
		Random Forest,
		and Naive Bayes.
		
		The CoxBenchmark compares: BSWiMS, LASSO, BeSS and Univariate Cox analysis in their ability to correctly predict the risk of event happening.
		It uses cox regression with the four alternatives, but BSWiMS, LASSO are also compared as Wrapper methods.
}

\value{
  \item{errorciTable}{
		the matrix of the balanced error with the 95 CI
	}
  \item{accciTable}{
		the matrix of the classification accuracy with the 95 CI
	}
  \item{aucTable}{
		the matrix of the ROC AUC with the 95 CI
	}
  \item{senTable}{
		the matrix of the sensitivity with the 95 CI
	}
  \item{speTable}{
		the matrix of the specificity with the 95 CI
	}
  \item{errorciTable_filter}{
		the matrix of the balanced error with the 95 CI for filter methods
	}
  \item{accciTable_filter}{
		the matrix of the classification accuracy with the 95 CI for filter methods
	}
  \item{senciTable_filter}{
		the matrix of the classification sensitivity with the 95 CI for filter methods
	}
  \item{speciTable_filter}{
		the matrix of the classification specificity with the 95 CI for filter methods
	}
  \item{aucTable_filter}{
		the matrix of the ROC AUC with the 95 CI for filter methods
	}
  \item{CorTable}{
		the matrix of the Pearson correlation with the 95 CI
	}
  \item{RMSETable}{
		the matrix of the root mean square error (RMSE) with the 95 CI
	}
  \item{BiasTable}{
		the matrix of the prediction bias with the 95 CI
	}
  \item{CorTable_filter}{
		the matrix of the Pearson correlation with the 95 CI for filter methods
	}
  \item{RMSETable_filter}{
		the matrix of the root mean square error (RMSE) with the 95 CI for filter methods
	}
  \item{BiasTable_filter}{
		the matrix of the prediction bias with the 95 CI for filter methods
	}
  \item{BMAETable}{
		the matrix of the balanced mean absolute error (MEA) with the 95 CI for filter methods
	}
  \item{KappaTable}{
		the matrix of the Kappa value with the 95 CI
	}
  \item{BiasTable}{
		the matrix of the prediction Bias with the 95 CI
	}
  \item{KendallTable}{
		the matrix of the Kendall correlation with the 95 CI
	}
  \item{MAETable_filter}{
		the matrix of the mean absolute error (MEA) with the 95 CI for filter methods
	}
  \item{KappaTable_filter}{
		the matrix of the Kappa value with the 95 CI for filter methods
	}
  \item{BiasTable_filter}{
		the matrix of the prediction Bias with the 95 CI for filter methods
	}
  \item{KendallTable_filter}{
		the matrix of the Kendall correlation with the 95 CI for filter methods
	}
	\item{CIRiskTable}{
  		the matrix of the concordance index on Risk with the 95 CI 
  }
  \item{CIFollowUpTable}{
  		the matrix of the concordance index on Follow-up times with the 95 CI 
  }
  \item{LogRankTable}{
  		the matrix of the LogRank Test with the 95 CI 
  }
  \item{CIRisksTable_filter}{
  	 	the matrix of the concordance index on Risk with the 95 CI for the filter methods
  }
  \item{CIFollowUpTable_filter}{
  	 	the matrix of the concordance index on Follow-up times with the 95 CI for the filter methods 
  }
  \item{LogRankTable_filter}{
  		the matrix of the LogRank Test with the 95 CI for the filter methods
  }
  \item{times}{
		The average CPU time used by the method
	}
  \item{jaccard_filter}{
		The average Jaccard Index of the feature selection methods
	}
  \item{TheCVEvaluations}{
		The output of the randomCV (\code{\link{randomCV}}) evaluations of the different methods
	}
  \item{testPredictions}{
		A matrix with all the test predictions
	}
  \item{featureSelectionFrequency}{
		The frequency of feature selection
	}
  \item{cpuElapsedTimes}{
		The mean elapsed times
	}
	cpuElapsedTimes
}
\seealso{\code{\link{randomCV}}}
\author{Jose G. Tamez-Pena}
\examples{
	\dontrun{

		### Binary Classification Example ####
		# Start the graphics device driver to save all plots in a pdf format
		pdf(file = "BinaryClassificationExample.pdf",width = 8, height = 6)
		# Get the stage C prostate cancer data from the rpart package

		data(stagec,package = "rpart")

		# Prepare the data. Create a model matrix without the event time
		stagec$pgtime <- NULL
		stagec$eet <- as.factor(stagec$eet)
		options(na.action = 'na.pass')
		stagec_mat <- cbind(pgstat = stagec$pgstat,
		as.data.frame(model.matrix(pgstat ~ .,stagec))[-1])

		# Impute the missing data
        dataCancerImputed <- nearestNeighborImpute(stagec_mat)
        dataCancerImputed[,1:ncol(dataCancerImputed)] <- sapply(dataCancerImputed,as.numeric)	

		# Cross validating a LDA classifier.
		# 80% Train and 20% test repeated 10 times
		cv <- randomCV(dataCancerImputed,"pgstat",MASS::lda,trainFraction = 0.8, 
		repetitions = 10,featureSelectionFunction = univariate_tstudent,
		featureSelection.control = list(limit = 0.5,thr = 0.975));

		# Compare the LDA classifier with other methods
		cp <- BinaryBenchmark(referenceCV = cv,referenceName = "LDA",
		                      referenceFilterName="t.Student")
		pl <- plot(cp,prefix = "StageC: ")

		# Default Benchmark classifiers method (BSWiMS) and filter methods. 
		# 80% Train Repeated 10 times
		cp <- BinaryBenchmark(theData = dataCancerImputed,
		theOutcome = "pgstat", reps = 10, fraction = 0.8)

		# plot the Cross Validation Metrics
		pl <- plot(cp,prefix = "Stagec:");

		# Shut down the graphics device driver
		dev.off()

		#### Regression Example ######
		# Start the graphics device driver to save all plots in a pdf format
		pdf(file = "RegressionExample.pdf",width=8, height=6)

		# Get the body fat data from the TH package

		data("bodyfat", package = "TH.data")

		# Benchmark regression methods and filter methods. 
		#80% Train Repeated 10 times
		cp <- RegresionBenchmark(theData = bodyfat, 
		theOutcome = "DEXfat", reps = 10, fraction = 0.8)

		# plot the Cross Validation Metrics
		pl <- plot(cp,prefix = "Body Fat:");
		# Shut down the graphics device driver
		dev.off()

		#### Ordinal Regression Example #####
		# Start the graphics device driver to save all plots in a pdf format
		pdf(file = "OrdinalRegressionExample.pdf",width=8, height=6)


		# Get the GBSG2 data
		data("GBSG2", package = "TH.data")

		# Prepare the model frame for benchmarking
		GBSG2$time <- NULL;
		GBSG2$cens <- NULL;
		GBSG2_mat <- cbind(tgrade = as.numeric(GBSG2$tgrade),
		as.data.frame(model.matrix(tgrade~.,GBSG2))[-1])

		# Benchmark regression methods and filter methods. 
		#30% Train Repeated 10 times
		cp <- OrdinalBenchmark(theData = GBSG2_mat, 
		theOutcome = "tgrade", reps = 10, fraction = 0.3)

		# plot the Cross Validation Metrics
		pl <- plot(cp,prefix = "GBSG:");

		# Shut down the graphics device driver
		dev.off()

	}

}
\keyword{Benchmarking}