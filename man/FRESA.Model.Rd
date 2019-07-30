\name{FRESA.Model}
\alias{FRESA.Model}
\title{Automated model selection}
\description{
	This function uses a wrapper procedure to select the best features of a non-penalized linear model that best predict the outcome, given the formula of an initial model template (linear, logistic, or Cox proportional hazards), an optimization procedure, and a data frame.
	A filter scheme may be enabled to reduce the search space of the wrapper procedure.
	The false selection rate may be empirically controlled by enabling bootstrapping, and model shrinkage can be evaluated by cross-validation. 
}
\usage{
	FRESA.Model(formula,
	            data,
	            OptType = c("Binary", "Residual"),
	            pvalue = 0.05,
	            filter.p.value = 0.10,
	            loops = 32,
	            maxTrainModelSize = 20,
	            elimination.bootstrap.steps = 100,
	            bootstrap.steps = 100,
	            print = FALSE,
	            plots = FALSE,
	            CVfolds = 1,
	            repeats = 1,
	            nk = 0,
	            categorizationType = c("Raw",
	                                   "Categorical",
	                                   "ZCategorical",
	                                   "RawZCategorical",
	                                   "RawTail",
	                                   "RawZTail",
	                                   "Tail",
	                                   "RawRaw"),
	            cateGroups = c(0.1, 0.9),
	            raw.dataFrame = NULL,
	            var.description = NULL,
	            testType = c("zIDI",
	                         "zNRI",
	                         "Binomial",
	                         "Wilcox",
	                         "tStudent",
	                         "Ftest"),
	            lambda="lambda.1se",
	            equivalent=FALSE,
	            bswimsCycles=20,
	            usrFitFun=NULL
	            )

}
\arguments{
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{OptType}{
		Optimization type: Based on the integrated discrimination improvement (Binary) index for binary classification ("Binary"), or based on the net residual improvement (NeRI) index for linear regression ("Residual") 
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the \code{testType}, allowed for a term in the model (it will control the false selection rate)
	}
	\item{filter.p.value}{
		The maximum \emph{p}-value, for a variable to be included to the feature selection procedure
	}
	\item{loops}{
		The number of bootstrap loops for the forward selection procedure
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
	\item{elimination.bootstrap.steps}{
		The number of bootstrap loops for the backwards elimination procedure
	}
	\item{bootstrap.steps}{
		The number of bootstrap loops for the bootstrap validation procedure
	}
	\item{print}{
		Logical. If \code{TRUE}, information will be displayed
	}
	\item{plots}{
		Logical. If \code{TRUE}, plots are displayed
	}
	\item{CVfolds}{
		The number of folds for the final cross-validation
	}
	\item{repeats}{
		The number of times that the cross-validation procedure will be repeated
	}
	\item{nk}{
		The number of neighbors used to generate a \emph{k}-nearest neighbors (KNN) classification. If zero, \emph{k} is set to the square root of the number of cases. If less than zero, it will not perform the KNN classification
	}
	\item{categorizationType}{
		How variables will be analyzed: As given in \code{data} ("Raw"); broken into the \emph{p}-value categories given by \code{cateGroups} ("Categorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, and weighted by the \emph{z}-score ("ZCategorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, weighted by the \emph{z}-score, plus the raw values ("RawZCategorical"); raw values, plus the tails ("RawTail"); or raw values, weighted by the \emph{z}-score, plus the tails ("RawZTail")
	}
	\item{cateGroups}{
		A vector of percentiles to be used for the categorization procedure
	}
	\item{raw.dataFrame}{
		A data frame similar to \code{data}, but with unadjusted data, used to get the means and variances of the unadjusted data
	}
	\item{var.description}{
		A vector of the same length as the number of columns of \emph{data}, containing a description of the variables
	}
	\item{testType}{
		For an Binary-based optimization, the type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-value of Binary or of NRI. For a NeRI-based optimization, the type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{lambda}{
		The passed value to the s parameter of the glmnet cross validation coefficient
	}
	\item{equivalent}{
		Is set to TRUE CV will compute the equivalent model
	}
	\item{bswimsCycles}{
		The maximum number of models to be returned by \code{BSWiMS.model}
	}
	\item{usrFitFun}{
		An optional user provided fitting function to be evaluated by the cross validation procedure: fitting: usrFitFun(formula,data), with a predict function
	}
}
\details{
This important function of FRESA.CAD will model or cross validate the models. Given an outcome formula, and a data.frame this function will do an univariate analysis of the data (\code{univariateRankVariables}),
then it will select the top ranked variables; after that it will select the model that best describes the outcome. At output it will return the bootstrapped performance of the model 
(\code{bootstrapValidation_Bin} or \code{bootstrapValidation_Res}). It can be set to report the cross-validation performance of the selection process which will return either 
a  \code{crossValidationFeatureSelection_Bin} or a \code{crossValidationFeatureSelection_Res} object.
}
\value{
	\item{BSWiMS.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the final model
	}
	\item{reducedModel}{
		The resulting object of the backward elimination procedure
	}
	\item{univariateAnalysis}{
		A data frame with the results from the univariate analysis
	}
	\item{forwardModel}{
		The resulting object of the feature selection function.
	}
	\item{updatedforwardModel}{
		The resulting object of the the update procedure
	}
	\item{bootstrappedModel}{
		The resulting object of the bootstrap procedure on \code{final.model}
	}
	\item{cvObject}{
		The resulting object of the cross-validation procedure
	}
	\item{used.variables}{
		The number of terms that passed the filter procedure
	}
	\item{call}{
		the function call
	}
	
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\examples{
	\dontrun{

		# Start the graphics device driver to save all plots in a pdf format
		pdf(file = "FRESA.Model.Example.pdf",width = 8, height = 6)
		# Get the stage C prostate cancer data from the rpart package
		data(stagec,package = "rpart")
		options(na.action = 'na.pass')
		stagec_mat <- cbind(pgstat = stagec$pgstat,
		    pgtime = stagec$pgtime,
		    as.data.frame(model.matrix(Surv(pgtime,pgstat) ~ .,stagec))[-1])

		data(cancerVarNames)
		dataCancerImputed <- nearestNeighborImpute(stagec_mat)

		# Get a Cox proportional hazards model using:
		# - The default parameters
		md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
						  data = dataCancerImputed,
						  var.description = cancerVarNames[,2])
		pt <- plot(md$bootstrappedModel)
		sm <- summary(md$BSWiMS.model)
		print(sm$coefficients)


		# Get a 10 fold CV Cox proportional hazards model using:
		# - Repeat 10 times de CV
		md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
						  data = dataCancerImputed, CVfolds = 10, 
						  repeats = 10,
						  var.description = cancerVarNames[,2])
		pt <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds = 10)
		print(pt$predictionTable)

		pt <- plotModels.ROC(md$cvObject$LASSO.testPredictions,theCVfolds = 10)
		pt <- plotModels.ROC(md$cvObject$KNN.testPrediction,theCVfolds = 10)

		# Get a  regression of the survival time

		timeSubjects <- dataCancerImputed
		timeSubjects$pgtime <- log(timeSubjects$pgtime)

		md <- FRESA.Model(formula = pgtime ~ 1,
						  data = timeSubjects,
						  var.description = cancerVarNames[,2])
		pt <- plot(md$bootstrappedModel)
		sm <- summary(md$BSWiMS.model)
		print(sm$coefficients)

		# Get a logistic regression model using
		# - The default parameters and removing time as possible predictor

		dataCancerImputed$pgtime <- NULL

		md <- FRESA.Model(formula = pgstat ~ 1,
						  data = dataCancerImputed,
						  var.description = cancerVarNames[,2])
		pt <- plot(md$bootstrappedModel)
		sm <- summary(md$BSWiMS.model)
		print(sm$coefficients)

		# Get a logistic regression model using:
		# - residual-based optimization
		md <- FRESA.Model(formula = pgstat ~ 1,
						  data = dataCancerImputed,
						  OptType = "Residual",
						  var.description = cancerVarNames[,2])
		pt <- plot(md$bootstrappedModel)
		sm <- summary(md$BSWiMS.model)
		print(sm$coefficients)


		# Shut down the graphics device driver
		dev.off()

	}
}
\keyword{Model_Generation}