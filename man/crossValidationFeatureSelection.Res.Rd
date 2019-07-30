\name{crossValidationFeatureSelection_Res}
\alias{crossValidationFeatureSelection_Res}
\title{NeRI-based selection of a linear, logistic, or Cox proportional hazards regression model from a set of candidate variables}
\description{
	This function performs a cross-validation analysis of a feature selection algorithm based on net residual improvement (NeRI) to return a predictive model.
	It is composed of a NeRI-based feature selection followed by an update procedure, ending with a bootstrapping backwards feature elimination.
	The user can control how many train and blind test sets will be evaluated.
}
\usage{
	crossValidationFeatureSelection_Res(size = 10,
	                                    fraction = 1.0,
	                                    pvalue = 0.05,
	                                    loops = 100,
	                                    covariates = "1",
	                                    Outcome,
	                                    timeOutcome = "Time",
	                                    variableList,
	                                    data,
	                                    maxTrainModelSize = 20,
	                                    type = c("LM", "LOGIT", "COX"),
	                                    testType = c("Binomial",
	                                                 "Wilcox",
	                                                 "tStudent",
	                                                 "Ftest"),
	                                    startOffset = 0,
	                                    elimination.bootstrap.steps = 100,
	                                    trainFraction = 0.67,
	                                    trainRepetition = 9,
	                                    setIntersect = 1,
	                                    unirank = NULL,
	                                    print=TRUE,
	                                    plots=TRUE,
	                                    lambda="lambda.1se",
	                                    equivalent=FALSE,
	                                    bswimsCycles=10,
	                                    usrFitFun=NULL,
	                                    featureSize=0)
}
\arguments{
	\item{size}{
		The number of candidate variables to be tested (the first \code{size} variables from \code{variableList})
	}
	\item{fraction}{
		The fraction of data (sampled with replacement) to be used as train
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{covariates}{
		A string of the type "1 + var1 + var2" that defines which variables will always be included in the models (as covariates)
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{startOffset}{
		Only terms whose position in the model is larger than the \code{startOffset} are candidates to be removed
	}
	\item{elimination.bootstrap.steps}{
		The number of bootstrap loops for the backwards elimination procedure
	}
	\item{trainFraction}{
		The fraction of data (sampled with replacement) to be used as train for the cross-validation procedure
	}
	\item{setIntersect}{
		The intersect of the model (To force a zero intersect, set this value to 0)
	}
	\item{trainRepetition}{
		The number of cross-validation folds (it should be at least equal to \eqn{1/trainFraction} for a complete cross-validation)
	}
	\item{unirank}{
		A list with the results yielded by the \code{uniRankVar} function, required only if the rank needs to be updated during the cross-validation procedure
	}
	\item{print}{
		Logical. If \code{TRUE}, information will be displayed
	}
	\item{plots}{
		Logical. If \code{TRUE}, plots are displayed
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
		A user fitting function to be evaluated by the cross validation procedure
	}
	\item{featureSize}{
		The original number of features to be explored in the data frame.
	}
}
\details{
	This function produces a set of data and plots that can be used to inspect the degree of over-fitting or shrinkage of a model.
	It uses bootstrapped data, cross-validation data, and, if possible, retrain data.
}
\value{
	\item{formula.list}{
		A list containing objects of class \code{formula} with the formulas used to fit the models found at each cycle
	}
	\item{Models.testPrediction}{
		A data frame with the blind test set predictions made at each fold of the cross validation (Full B:SWiMS,Median,Bagged,Forward,Backward Elimination), where the models used to generate such predictions (\code{formula.list}) were generated via a feature selection process which included only the train set.
		It also includes a column with the \code{Outcome} of each prediction, and a column with the number of the fold at which the prediction was made.
	}
	\item{FullBSWiMS.testPrediction}{
		A data frame similar to \code{Models.testPrediction}, but where the model used to generate the predictions was the Full model, generated via a feature selection process which included all data.
	}
	\item{BSWiMS}{
		A list containing the values returned by \code{bootstrapVarElimination_Res} using all data and the model from \code{updatedforwardModel}
	}
	\item{forwardSelection}{
		A list containing the values returned by \code{ForwardSelection.Model.Res} using all data
	}
	\item{updatedforwardModel}{
		A list containing the values returned by \code{updateModel.Res} using all data and the model from \code{forwardSelection}
	}
	\item{testRMSE}{
		The global blind test root-mean-square error (RMSE) of the cross-validation procedure
	}
	\item{testPearson}{
		The global blind test Pearson \emph{r} product-moment correlation coefficient of the cross-validation procedure
	}
	\item{testSpearman}{
		The global blind test Spearman \eqn{\rho} rank correlation coefficient of the cross-validation procedure
	}
	\item{FulltestRMSE}{
		The global blind test RMSE of the Full model
	}
	\item{FullTestPearson}{
		The global blind test Pearson \emph{r} product-moment correlation coefficient of the Full model
	}
	\item{FullTestSpearman}{
		The global blind test Spearman \eqn{\rho} rank correlation coefficient of the Full model
	}
	\item{trainRMSE}{
		The train RMSE at each fold of the cross-validation procedure
	}
	\item{trainPearson}{
		The train Pearson \emph{r} product-moment correlation coefficient at each fold of the cross-validation procedure
	}
	\item{trainSpearman}{
		The train Spearman \eqn{\rho} rank correlation coefficient at each fold of the cross-validation procedure
	}
	\item{FullTrainRMSE}{
		The train RMSE of the Full model at each fold of the cross-validation procedure
	}
	\item{FullTrainPearson}{
		The train Pearson \emph{r} product-moment correlation coefficient of the Full model at each fold of the cross-validation procedure
	}
	\item{FullTrainSpearman}{
		The train Spearman \eqn{\rho} rank correlation coefficient of the Full model at each fold of the cross-validation procedure
	}
	\item{testRMSEAtFold}{
		The blind test RMSE at each fold of the cross-validation procedure
	}
	\item{FullTestRMSEAtFold}{
		The blind test RMSE of the Full model at each fold of the cross-validation procedure
	}
	\item{Fullenet}{
		An object of class \code{cv.glmnet} containing the results of an elastic net cross-validation fit
	}
	\item{LASSO.testPredictions}{
		A data frame similar to \code{Models.testPrediction}, but where the predictions were made by the elastic net model
	}
	\item{LASSOVariables}{
		A list with the elastic net Full model and the models found at each cross-validation fold
	}
	\item{byFoldTestMS}{
		A vector with the Mean Square error for each blind fold
	}
	\item{byFoldTestSpearman}{
		A vector with the Spearman correlation between prediction and outcome for each blind fold
	}
	\item{byFoldTestPearson}{
		A vector with the Pearson correlation between prediction and outcome for each blind fold
	}
	\item{byFoldCstat}{
		A vector with the C-index (Somers' Dxy rank correlation :\code{rcorr.cens}) between prediction and outcome for each blind fold
	}
	\item{CVBlindPearson}{
		A vector with the Pearson correlation between the outcome and prediction for each repeated experiment
	}
	\item{CVBlindSpearman}{
		A vector with the Spearm correlation between the outcome and prediction for each repeated experiment
	}
	\item{CVBlindRMS}{
		A vector with the RMS between the outcome and prediction for each repeated experiment
	}
	\item{Models.trainPrediction}{
		A data frame with the outcome and the train prediction of every model
	}
	\item{FullBSWiMS.trainPrediction}{
		A data frame with the outcome and the train prediction at each CV fold for the main model
	}
	\item{LASSO.trainPredictions}{
		A data frame with the outcome and the prediction of each enet lasso model
	}	
	\item{uniTrainMSS}{
		A data frame with mean square of the train residuals from the univariate models of the model terms
	}	
	\item{uniTestMSS}{
		A data frame with mean square of the test residuals of the univariate models of the model terms
	}
	\item{BSWiMS.ensemble.prediction}{
		The ensemble prediction by all models on the test data
	}
	\item{AtOptFormulas.list}{
		The list of formulas with "optimal" performance
	}
	\item{ForwardFormulas.list}{
		The list of formulas produced by the forward procedure
	}
	\item{baggFormulas.list}{
		The list of the bagged models
	}
	\item{LassoFilterVarList}{
		The list of variables used by LASSO fitting
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{crossValidationFeatureSelection_Bin},
				 \link{improvedResiduals},
				 \link{bootstrapVarElimination_Res}}}

\keyword{Model_Generation}