\name{crossValidationFeatureSelection_Bin}
\alias{crossValidationFeatureSelection_Bin}
\title{IDI/NRI-based selection of a linear, logistic, or Cox proportional hazards regression model from a set of candidate variables}
\description{
	This function performs a cross-validation analysis of a feature selection algorithm based on the integrated discrimination improvement (IDI) or the net reclassification improvement (NRI) to return a predictive model.
	It is composed of an IDI/NRI-based feature selection followed by an update procedure, ending with a bootstrapping backwards feature elimination.
	The user can control how many train and blind test sets will be evaluated.
}
\usage{
	crossValidationFeatureSelection_Bin(size = 10,
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
	                                selectionType = c("zIDI", "zNRI"),
	                                startOffset = 0,
	                                elimination.bootstrap.steps = 100,
	                                trainFraction = 0.67,
	                                trainRepetition = 9,
	                                bootstrap.steps = 100,
	                                nk = 0,
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
		The maximum \emph{p}-value, associated to either IDI or NRI, allowed for a term in the model
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
	\item{selectionType}{
		The type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-score of IDI or of NRI
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
	\item{trainRepetition}{
		The number of cross-validation folds (it should be at least equal to \eqn{1/}\code{trainFraction} for a complete cross-validation)
	}
	\item{bootstrap.steps}{
		The number of bootstrap loops for the confidence intervals estimation
	}
	\item{nk}{
		The number of neighbours used to generate a \emph{k}-nearest neighbours (KNN) classification. If zero, \emph{k} is set to the square root of the number of cases. If less than zero, it will not perform the KNN classification
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
	During each cycle, a train and a test ROC will be generated using bootstrapped data.
	At the end of the cross-validation feature selection procedure, a set of three plots may be produced depending on the specifications of the analysis.
	The first plot shows the ROC for each cross-validation blind test.
	The second plot, if enough samples are given, shows the ROC of each model trained and tested in the blind test partition.
	The final plot shows ROC curves generated with the train, the bootstrapped blind test, and the cross-validation test data.
	Additionally, this plot will also contain the ROC of the cross-validation mean test data, and of the cross-validation coherence.
	These set of plots may be used to get an overall perspective of the expected model shrinkage.
	Along with the plots, the function provides the overall performance of the system (accuracy, sensitivity, and specificity).
	The function also produces a report of the expected performance of a KNN algorithm trained with the selected features of the model, and an elastic net algorithm.
	The test predictions obtained with these algorithms can then be compared to the predictions generated by the logistic, linear, or Cox proportional hazards regression model.
}
\value{
	\item{formula.list}{
		A list containing objects of class \code{formula} with the formulas used to fit the models found at each cycle
	}
	\item{Models.testPrediction}{
		A data frame with the blind test set predictions (Full B:SWiMS,Median,Bagged,Forward,Backwards Eliminations) made at each fold of the cross validation, where the models used to generate such predictions (\code{formula.list}) were generated via a feature selection process which included only the train set.
		It also includes a column with the \code{Outcome} of each prediction, and a column with the number of the fold at which the prediction was made.
	}
	\item{FullBSWiMS.testPrediction}{
		A data frame similar to \code{Models.testPrediction}, but where the model used to generate the predictions was the Full model, generated via a feature selection process which included all data.
	}
	\item{TestRetrained.blindPredictions}{
		A data frame similar to \code{Models.testPrediction}, but where the models were retrained on an independent set of data (only if enough samples are given at each fold)
	}
	\item{LastTrainBSWiMS.bootstrapped}{
		An object of class \code{bootstrapValidation_Bin} containing the results of the bootstrap validation in the last trained model
	}
	\item{Test.accuracy}{
		The global blind test accuracy of the cross-validation procedure
	}
	\item{Test.sensitivity}{
		The global blind test sensitivity of the cross-validation procedure
	}
	\item{Test.specificity}{
		The global blind test specificity of the cross-validation procedure
	}
	\item{Train.correlationsToFull}{
		The Spearman \eqn{\rho} rank correlation coefficient between the predictions made with each model from \code{formula.list} and the Full model in the train set
	}
	\item{Blind.correlationsToFull}{
		The Spearman \eqn{\rho} rank correlation coefficient between the predictions made with each model from \code{formula.list} and the Full model in the test set
	}
	\item{FullModelAtFoldAccuracies}{
		The blind test accuracy for the Full model at each cross-validation fold
	}
	\item{FullModelAtFoldSpecificties}{
		The blind test specificity for the Full model at each cross-validation fold
	}
	\item{FullModelAtFoldSensitivities}{
		The blind test sensitivity for the Full model at each cross-validation fold
	}
	\item{FullModelAtFoldAUC}{
		The blind test ROC AUC for the Full model at each cross-validation fold
	}
	\item{AtCVFoldModelBlindAccuracies}{
		The blind test accuracy for the Full model at each final cross-validation fold
	}
	\item{AtCVFoldModelBlindSpecificities}{
		The blind test specificity for the Full model at each final cross-validation fold
	}
	\item{AtCVFoldModelBlindSensitivities}{
		The blind test sensitivity for the Full model at each final cross-validation fold
	}
	\item{CVTrain.Accuracies}{
		The train accuracies at each fold
	}
	\item{CVTrain.Sensitivity}{
		The train sensitivity at each fold
	}
	\item{CVTrain.Specificity}{
		The train specificity at each fold
	}
	\item{CVTrain.AUCs}{
		The train ROC AUC for each fold
	}
	\item{forwardSelection}{
		A list containing the values returned by \code{ForwardSelection.Model.Bin} using all data
	}
	\item{updateforwardSelection}{
		A list containing the values returned by \code{updateModel.Bin} using all data and the model from \code{forwardSelection}
	}
	\item{BSWiMS}{
		A list containing the values returned by \code{bootstrapVarElimination_Bin} using all data and the model from \code{updateforwardSelection}
	}
	\item{FullBSWiMS.bootstrapped}{
		An object of class \code{bootstrapValidation_Bin} containing the results of the bootstrap validation in the Full model
	}
	\item{Models.testSensitivities }{
		A matrix with the mean ROC sensitivities at certain specificities for each train and all test cross-validation folds using the cross-validation models (i.e. 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, and 0.05)
	}
	\item{FullKNN.testPrediction}{
		A data frame similar to \code{Models.testPrediction}, but where a KNN classifier with the same features as the Full model was used to generate the predictions
	}
	\item{KNN.testPrediction}{
		A data frame similar to \code{Models.testPrediction}, but where KNN classifiers with the same features as the cross-validation models were used to generate the predictions at each cross-validation fold
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
	\item{uniTrain.Accuracies}{
		The list of accuracies of an univariate analysis on each one of the model variables in the train sets
	}
	\item{uniTest.Accuracies}{
		The list of accuracies of an univariate analysis on each one of the model variables in the test sets
	}
	\item{uniTest.TopCoherence}{
		The accuracy coherence of the top ranked variable on the test set
	}
	\item{uniTrain.TopCoherence}{
		The accuracy coherence of the top ranked variable on the train set
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
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{crossValidationFeatureSelection_Res},
				 \link{ForwardSelection.Model.Bin},
				 \link{ForwardSelection.Model.Res}}}
\keyword{Model_Generation}