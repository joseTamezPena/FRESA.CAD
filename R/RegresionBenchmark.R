RegresionBenchmark <-  function(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,referenceCV = NULL,referenceName = "Reference",referenceFilterName="Reference")
{
  if (!requireNamespace("e1071", quietly = TRUE)) {
	  install.packages("e1071", dependencies = TRUE)
	  }
  if (!requireNamespace("randomForest", quietly = TRUE)) {
	  install.packages("randomForest", dependencies = TRUE)
	  }
  if (!requireNamespace("rpart", quietly = TRUE)) {
	  install.packages("rpart", dependencies = TRUE)
	  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
	  install.packages("MASS", dependencies = TRUE)
	  }

	if (is.null(theData))
	{
		if (exists("theDataSet", envir=FRESAcacheEnv))
		{
			theData <- get("theDataSet", envir=FRESAcacheEnv);
			theOutcome <- get("theDataOutcome", envir=FRESAcacheEnv);
		}	
	}
	else
	{
		assign("theDataSet",theData,FRESAcacheEnv);
		assign("theDataOutcome",theOutcome,FRESAcacheEnv);
	}
	  
	  
	  RMSETable <- NULL 
	  CorTable <- NULL
	  BiasTable <- NULL
	  MAETable <- NULL
	  CorSpearman <- NULL
	  
	  RMSETable_filter <- NULL 
	  CorTable_filter <- NULL
	  BiasTable_filter <- NULL
	  MAETable_filter <- NULL
	  CorSpearman_filter <- NULL
	  fmeth_0 <- NULL; 
#	  par(mfrow = c(1,1));

	  FilterMethod <-  function(regresionfun = e1071::svm, regnamefunc = "",...)
	  {
		RMSETable_f <- NULL 
		CorTable_f <- NULL
		MAETable_f <- NULL
		CorSpearman_f<- NULL
		BiasTable_f <- NULL
		
		parm <- list(...);
		rcvFilter_Reference <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,...);
		
		stats <- predictionStats_regression(rcvFilter_Reference$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
				
		rcvFilter_LASSO <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvLASSO$selectedFeaturesSet,...);
		stats <- predictionStats_regression(rcvFilter_LASSO$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
		
		rcvFilter_RPART <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvRPART$selectedFeaturesSet,...);
		stats <- predictionStats_regression(rcvFilter_RPART$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);

		selectedFeaturesSet <- rcvRF$selectedFeaturesSet
		for (i in 1:length(selectedFeaturesSet))
		{
		  if (length(referenceCV$selectedFeaturesSet[[i]]) > 1)
		  {
			if (length(selectedFeaturesSet[[i]]) > length(referenceCV$selectedFeaturesSet[[i]]))
			{
			  selectedFeaturesSet[[i]] <- selectedFeaturesSet[[i]][1:length(referenceCV$selectedFeaturesSet[[i]])];
			}
		  }
		  else # the top five or RF
		  {
			warning ("Less than 2 features, then will keep the top five of RF\n")
			if (length(selectedFeaturesSet[[i]]) > 5)
			{
			  selectedFeaturesSet[[i]] <- selectedFeaturesSet[[i]][1:5];
			}
		  }
		}
		rcvFilter_RF <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = selectedFeaturesSet,...);
		stats <- predictionStats_regression(rcvFilter_RF$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);

		
		if (is.null(fmeth_0))
		{
			rcvFilter_FT <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_residual,featureSelection.control = list(uniTest = "Ftest",limit = 0.35,thr = 0.975),...);
		}
		else
		{
			rcvFilter_FT <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_FT$selectedFeaturesSet,...);
		}
		stats <- predictionStats_regression(rcvFilter_FT$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
				
		if (is.null(fmeth_0))
		{
				rcvFilter_Wt <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_residual,featureSelection.control = list(uniTest = "Wilcox",limit = 0.35,thr = 0.975),...);
		}
		else
		{
			rcvFilter_Wt <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_Wt$selectedFeaturesSet,...);
		}
		stats <- predictionStats_regression(rcvFilter_Wt$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
		
		
		if (is.null(fmeth_0))
		{
				rcvFilter_pearson <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_correlation,featureSelection.control = list(method = "pearson",limit = 0.35,thr = 0.975),...);
		}
		else
		{
			rcvFilter_pearson <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_pearson$selectedFeaturesSet,...);
		}
		stats <- predictionStats_regression(rcvFilter_pearson$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
		
		if (is.null(fmeth_0))
		{
			rcvFilter_kendall <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_correlation,featureSelection.control = list(method = "kendall",limit = 0.35,thr = 0.975),...);
		}
		else
		{
			rcvFilter_kendall <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_kendall$selectedFeaturesSet,...);
		}
		stats <- predictionStats_regression(rcvFilter_kendall$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
		
		if (is.null(fmeth_0))
		{
			rcvFilter_mRMR <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,...);
		}
		else
		{
			rcvFilter_mRMR <- randomCV(theData,theOutcome,regresionfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_mRMR$selectedFeaturesSet,...);
		}
		stats <- predictionStats_regression(rcvFilter_mRMR$medianTest);
		CorTable_f <- rbind(CorTable_f,stats$corci);
		BiasTable_f <- rbind(BiasTable_f,stats$biasci);
		RMSETable_f <- rbind(RMSETable_f,stats$RMSEci);
		CorSpearman_f <- rbind(CorSpearman_f,stats$spearmanci);
		MAETable_f <- rbind(MAETable_f,stats$MAEci);
		
		result <- list(RMSETable_filter = RMSETable_f,
					   CorTable_filter = CorTable_f,
					   BiasTable_filter = BiasTable_f,
					   MAE_filter = MAETable_f,
					   Spearam_filter = CorSpearman_f,
					   rcvFilter_Reference = rcvFilter_Reference,
					   rcvFilter_LASSO = rcvFilter_LASSO,
					   rcvFilter_RPART = rcvFilter_RPART,
					   rcvFilter_RF = rcvFilter_RF,
					   rcvFilter_FT = rcvFilter_FT,
					   rcvFilter_Wt = rcvFilter_Wt,
					   rcvFilter_pearson = rcvFilter_pearson,
					   rcvFilter_kendall = rcvFilter_kendall,
					   rcvFilter_mRMR = rcvFilter_mRMR
		)

		return(result);
	  }
	  
	  
	  
	  ######################Regression Algorithms####################################  
	  
	  if (is.null(referenceCV))
	  {
		referenceCV <- randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
		referenceFilterName = "BSWiMS";
		referenceName = "BSWiMS";
	  }
	  else
	  {
			reps <- referenceCV$repetitions;
	  }

	  stats <- predictionStats_regression(referenceCV$medianTest,referenceName);
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);
	  	  
	  rcvRF <- randomCV(theData,theOutcome,randomForest::randomForest,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self");
	  
	  stats <- predictionStats_regression(rcvRF$medianTest,"RF");
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);
	  
	  rcvRPART <- randomCV(theData,theOutcome,rpart::rpart,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self");
	  
	  stats <- predictionStats_regression(rcvRPART$medianTest,"RPART");
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);

	  rcvLASSO <- randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self");
	  
	  stats <- predictionStats_regression(rcvLASSO$medianTest,"LASSO");
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);
	  
	  rcvSVM <- randomCV(theData,theOutcome,e1071::svm,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA);
	  
	  stats <- predictionStats_regression(rcvSVM$medianTest,"SVM");
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);
	  
	#  Method Meta Ensemble
	  ens <- rowMeans(cbind(referenceCV$medianTest[,2],rcvRF$medianTest[,2],rcvLASSO$medianTest[,2],rcvSVM$medianTest[,2]))

	  stats <- predictionStats_regression(cbind(rcvSVM$medianTest[,1],ens),"Ensemble");
	  CorTable <- rbind(CorTable,stats$corci);
	  BiasTable <- rbind(BiasTable,stats$biasci);
	  RMSETable <- rbind(RMSETable,stats$RMSEci);
	  CorSpearman <- rbind(CorSpearman,stats$spearmanci);
	  MAETable <- rbind(MAETable,stats$MAEci);
	  
	  
	  ######################Filters ####################################  
	  
	  cat("Robust Regression\n")
	  fmeth <- FilterMethod(MASS::rlm,"Robust Regression")
	  fmeth_0 <- fmeth;
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);
	  
	  cat("LASSO\n")
	  fmeth <- FilterMethod(LASSO_MIN,"LASSO")
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);
	  
	  cat("Filtered SVM\n")
	  fmeth <- FilterMethod(e1071::svm,"SVM")
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);
	  
	  cat("Random Forest\n")
	  fmeth <- FilterMethod(randomForest::randomForest,"Random Forest")
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);
  
	  cat("Linear Regression\n")
	  fmeth <- FilterMethod(lm,"Linear Regression")
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);

	  cat("Ridge Regression\n")
	  fmeth <- FilterMethod(LM_RIDGE_MIN,"Ridge Regression")
	  CorTable_filter <- rbind(CorTable_filter,fmeth$CorTable_filter);
	  BiasTable_filter <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
	  RMSETable_filter <- rbind(RMSETable_filter,fmeth$RMSETable_filter);
	  CorSpearman_filter <- rbind(CorSpearman_filter,fmeth$Spearam_filter);
	  MAETable_filter <- rbind(MAETable_filter,fmeth$MAE_filter);

	
	test_Predictions <- referenceCV$medianTest;
	tnames <- rownames(test_Predictions);
	testres <- rep(NA,nrow(test_Predictions));
	names(testres) <- tnames;
	testresnamed <- testres;
	testres[tnames %in% rownames(rcvRF$medianTest)] <- rcvRF$medianTest[tnames[tnames %in% rownames(rcvRF$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(rcvLASSO$medianTest)] <- rcvLASSO$medianTest[tnames[tnames %in% rownames(rcvLASSO$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(rcvRPART$medianTest)] <- rcvRPART$medianTest[tnames[tnames %in% rownames(rcvRPART$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(rcvSVM$medianTest)] <- rcvSVM$medianTest[tnames[tnames %in% rownames(rcvSVM$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
 	test_Predictions <- cbind(test_Predictions,ens);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_Reference$medianTest)] <- fmeth$rcvFilter_Reference$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_Reference$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_LASSO$medianTest)] <- fmeth$rcvFilter_LASSO$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_LASSO$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_RPART$medianTest)] <- fmeth$rcvFilter_RPART$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_RPART$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_RF$medianTest)] <- fmeth$rcvFilter_RF$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_RF$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_FT$medianTest)] <- fmeth$rcvFilter_FT$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_FT$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_Wt$medianTest)] <- fmeth$rcvFilter_Wt$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_Wt$medianTest)],2];
	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_pearson$medianTest)] <- fmeth$rcvFilter_pearson$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_pearson$medianTest)],2];
 	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_kendall$medianTest)] <- fmeth$rcvFilter_kendall$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_kendall$medianTest)],2];
 	test_Predictions <- cbind(test_Predictions,testres);
	testres <- testresnamed;
	testres[tnames %in% rownames(fmeth$rcvFilter_mRMR$medianTest)] <- fmeth$rcvFilter_mRMR$medianTest[tnames[tnames %in% rownames(fmeth$rcvFilter_mRMR$medianTest)],2];
 	test_Predictions <- cbind(test_Predictions,testres)
	
	
	colnames(test_Predictions) <- c("Outcome",referenceName,"RF","LASSO","RPART","SVM.mRMR","Ensemble",paste("RIDGE.",referenceFilterName,sep=""),"RIDGE.LASSO","RIDGE.RPART","RIDGE.RF.ref","RIDGE.FT","RIDGE.Wt","RIDGE.Pearson","RIDGE.Kendall","RIDGE.mRMR");
	test_Predictions <- as.data.frame(test_Predictions)


	thesets <- c("Regression Algorithm")
	theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","ENS")

	theRegressMethod <- c("Robust Regression","LASSO","SVM","Random Forest","Linear Regression","Ridge Regression")
	theFiltersets <- c(referenceFilterName,"LASSO","RPART","RF.ref","F-Test","W-Test","Pearson","Kendall","mRMR")

	rownames(CorTable) <- theMethod;
	rownames(BiasTable) <- theMethod;
	rownames(RMSETable) <- theMethod;
	rownames(CorSpearman) <- theMethod
	rownames(MAETable) <- theMethod

	
    ff <- names(referenceCV$featureFrequency)
	ff <- c(ff,names(rcvLASSO$featureFrequency))
	ff <- c(ff,names(rcvRPART$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_RF$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_FT$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_pearson$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_kendall$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_mRMR$featureFrequency))
	ff <- unique(ff)

	Nvar <- length(ff);
	selFrequency <- matrix(0,nrow = Nvar,ncol = length(theFiltersets))
#	rownames(selFrequency) <- ff[1:Nvar]
	rownames(selFrequency) <- ff
	selnames <- rownames(selFrequency)
	colnames(selFrequency) <- theFiltersets
	ff <- referenceCV$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,referenceFilterName] <- ff[selnames[fnames]]
	ff <- rcvLASSO$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"LASSO"] <- ff[selnames[fnames]]
	ff <- rcvRPART$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"RPART"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_RF$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"RF.ref"] <- ff[selnames[fnames]]

	ff <- fmeth$rcvFilter_FT$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"F-Test"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_Wt$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"W-Test"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_pearson$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"Pearson"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_kendall$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"Kendall"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_mRMR$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"mRMR"] <- ff[selnames[fnames]]
	selFrequency <- selFrequency/reps


	theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","ENS")

	elapcol <- names(referenceCV$theTimes) == "elapsed"
	cputimes <- list(Reference = mean(referenceCV$theTimes[ elapcol ]),RF = mean(rcvRF$theTimes[ elapcol ]),RPART = mean(rcvRPART$theTimes[ elapcol ]),LASSO = mean(rcvLASSO$theTimes[ elapcol ]),SVM = mean(rcvSVM$theTimes[ elapcol ]));


	jaccard_filter = list(Reference = referenceCV$jaccard,
                                       LASSO = rcvLASSO$jaccard,
                                       rpart = rcvRPART$jaccard,
                                       RF.ref = fmeth$rcvFilter_RF$jaccard,
                                       FT = fmeth$rcvFilter_FT$jaccard,
                                       WT = fmeth$rcvFilter_Wt$jaccard,
                                       pearson = fmeth$rcvFilter_pearson$jaccard,
                                       kendall = fmeth$rcvFilter_kendall$jaccard,
                                       mRMR = fmeth$rcvFilter_mRMR$jaccard
                 );
	featsize <- unlist(lapply(jaccard_filter, `[`, c('averageLength')))
	names(featsize) <- theFiltersets;
	jaccard <- unlist(lapply(jaccard_filter, `[`, c('Jaccard.SM')))
	names(jaccard) <- theFiltersets;
	cputimes <- unlist(cputimes);
	cputimes <- c(cputimes,sum(cputimes));
	names(cputimes) <- theMethod;


    result <- list(CorTable = CorTable,BiasTable = BiasTable,RMSETable = RMSETable,
                 CorTable_filter = CorTable_filter,BiasTable_filter = BiasTable_filter,RMSETable_filter = RMSETable_filter,
				 CorSpearman = CorSpearman, MAETable = MAETable,
				 CorSpearman_filter = CorSpearman_filter, MAETable_filter = MAETable_filter,
                 times = list(Reference = referenceCV$theTimes,RF = rcvRF$theTimes,rpart = rcvRPART$theTimes,LASSO = rcvLASSO$theTimes,SVM = rcvSVM$theTimes), 
                 featsize = featsize,
				 jaccard = jaccard,
                 TheCVEvaluations = list(Reference = referenceCV,
                                         RF = rcvRF,
                                         LASSO = rcvLASSO,
                                         RPART = rcvRPART,
                                         SVM = rcvSVM,
                                         RF.ref = fmeth$rcvFilter_RF,
                                         FT = fmeth$rcvFilter_FT,
                                         WT = fmeth$rcvFilter_Wt,
                                         pearson = fmeth$rcvFilter_pearson,
                                         kendall = fmeth$rcvFilter_kendall,
                                         mRMR = fmeth$rcvFilter_mRMR
                 ),
				 thesets = thesets,
				 theMethod = theMethod,
				 theRegressMethod = theRegressMethod,
				 theFiltersets = theFiltersets,
				 testPredictions = test_Predictions,
				 featureSelectionFrequency = selFrequency,
				 cpuElapsedTimes=cputimes				 
  )
  
  class(result) <- c("FRESA_benchmark","Regression");

  
  return(result)
}