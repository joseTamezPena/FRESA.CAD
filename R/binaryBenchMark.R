BinaryBenchmark <-  function(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,referenceCV = NULL,referenceName = "Reference",referenceFilterName="Reference")
{
  if (!requireNamespace("epiR", quietly = TRUE)) {
	  install.packages("epiR", dependencies = TRUE)
	  }
  if (!requireNamespace("e1071", quietly = TRUE)) {
	  install.packages("e1071", dependencies = TRUE)
	  }
  if (!requireNamespace("randomForest", quietly = TRUE)) {
	  install.packages("randomForest", dependencies = TRUE)
	  }
  if (!requireNamespace("rpart", quietly = TRUE)) {
	  install.packages("rpart", dependencies = TRUE)
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
	  
  aucTable <- NULL;
  accciTable <- NULL;
  errorciTable <- NULL;
  senTable <- NULL;
  speTable <- NULL;
  cidxTable <- NULL;

  aucTable_filter <- NULL;
  accciTable_filter <- NULL;
  errorciTable_filter <- NULL;
  senciTable_filter <- NULL;
  speciTable_filter <- NULL;
  cindexTable_filter <- NULL;
  fmeth_0 <- NULL;
  
#   par(mfrow = c(1,1));
  FilterMethod <-  function(clasfun = e1071::svm, classname = "", center = FALSE, ...)
  {
    errorciTable_filter <- NULL;
    aucTable_filter <- NULL;
    accciTable_filter <- NULL;
    senTable_filter <- NULL;
    speTable_filter <- NULL;
	cindexTable_filter <- NULL;

    rcvFilter_reference <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,...);
    cStats <- predictionStats_binary(rcvFilter_reference$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
    rcvFilter_LASSO <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvLASSO$selectedFeaturesSet,...);
    cStats <- predictionStats_binary(rcvFilter_LASSO$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
    rcvFilter_RPART <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvRPART$selectedFeaturesSet,...);
    cStats <- predictionStats_binary(rcvFilter_RPART$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    

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
        warning ("Reference features less than 2, then will keep the top five of RF\n")
        if (length(selectedFeaturesSet[[i]]) > 5)
        {
          selectedFeaturesSet[[i]] <- selectedFeaturesSet[[i]][1:5];
        }
		if (length(selectedFeaturesSet[[i]]) < 2)
		{
			selectedFeaturesSet <- "Self";
		}
      }
    }
    rcvFilter_RF <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = selectedFeaturesSet,...);
    cStats <- predictionStats_binary(rcvFilter_RF$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
	if (is.null(fmeth_0))
	{
		cat("zIDI Feature Selection: ");
	    rcvFilter_IDI <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_Logit,featureSelection.control = list(uniTest = "zIDI",limit = 0.9,thr = 0.975),...);
	}
	else
	{
	    rcvFilter_IDI <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_IDI$selectedFeaturesSet,...);
	}
    cStats <- predictionStats_binary(rcvFilter_IDI$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
	if (is.null(fmeth_0))
	{
		cat("zNRI Feature Selection: ");
	    rcvFilter_NRI <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_Logit,featureSelection.control = list(uniTest = "zNRI",limit = 0.9,thr = 0.975),...);
	}
	else
	{
	    rcvFilter_NRI <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_NRI$selectedFeaturesSet,...);
	}
    cStats <- predictionStats_binary(rcvFilter_NRI$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
	if (is.null(fmeth_0))
	{
		cat("t student Feature Selection: ");
	    rcvFilter_tStudent <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_tstudent,featureSelection.control = list(limit = 0.9,thr = 0.975),...);
	}
	else
	{
	    rcvFilter_tStudent <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_tStudent$selectedFeaturesSet,...);
	}
    cStats <- predictionStats_binary(rcvFilter_tStudent$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
	if (is.null(fmeth_0))
	{
		cat("Wilcoxon Feature Selection: ");
	    rcvFilter_wilcox <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_Wilcoxon,featureSelection.control = list(limit = 0.9,thr = 0.975),...);
 	}
	else
	{
	    rcvFilter_wilcox <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_wilcox$selectedFeaturesSet,...);
	}
    cStats <- predictionStats_binary(rcvFilter_wilcox$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
        
	if (is.null(fmeth_0))
	{
		cat("Kendall Feature Selection: ");
	    rcvFilter_kendall <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_correlation,featureSelection.control = list(method = "kendall",limit = 0.9,thr = 0.975),...);
 	}
	else
	{
	    rcvFilter_kendall <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_kendall$selectedFeaturesSet,...);
	}
   cStats <- predictionStats_binary(rcvFilter_kendall$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
	if (is.null(fmeth_0))
	{
		cat("Classic mRMR Feature Selection: ");
	    rcvFilter_mRMR <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,...);
 	}
	else
	{
	    rcvFilter_mRMR <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_mRMR$selectedFeaturesSet,...);
	}
    cStats <- predictionStats_binary(rcvFilter_mRMR$testPredictions,"", center = center);
    accciTable_filter <- rbind(accciTable_filter,cStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,cStats$berror)
    aucTable_filter <- rbind(aucTable_filter,cStats$aucs)
    senTable_filter <- rbind(senTable_filter,cStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,cStats$specificity)
	cindexTable_filter <- rbind(cindexTable_filter,cStats$cIndexCI)
    
    result <- list(errorciTable_filter = errorciTable_filter,
                   accciTable_filter = accciTable_filter,
                   aucTable_filter = aucTable_filter,
				   senTable_filter = senTable_filter,
				   speTable_filter = speTable_filter,
				   cindexTable_filter = cindexTable_filter,
				   rcvFilter_reference = rcvFilter_reference,
				   rcvFilter_LASSO = rcvFilter_LASSO,
				   rcvFilter_RPART = rcvFilter_RPART,
                   rcvFilter_IDI = rcvFilter_IDI,
                   rcvFilter_NRI = rcvFilter_NRI,
                   rcvFilter_RF = rcvFilter_RF,
                   rcvFilter_tStudent = rcvFilter_tStudent,
                   rcvFilter_wilcox = rcvFilter_wilcox,
                   rcvFilter_kendall = rcvFilter_kendall,
                   rcvFilter_mRMR = rcvFilter_mRMR
                   )
    
    return(result);
  }
  
  
  
######################Classification Algorithms####################################  
  
  
  if (is.null(referenceCV))
  {
 	  cat("Modeling BSWiMS: + Model found, - No Model \n"); 
	  referenceCV <- randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
	  referenceName = "BSWiMS";
	  referenceFilterName = "BSWiMS";
  }
  else
  {
		reps <- referenceCV$repetitions;
  }
  cStats <- predictionStats_binary(referenceCV$testPredictions,plotname = referenceName,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)

  rcvRF <- randomCV(theData,theOutcome,randomForest::randomForest,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
  cStats <- predictionStats_binary(rcvRF$testPredictions,plotname = "Random Forest",center = TRUE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)
  
  rcvRPART <- randomCV(theData,theOutcome,rpart::rpart,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
  cStats <- predictionStats_binary(rcvRPART$testPredictions,plotname = "RPART",center = TRUE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)
  
  rcvLASSO <- randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",family = "binomial");
  cStats <- predictionStats_binary(rcvLASSO$testPredictions,plotname = "LASSO",center = FALSE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)
  
  rcvSVM <- randomCV(theData,theOutcome,e1071::svm,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,asFactor=TRUE);
  cStats <- predictionStats_binary(rcvSVM$testPredictions,plotname = "SVM",center = TRUE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)
  
  rcvKNN <- randomCV(theData,theOutcome,KNN_method,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,scaleMethod = "Order");
  cStats <- predictionStats_binary(rcvKNN$testPredictions,plotname = "KNN",center = TRUE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)

# Method Meta Ensemble  

  reftest <- referenceCV$medianTest[,2];
  if (min(reftest) < 0) reftest <- 1.0*(reftest > 0);
  lasstest <- rcvLASSO$medianTest[,2];
  if (min(lasstest) < 0) lasstest <- 1.0*(lasstest > 0);
  ens <- cbind(referenceCV$medianTest[,1],rowMeans(cbind(reftest,lasstest,rcvRF$medianTest[,2],rcvKNN$medianTest[,2],rcvSVM$medianTest[,2])));
  cStats <- predictionStats_binary(ens,plotname = "Ensemble",center = TRUE,cex=0.8);
  accciTable <- rbind(accciTable,cStats$accc)
  errorciTable <- rbind(errorciTable,cStats$berror)
  aucTable <- rbind(aucTable,cStats$aucs)
  senTable <- rbind(senTable,cStats$sensitivity)
  speTable <- rbind(speTable,cStats$specificity)
  cidxTable <- rbind(cidxTable,cStats$cIndexCI)
  
######################Filters  ####################################  
  
  cat("KNN\n")
  fmeth_0 <- FilterMethod(KNN_method,"KNN",center = TRUE,scaleMethod = "Order")
  fmeth <- fmeth_0;
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);

  cat("Naive Bayes\n")
  fmeth <- FilterMethod(NAIVE_BAYES,"Naive Bayes",center = TRUE,asFactor = TRUE,usekernel = TRUE)
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);
  
  cat("Filtered Signature RSS\n")
  fmeth <- FilterMethod(CVsignature,"Nearest Centroid (RSS)",center = FALSE,method = "RSS")
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);

  cat("Filtered Signature Spearman\n")
  fmeth <- FilterMethod(CVsignature,"Nearest Centroid (Spearman)",center = FALSE,method = "spearman")
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);
  
  cat("Filtered RF\n")
  fmeth <- FilterMethod(randomForest::randomForest,"RF",center = TRUE,asFactor=TRUE)
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);

  cat("Filtered SVM\n")
  
  fmeth <- FilterMethod(e1071::svm,"SVM",center = TRUE,asFactor=TRUE)
  aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
  accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
  errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
  senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
  speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
  cindexTable_filter  <- rbind(cindexTable_filter,fmeth$cindexTable_filter);

	test_Predictions <- referenceCV$medianTest;
	tnames <- rownames(test_Predictions)
	test_Predictions[,2] <- 1.0/(1.0+exp(-test_Predictions[,2]));
	test_Predictions <- cbind(test_Predictions,rcvRF$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvLASSO$medianTest[tnames,2])
	test_Predictions[,4] <- 1.0/(1.0+exp(-test_Predictions[,4]));
	test_Predictions <- cbind(test_Predictions,rcvRPART$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvKNN$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvSVM$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,ens[tnames,2])
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_reference$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_LASSO$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RPART$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RF$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_IDI$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_NRI$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_tStudent$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_wilcox$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_kendall$medianTest[tnames,2]);

	
	colnames(test_Predictions) <- c("Outcome",referenceName,"RF","LASSO","RPART","KNN","SVM.mRMR","ENS",
	paste("SVM.",referenceFilterName,sep=""),"SVM.LASSO","SVM.RPART","SVM.RF","SVM.IDI","SVM.NRI","SVM.tStudent","SVM.Wilcox","SVM.Kendall");
	test_Predictions <- as.data.frame(test_Predictions)

	thesets <- c("Classifier Algorithm")
	theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","KNN","ENS")

	theClassMethod <- c("KNN","Naive Bayes","NC RSS","NC Spearman","RF","SVM")
	theFiltersets <- c(referenceFilterName,"LASSO","RPART","RF.ref","IDI","NRI","t-test","Wilcoxon","Kendall","mRMR")

	
	rownames(accciTable) <- theMethod;
	rownames(errorciTable) <- theMethod;
	rownames(aucTable) <- theMethod;
	rownames(senTable) <- theMethod;
	rownames(speTable) <- theMethod;
	rownames(cidxTable) <- theMethod;

    ff <- names(referenceCV$featureFrequency)
	ff <- c(ff,names(rcvLASSO$featureFrequency))
	ff <- c(ff,names(rcvRPART$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_RF$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_IDI$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_NRI$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_tStudent$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_wilcox$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_kendall$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_mRMR$featureFrequency))
	ff <- unique(ff)

#	Nvar <- min(c(1000,length(ff)))
#	selFrequency <- matrix(0,nrow = Nvar,ncol = length(theFiltersets))
#	rownames(selFrequency) <- ff[1:Nvar]

	Nvar <- length(ff);
	selFrequency <- matrix(0,nrow = Nvar,ncol = length(theFiltersets))
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
	ff <- fmeth$rcvFilter_IDI$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"IDI"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_NRI$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"NRI"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_wilcox$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"Wilcoxon"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_tStudent$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"t-test"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_kendall$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"Kendall"] <- ff[selnames[fnames]]
	ff <- fmeth$rcvFilter_mRMR$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"mRMR"] <- ff[selnames[fnames]]
	selFrequency <- selFrequency/reps

	elapcol <- names(referenceCV$theTimes) == "elapsed"
	theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","KNN","ENS")
	cputimes <- list(Reference = mean(referenceCV$theTimes[ elapcol ]),RF = mean(rcvRF$theTimes[ elapcol ]),RPART = mean(rcvRPART$theTimes[ elapcol ]),LASSO = mean(rcvLASSO$theTimes[ elapcol ]),SVM = mean(rcvSVM$theTimes[ elapcol ]),KNN = mean(rcvKNN$theTimes[ elapcol ]))

	
	jaccard_filter = list(Reference = referenceCV$jaccard,
                      LASSO = rcvLASSO$jaccard,
                      rpart = rcvRPART$jaccard,
                      RF = fmeth$rcvFilter_RF$jaccard,
                      IDI = fmeth$rcvFilter_IDI$jaccard,
                      NRI = fmeth$rcvFilter_NRI$jaccard,
                      ttest = fmeth$rcvFilter_tStudent$jaccard,
                      wtest = fmeth$rcvFilter_wilcox$jaccard,
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

	
  result <- list(errorciTable = errorciTable,accciTable = accciTable,aucTable = aucTable,senTable = senTable,speTable = speTable,cidxTable=cidxTable,
                 errorciTable_filter = errorciTable_filter,accciTable_filter = accciTable_filter,aucTable_filter = aucTable_filter,
				 senciTable_filter = senciTable_filter,speciTable_filter = speciTable_filter,cindexTable_filter=cindexTable_filter,
                 times = list(Reference = referenceCV$theTimes,RF = rcvRF$theTimes,rpart = rcvRPART$theTimes,LASSO = rcvLASSO$theTimes,SVM = rcvSVM$theTimes,KNN = rcvKNN$theTimes),
				 jaccard = jaccard,
				 featsize = featsize,
                 TheCVEvaluations = list(Reference = referenceCV,
                                         RF = rcvRF,
                                         LASSO = rcvLASSO,
                                         RPART = rcvRPART,
                                         KNN = rcvKNN,
                                         SVM = rcvSVM,
                                         FRF = fmeth$rcvFilter_RF,
                                         IDI = fmeth$rcvFilter_IDI,
                                         NRI = fmeth$rcvFilter_NRI,
                                         ttest = fmeth$rcvFilter_tStudent,
                                         wtest = fmeth$rcvFilter_wilcox,
                                         kendall = fmeth$rcvFilter_kendall,
                                         mRMR = fmeth$rcvFilter_mRMR
                                         ),
				 thesets = thesets,
				 theMethod = theMethod,
				 theClassMethod = theClassMethod,
				 theFiltersets = theFiltersets,
				 testPredictions = test_Predictions,
				 featureSelectionFrequency = selFrequency,
				 cpuElapsedTimes=cputimes				 
                )
  class(result) <- c("FRESA_benchmark","Binary");
  return(result)
}


