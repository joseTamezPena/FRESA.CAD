BinaryBenchmark <-	function(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,referenceCV = NULL,referenceName = "Reference",referenceFilterName="Reference")
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
	
#	 par(mfrow = c(1,1));
	FilterMethod <-	function(clasfun = e1071::svm, classname = "", center = FALSE, ...)
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
	
	reftest <- numeric();
	theFiltersets <- character();
	theClassMethod <- character();
	elapcol <- character();
	cputimes <- list();
	jaccard <- NULL;
	featsize <- NULL;
	TheCVEvaluations = list();
	times <- list();
	jaccard_filter <- list();
	selFrequency <- data.frame(colnames(theData));
	rownames(selFrequency) <- colnames(theData);
	if (is.null(referenceCV))
	{
 		cat("Modeling BSWiMS: + Model found, - No Model \n"); 
		referenceCV <- randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
		referenceName = "BSWiMS";
		referenceFilterName = "BSWiMS";
	}
	if (inherits(referenceCV,"list"))
	{
		elapcol <- names(referenceCV[[1]]$theTimes) == "elapsed"
		TheCVEvaluations <- referenceCV;
		for (i in 1:length(referenceCV))
		{
			cStats <- predictionStats_binary(referenceCV[[i]]$testPredictions,plotname = names(referenceCV)[i],cex=0.8);
			accciTable <- rbind(accciTable,cStats$accc)
			errorciTable <- rbind(errorciTable,cStats$berror)
			aucTable <- rbind(aucTable,cStats$aucs)
			senTable <- rbind(senTable,cStats$sensitivity)
			speTable <- rbind(speTable,cStats$specificity)
			cidxTable <- rbind(cidxTable,cStats$cIndexCI)
			preftest <- referenceCV[[i]]$medianTest[,2];
			if ((min(preftest) < -0.1) || (max(preftest) > 1.1))
			{
				preftest <- 1.0/(1.0 + exp(-preftest));
			}
			reftest <- cbind(reftest,preftest);
			cputimes[[i]] = mean(referenceCV[[i]]$theTimes[ elapcol ]);
			times[[i]] <- referenceCV[[i]]$theTimes;
			selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
			selFrequency[names(referenceCV[[i]]$featureFrequency),ncol(selFrequency)] <- referenceCV[[i]]$featureFrequency;
			jaccard_filter[[i]] <- referenceCV[[i]]$jaccard;
		}
		referenceName <- names(referenceCV);
		referenceFilterName <- paste("FS",names(referenceCV),sep="_");
		referenceCV <- referenceCV[[1]];
		class(referenceCV) <- "list"
	}
	else
	{
		cStats <- predictionStats_binary(referenceCV$testPredictions,plotname = referenceName,cex=0.8);
		accciTable <- rbind(accciTable,cStats$accc)
		errorciTable <- rbind(errorciTable,cStats$berror)
		aucTable <- rbind(aucTable,cStats$aucs)
		senTable <- rbind(senTable,cStats$sensitivity)
		speTable <- rbind(speTable,cStats$specificity)
		cidxTable <- rbind(cidxTable,cStats$cIndexCI)
		reftest <- referenceCV$medianTest[,2];
		if ((min(reftest) < -0.1) || (max(reftest) > 1.1))
		{
			reftest <- 1.0/(1.0 + exp(-reftest));
		}
		TheCVEvaluations$Reference <- referenceCV;
		times[[1]] <- referenceCV$theTimes;
		elapcol <- names(referenceCV$theTimes) == "elapsed"
		cputimes[[1]] = mean(referenceCV$theTimes[ elapcol ]);
		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(referenceCV$featureFrequency),ncol(selFrequency)] <- referenceCV$featureFrequency;
		jaccard_filter[[1]] <- referenceCV$jaccard;
	}
	reps <- referenceCV$repetitions;

	rcvRF <- randomCV(theData,theOutcome,randomForest::randomForest,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
	cStats <- predictionStats_binary(rcvRF$testPredictions,plotname = "Random Forest",center = TRUE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)
	TheCVEvaluations$RF <- rcvRF;
	times$RF <- rcvRF$theTimes
	selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
	freq <- rcvRF$featureFrequency[names(rcvRF$featureFrequency) %in% rownames(selFrequency)];
	selFrequency[names(freq),ncol(selFrequency)] <- freq;
	theFiltersets <- c(referenceFilterName,"RF");
	jaccard_filter$RF <- rcvRF$jaccard;

	
	rcvRPART <- randomCV(theData,theOutcome,rpart::rpart,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
	cStats <- predictionStats_binary(rcvRPART$testPredictions,plotname = "RPART",center = TRUE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)
	TheCVEvaluations$RPART <- rcvRPART;
	times$RPART <- rcvRPART$theTimes
	selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
	selFrequency[names(rcvRPART$featureFrequency),ncol(selFrequency)] <- rcvRPART$featureFrequency;
	theFiltersets <- c(theFiltersets,"RPART");
	jaccard_filter$RPART <- rcvRPART$jaccard;


	rcvLASSO <- randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",family = "binomial");
	cStats <- predictionStats_binary(rcvLASSO$testPredictions,plotname = "LASSO",center = FALSE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)
	TheCVEvaluations$LASSO <- rcvLASSO;
	times$LASSO <- rcvLASSO$theTimes
	selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
	selFrequency[names(rcvLASSO$featureFrequency),ncol(selFrequency)] <- rcvLASSO$featureFrequency;
	theFiltersets <- c(theFiltersets,"LASSO_MIN");
	jaccard_filter$LASSO <- rcvLASSO$jaccard;
	
	rcvSVM <- randomCV(theData,theOutcome,e1071::svm,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,asFactor=TRUE,probability = TRUE);
	cStats <- predictionStats_binary(rcvSVM$testPredictions,plotname = "SVM",center = TRUE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)
	TheCVEvaluations$SVM <- rcvSVM;
	times$SVM <- rcvSVM$theTimes
	selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
	selFrequency[names(rcvSVM$featureFrequency),ncol(selFrequency)] <- rcvSVM$featureFrequency;
	theFiltersets <- c(theFiltersets,"mRMR.classic");
	jaccard_filter$mRMR <- rcvSVM$jaccard;

	rcvKNN <- randomCV(theData,theOutcome,KNN_method,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,scaleMethod = "Order");
	cStats <- predictionStats_binary(rcvKNN$testPredictions,plotname = "KNN",center = TRUE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)
	TheCVEvaluations$KNN <- rcvKNN;
	times$KNN <- rcvKNN$theTimes;

# Method Meta Ensemble	

	lasstest <- rcvLASSO$medianTest[,2];
	if ((min(lasstest) < 0) || (max(lasstest) > 1.1)) 
	{
		lasstest <- 1.0/(1.0 + exp(-lasstest));
	}

	ens <- cbind(referenceCV$medianTest[,1],rowMeans(cbind(reftest,lasstest,rcvRF$medianTest[,2],rcvKNN$medianTest[,2],rcvSVM$medianTest[,2])));
	cStats <- predictionStats_binary(ens,plotname = "Ensemble",center = TRUE,cex=0.8);
	accciTable <- rbind(accciTable,cStats$accc)
	errorciTable <- rbind(errorciTable,cStats$berror)
	aucTable <- rbind(aucTable,cStats$aucs)
	senTable <- rbind(senTable,cStats$sensitivity)
	speTable <- rbind(speTable,cStats$specificity)
	cidxTable <- rbind(cidxTable,cStats$cIndexCI)

	test_Predictions <- cbind(referenceCV$medianTest[,1],reftest);
	tnames <- rownames(referenceCV$medianTest)
	test_Predictions <- cbind(test_Predictions,rcvRF$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvLASSO$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvRPART$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvKNN$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvSVM$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,ens[tnames,2])

	colnames(test_Predictions) <- c("Outcome",referenceName,"RF","LASSO","RPART","KNN","SVM.mRMR","ENS")

	thesets <- c("Classifier Algorithm")
	theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","KNN","ENS")


	
	rownames(accciTable) <- theMethod;
	rownames(errorciTable) <- theMethod;
	rownames(aucTable) <- theMethod;
	rownames(senTable) <- theMethod;
	rownames(speTable) <- theMethod;
	rownames(cidxTable) <- theMethod;

	cputimes$RF = mean(rcvRF$theTimes[ elapcol ])
	cputimes$RPART = mean(rcvRPART$theTimes[ elapcol ])
	cputimes$LASSO = mean(rcvLASSO$theTimes[ elapcol ])
	cputimes$SVM = mean(rcvSVM$theTimes[ elapcol ])
	cputimes$KNN = mean(rcvKNN$theTimes[ elapcol ])
	cputimes <- unlist(cputimes);
	cputimes <- c(cputimes,sum(cputimes)-mean(rcvRPART$theTimes[ elapcol ]));
	names(cputimes) <- theMethod;

######################Filters	####################################	
	
	if (!inherits(referenceCV,"list"))
	{
		classnames <- colnames(test_Predictions);
		cat("KNN\n")
		fmeth_0 <- FilterMethod(KNN_method,"KNN",center = TRUE,scaleMethod = "Order")
		fmeth <- fmeth_0;
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);

		cat("Naive Bayes\n")
		fmeth <- FilterMethod(NAIVE_BAYES,"Naive Bayes",center = TRUE,asFactor = TRUE,usekernel = TRUE)
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);
		
		cat("Filtered Signature RSS\n")
		fmeth <- FilterMethod(CVsignature,"Nearest Centroid (RSS)",center = FALSE,method = "RSS")
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);

		cat("Filtered Signature Spearman\n")
		fmeth <- FilterMethod(CVsignature,"Nearest Centroid (Spearman)",center = FALSE,method = "spearman")
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);
		
		cat("Filtered RF\n")
		fmeth <- FilterMethod(randomForest::randomForest,"RF",center = TRUE,asFactor=TRUE)
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);

		cat("Filtered SVM\n")
		
		fmeth <- FilterMethod(e1071::svm,"SVM",center = TRUE,asFactor=TRUE,probability = TRUE)
		aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter);
		accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter);
		errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter);
		senciTable_filter <- rbind(senciTable_filter,fmeth$senTable_filter);
		speciTable_filter <- rbind(speciTable_filter,fmeth$speTable_filter);
		cindexTable_filter	<- rbind(cindexTable_filter,fmeth$cindexTable_filter);
		TheCVEvaluations$IDI.SVM <- fmeth$rcvFilter_IDI;
		TheCVEvaluations$NRI.SVM <- fmeth$rcvFilter_NRI;
		TheCVEvaluations$tStudent.SVM <- fmeth$rcvFilter_tStudent;
		TheCVEvaluations$wilcox.SVM <- fmeth$rcvFilter_wilcox;
		TheCVEvaluations$kendall.SVM <- fmeth$rcvFilter_kendall;


		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_reference$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_LASSO$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RPART$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RF$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_IDI$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_NRI$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_tStudent$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_wilcox$medianTest[tnames,2]);
		test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_kendall$medianTest[tnames,2]);
		colnames(test_Predictions) <-	c(classnames,paste("SVM.",referenceFilterName,sep=""),"SVM.LASSO","SVM.RPART","SVM.RF","SVM.IDI","SVM.NRI","SVM.tStudent","SVM.Wilcox","SVM.Kendall");
		theClassMethod <- c("KNN","Naive Bayes","NC RSS","NC Spearman","RF","SVM")

		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(fmeth$rcvFilter_IDI$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_IDI$featureFrequency;
		theFiltersets <- c(theFiltersets,"IDI");
		jaccard_filter$IDI <- fmeth$rcvFilter_IDI$jaccard;

		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(fmeth$rcvFilter_NRI$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_NRI$featureFrequency;
		theFiltersets <- c(theFiltersets,"NRI");
		jaccard_filter$NRI <- fmeth$rcvFilter_NRI$jaccard;

		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(fmeth$rcvFilter_tStudent$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_tStudent$featureFrequency;
		theFiltersets <- c(theFiltersets,"tStudent");
		jaccard_filter$tStudent <- fmeth$rcvFilter_tStudent$jaccard;

		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(fmeth$rcvFilter_wilcox$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_wilcox$featureFrequency;
		theFiltersets <- c(theFiltersets,"wilcox");
		jaccard_filter$wilcox <- fmeth$rcvFilter_wilcox$jaccard;

		selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
		selFrequency[names(fmeth$rcvFilter_kendall$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_kendall$featureFrequency;
		theFiltersets <- c(theFiltersets,"kendall");
		jaccard_filter$kendall <- fmeth$rcvFilter_kendall$jaccard;

	}

	featsize <- unlist(lapply(jaccard_filter, `[`, c('averageLength')))
	names(featsize) <- theFiltersets;
	jaccard <- unlist(lapply(jaccard_filter, `[`, c('Jaccard.SM')))
	names(jaccard) <- theFiltersets;
	selFrequency <- as.data.frame(selFrequency[,-1])
	selFrequency <- selFrequency/reps;
	colnames(selFrequency) <- theFiltersets;
	totsum <- apply(selFrequency,1,sum);
	selFrequency <- selFrequency[order(-totsum),];
	totsum <- totsum[order(-totsum)];
	selFrequency <- selFrequency[totsum>0,];

	test_Predictions <- as.data.frame(test_Predictions)

	
	result <- list(errorciTable = errorciTable,accciTable = accciTable,aucTable = aucTable,senTable = senTable,speTable = speTable,cidxTable=cidxTable,
				  errorciTable_filter = errorciTable_filter,accciTable_filter = accciTable_filter,aucTable_filter = aucTable_filter,
				 senciTable_filter = senciTable_filter,speciTable_filter = speciTable_filter,cindexTable_filter=cindexTable_filter,
				 times = times,
				 jaccard = jaccard,
				 featsize = featsize,
				 TheCVEvaluations = TheCVEvaluations,
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


