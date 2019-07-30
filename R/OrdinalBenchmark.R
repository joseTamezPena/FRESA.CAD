OrdinalBenchmark <-  function(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,referenceCV = NULL,referenceName = "Reference",referenceFilterName="Reference")
{

  if (!requireNamespace("DescTools", quietly = TRUE)) {
	  install.packages("DescTools", dependencies = TRUE)
	  }
  if (!requireNamespace("irr", quietly = TRUE)) {
	  install.packages("irr", dependencies = TRUE)
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

	  
  BMAETable <- NULL 
  KappaTable <- NULL
  BiasTable <- NULL
  KendallTable <- NULL
  AUCTable <- NULL;
  SENTable <- NULL;
  ACCTable <- NULL;

  
  
  BMAETable_filter <- NULL 
  KappaTable_filter <- NULL
  BiasTable_filter <- NULL
  KendallTable_filter <- NULL
  AUCTable_filter <- NULL;
  SENTable_filter <- NULL;
  ACCTable_filter <- NULL;
  fmeth_0 <- NULL;
#  par(mfrow = c(1,1));

  FilterMethod <-  function(ordFun = e1071::svm, classname = "",...)
  {
    BMAETable_filter <- NULL 
    KappaTable_filter <- NULL
    BiasTable_filter <- NULL
    KendallTable_filter <- NULL
	AUCTable_filter <- NULL    
	SENTable_filter <- NULL    
	ACCTable_filter <- NULL    
	
    parm <- list(...);
    rcvFilter_REFERENCE <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,...);
    sta <- predictionStats_ordinal(rcvFilter_REFERENCE$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    

    rcvFilter_LASSO <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvLASSO$selectedFeaturesSet,...);
    sta <- predictionStats_ordinal(rcvFilter_LASSO$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    
	rcvFilter_RPART <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvRPART$selectedFeaturesSet,...);
    sta <- predictionStats_ordinal(rcvFilter_RPART$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    
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
        warning ("Less than 2, then will keep the top five of RF\n")
        if (length(selectedFeaturesSet[[i]]) > 5)
        {
          selectedFeaturesSet[[i]] <- selectedFeaturesSet[[i]][1:5];
        }
      }
      #      cat("RF:", length(selectedFeaturesSet[[i]])," BSWIMS:",length(referenceCV$selectedFeaturesSet[[i]]),"\n")
    }

	
    rcvFilter_RF <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = selectedFeaturesSet,...);
    sta <- predictionStats_ordinal(rcvFilter_RF$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    
	if (is.null(fmeth_0))
	{
		rcvFilter_Ft <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_residual,featureSelection.control = list(uniTest = "Ftest",limit = 0.33,thr = 0.975),...);
	}
	else
	{
		rcvFilter_Ft <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_Ft$selectedFeaturesSet,...);
	}
    sta <- predictionStats_ordinal(rcvFilter_Ft$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    
        
	if (is.null(fmeth_0))
	{
	    rcvFilter_kendall <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_correlation,featureSelection.control = list(method = "kendall",limit = 0.33,thr = 0.975),...);
	}
	else
	{
		rcvFilter_kendall <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_kendall$selectedFeaturesSet,...);
	}
    sta <- predictionStats_ordinal(rcvFilter_kendall$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);

	if (is.null(fmeth_0))
	{
	    rcvFilter_mRMR <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,...);
	}
	else
	{
		rcvFilter_mRMR <- randomCV(theData,theOutcome,ordFun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = fmeth_0$rcvFilter_mRMR$selectedFeaturesSet,...);
	}
    sta <- predictionStats_ordinal(rcvFilter_mRMR$medianTest);
    BMAETable_filter <- rbind(BMAETable_filter,sta$BMAE);
    KappaTable_filter <- rbind(KappaTable_filter,sta$Kapp);
    BiasTable_filter <- rbind(BiasTable_filter,sta$Bias);
    KendallTable_filter <- rbind(KendallTable_filter,sta$Kendall);
	AUCTable_filter <- rbind(AUCTable_filter,sta$class95ci$aucci);
	SENTable_filter <- rbind(SENTable_filter,sta$class95ci$senci);
	ACCTable_filter <- rbind(ACCTable_filter,sta$class95ci$accci);
    
    result <- list(BMAETable_filter = BMAETable_filter,
                   KappaTable_filter = KappaTable_filter,
                   BiasTable_filter = BiasTable_filter,
                   KendallTable_filter = KendallTable_filter,
				   AUCTable_filter = AUCTable_filter,
				   SENTable_filter = SENTable_filter,
				   ACCTable_filter = ACCTable_filter,
				   rcvFilter_REFERENCE = rcvFilter_REFERENCE,
				   rcvFilter_LASSO = rcvFilter_LASSO,
				   rcvFilter_RPART = rcvFilter_RPART,
                   rcvFilter_RF = rcvFilter_RF,
                   rcvFilter_Ft = rcvFilter_Ft,
                   rcvFilter_kendall = rcvFilter_kendall,
                   rcvFilter_mRMR = rcvFilter_mRMR)
    
    return(result);
  }
  
  
  
  ######################Clasification Algorithms####################################  
  
  if (is.null(referenceCV))
  {
	  referenceCV <- randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
	  referenceName = "BSWiMS";
	  referenceFilterName = "BSWiMS";
  }
  else
  {
		reps <- referenceCV$repetitions;
  }
  
  sta <- predictionStats_ordinal(referenceCV$medianTest,referenceName);
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);

  rcvRF <- randomCV(theData,theOutcome,randomForest::randomForest,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
  sta <- predictionStats_ordinal(rcvRF$medianTest,"Random Forest");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  rcvRPART <- randomCV(theData,theOutcome,rpart::rpart,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",asFactor = TRUE);
  sta <- predictionStats_ordinal(rcvRPART$medianTest,"RPART");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  rcvLASSO <- randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",family = "multinomial");
  sta <- predictionStats_ordinal(rcvLASSO$medianTest,"LASSO");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  rcvSVM <- randomCV(theData,theOutcome,e1071::svm,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = mRMR.classic_FRESA,asFactor = TRUE);
  sta <- predictionStats_ordinal(rcvSVM$medianTest,"SVM");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  #  Filtered KNN  
  rcvKNN <- randomCV(theData,theOutcome,KNN_method,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet,scaleMethod = "Order");
  sta <- predictionStats_ordinal(rcvKNN$medianTest,"KNN");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  #  Method Meta Ensemble  
  ens <- cbind(referenceCV$medianTest[,1],rowMedians(cbind(referenceCV$medianTest[,2],rcvLASSO$medianTest[,2],rcvRF$medianTest[,2],rcvKNN$medianTest[,2],rcvSVM$medianTest[,2])))
  sta <- predictionStats_ordinal(ens,"Ensemble");
  BMAETable <- rbind(BMAETable,sta$BMAE);
  KappaTable <- rbind(KappaTable,sta$Kapp);
  BiasTable <- rbind(BiasTable,sta$Bias);
  KendallTable <- rbind(KendallTable,sta$Kendall);
  AUCTable <- rbind(AUCTable,sta$class95ci$aucci);
  SENTable <- rbind(SENTable,sta$class95ci$senci);
  ACCTable <- rbind(ACCTable,sta$class95ci$accci);
  
  ######################Filters for SVM ####################################  

  cat("Ordinal Logit\n")
  fmeth <- FilterMethod(MASS::polr,"Ordinal Logit",asFactor=TRUE)
  fmeth_0 <- fmeth
  BMAETable_filter <- rbind(BMAETable_filter,fmeth$BMAETable_filter);
  KappaTable_filter  <- rbind(KappaTable_filter,fmeth$KappaTable_filter);
  BiasTable_filter  <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
  KendallTable_filter  <- rbind(KendallTable_filter,fmeth$KendallTable_filter);
  AUCTable_filter <- rbind(AUCTable_filter,fmeth$AUCTable_filter);
  SENTable_filter <- rbind(SENTable_filter,fmeth$SENTable_filter);
  ACCTable_filter <- rbind(ACCTable_filter,fmeth$ACCTable_filter);

  cat("KNN\n")
  fmeth <- FilterMethod(KNN_method,"KNN",scaleMethod = "Order")
  BMAETable_filter <- rbind(BMAETable_filter,fmeth$BMAETable_filter);
  KappaTable_filter  <- rbind(KappaTable_filter,fmeth$KappaTable_filter);
  BiasTable_filter  <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
  KendallTable_filter  <- rbind(KendallTable_filter,fmeth$KendallTable_filter);
  AUCTable_filter <- rbind(AUCTable_filter,fmeth$AUCTable_filter);
  SENTable_filter <- rbind(SENTable_filter,fmeth$SENTable_filter);
  ACCTable_filter <- rbind(ACCTable_filter,fmeth$ACCTable_filter);

  cat("Naive Bayes\n")
  fmeth <- FilterMethod(NAIVE_BAYES,"Naive Bayes",asFactor = TRUE,usekernel = TRUE)
  BMAETable_filter <- rbind(BMAETable_filter,fmeth$BMAETable_filter);
  KappaTable_filter  <- rbind(KappaTable_filter,fmeth$KappaTable_filter);
  BiasTable_filter  <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
  KendallTable_filter  <- rbind(KendallTable_filter,fmeth$KendallTable_filter);
  AUCTable_filter <- rbind(AUCTable_filter,fmeth$AUCTable_filter);
  SENTable_filter <- rbind(SENTable_filter,fmeth$SENTable_filter);
  ACCTable_filter <- rbind(ACCTable_filter,fmeth$ACCTable_filter);

  cat("RF\n")
  fmeth <- FilterMethod(randomForest::randomForest,"RF",asFactor = TRUE)
  BMAETable_filter <- rbind(BMAETable_filter,fmeth$BMAETable_filter);
  KappaTable_filter  <- rbind(KappaTable_filter,fmeth$KappaTable_filter);
  BiasTable_filter  <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
  KendallTable_filter  <- rbind(KendallTable_filter,fmeth$KendallTable_filter);
  AUCTable_filter <- rbind(AUCTable_filter,fmeth$AUCTable_filter);
  SENTable_filter <- rbind(SENTable_filter,fmeth$SENTable_filter);
  ACCTable_filter <- rbind(ACCTable_filter,fmeth$ACCTable_filter);

  cat("SVM\n")
  fmeth <- FilterMethod(e1071::svm,"SVM",asFactor = TRUE)
  BMAETable_filter <- rbind(BMAETable_filter,fmeth$BMAETable_filter);
  KappaTable_filter  <- rbind(KappaTable_filter,fmeth$KappaTable_filter);
  BiasTable_filter  <- rbind(BiasTable_filter,fmeth$BiasTable_filter);
  KendallTable_filter  <- rbind(KendallTable_filter,fmeth$KendallTable_filter);
  AUCTable_filter <- rbind(AUCTable_filter,fmeth$AUCTable_filter);
  SENTable_filter <- rbind(SENTable_filter,fmeth$SENTable_filter);
  ACCTable_filter <- rbind(ACCTable_filter,fmeth$ACCTable_filter);

  
  	test_Predictions <- referenceCV$medianTest;
	tnames <- rownames(test_Predictions)
	test_Predictions <- cbind(test_Predictions,rcvRF$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvLASSO$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvRPART$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvKNN$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,rcvSVM$medianTest[tnames,2])
	test_Predictions <- cbind(test_Predictions,ens[tnames,2])
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_REFERENCE$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_LASSO$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RPART$medianTest[tnames,2]);
	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_RF$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_Ft$medianTest[tnames,2]);
 	test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_kendall$medianTest[tnames,2]);

	
	colnames(test_Predictions) <- c("Outcome",referenceName,"RF","LASSO","RPART","KNN","SVM.mRMR","ENS",
	paste("SVM.",referenceFilterName,sep=""),"SVM.LASSO","SVM.RPART","SVM.RF","SVM.FT","SVM.Kendall");
	test_Predictions <- as.data.frame(test_Predictions)
  
  ff <- names(referenceCV$featureFrequency)
	ff <- c(ff,names(rcvLASSO$featureFrequency))
	ff <- c(ff,names(rcvRPART$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_RF$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_Ft$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_kendall$featureFrequency))
	ff <- c(ff,names(fmeth$rcvFilter_mRMR$featureFrequency))
	ff <- unique(ff)

  theOrdinalMethod <- c("Ordinal","KNN","Naive Bayes","RF","SVM")
  theFiltersets <- c(referenceFilterName,"LASSO","RPART","RF.ref","F.Test","Kendall","mRMR")

#  Nvar <- min(c(1000,length(ff)))
#	selFrequency <- matrix(0,nrow = Nvar,ncol = length(theFiltersets))
#	rownames(selFrequency) <- names(rcvRF$featureFrequency)[1:Nvar]

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

	ff <- fmeth$rcvFilter_Ft$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"F.Test"] <- ff[selnames[fnames]]

	ff <- fmeth$rcvFilter_kendall$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"Kendall"] <- ff[selnames[fnames]]

	ff <- fmeth$rcvFilter_mRMR$featureFrequency
	fnames <- selnames %in% names(ff)
	selFrequency[fnames,"mRMR"] <- ff[selnames[fnames]]

	selFrequency <- selFrequency/reps


  thesets <- c("Ordinal Algorithm")
  theMethod <- c(referenceName,"RF","RPART","LASSO","SVM","KNN","ENS")

  	rownames(BMAETable) <- theMethod;
	rownames(KappaTable) <- theMethod;
	rownames(BiasTable) <- theMethod;
	rownames(KendallTable) <- theMethod;

	elapcol <- names(referenceCV$theTimes) == "elapsed"
	cputimes <- list(Reference = mean(referenceCV$theTimes[ elapcol ]),RF = mean(rcvRF$theTimes[ elapcol ]),RPART = mean(rcvRPART$theTimes[ elapcol ]),
	LASSO = mean(rcvLASSO$theTimes[ elapcol ]),SVM = mean(rcvSVM$theTimes[ elapcol ]),KNN = mean(rcvKNN$theTimes[ elapcol ]))

	jaccard_filter = list(reference = referenceCV$jaccard,
					   LASSO = rcvLASSO$jaccard,
					   rpart = rcvRPART$jaccard,
					   RF = fmeth$rcvFilter_RF$jaccard,
					   Ftest = fmeth$rcvFilter_Ft$jaccard,
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

  
  result <- list(BMAETable = BMAETable,KappaTable = KappaTable,
				 BiasTable = BiasTable,KendallTable = KendallTable,
				 AUCTable = AUCTable,SENTable = SENTable,ACCTable = ACCTable,
                 BMAETable_filter = BMAETable_filter,KappaTable_filter = KappaTable_filter,
				 BiasTable_filter = BiasTable_filter,KendallTable_filter = KendallTable_filter,
                 AUCTable_filter = AUCTable_filter,SENTable_filter = SENTable_filter,ACCTable_filter = ACCTable_filter,
                 times = list(Reference = referenceCV$theTimes,RF = rcvRF$theTimes,rpart = rcvRPART$theTimes,LASSO = rcvLASSO$theTimes,SVM = rcvSVM$theTimes,KNN = rcvKNN$theTimes), 
				 jaccard = jaccard,
				 featsize = featsize,
                 TheCVEvaluations = list(BSWIMS = referenceCV,
                                         RF = rcvRF,
                                         LASSO = rcvLASSO,
                                         RPART = rcvRPART,
                                         KNN = rcvKNN,
                                         FRF = fmeth$rcvFilter_RF,
                                         Ftest = fmeth$rcvFilter_Ft,
                                         kendall = fmeth$rcvFilter_kendall,
                                         mRMR = fmeth$rcvFilter_mRMR
                 ),
				thesets = thesets,
				theMethod = theMethod,
				theOrdinalMethod = theOrdinalMethod,
				theFiltersets = theFiltersets,
				testPredictions = test_Predictions,
				featureSelectionFrequency = selFrequency,
				cpuElapsedTimes=cputimes				 
  )
  
  class(result) <- c("FRESA_benchmark","Ordinal");

  return(result)
}
