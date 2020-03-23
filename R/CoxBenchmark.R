CoxBenchmark <-  function(theData = NULL, theOutcome = "Class", reps = 100, trainFraction = 0.5,referenceCV = NULL,referenceName = "Reference",referenceFilterName="COX.BSWiMS")
{
  if (!requireNamespace("BeSS", quietly = TRUE)) {
    install.packages("BeSS", dependencies = TRUE)
  }
  if (!requireNamespace("survminer", quietly = TRUE)) {
    install.packages("survminer", dependencies = TRUE)
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
  
  aucTable <- NULL 
  accciTable <- NULL
  errorciTable <- NULL
  senTable <- NULL
  speTable <- NULL
  
  aucTable_filter <- NULL 
  accciTable_filter <- NULL
  errorciTable_filter <- NULL
  senTable_filter <- NULL
  speTable_filter <- NULL
  
  CIFollowUPTable <- NULL 
  CIRisksTable <- NULL 
  LogRankTable <- NULL
  
  CIFollowUPTable_filter <- NULL 
  CIRisksTable_filter <- NULL 
  LogRankTable_filter <- NULL
  fmeth_0 <- NULL;
  
  FilterMethod <-  function(clasfun = survival::coxph, classname = "", center = FALSE, ...)
  {
    #rcvFilter_reference <- cpFinal$TheCVEvaluations$Reference$testPredictions
    rcvFilter_reference <- FRESA.CAD::randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet);
    cStats <- predictionStats_survival(rcvFilter_reference$survMedianTest,plotname="Cox with BSWiMS");
    CIFollowUPTable_filter <- rbind(CIFollowUPTable_filter,cStats$CIFollowUp);
    CIRisksTable_filter <- rbind(CIRisksTable_filter,cStats$CIRisk);
    LogRankTable_filter <- rbind(LogRankTable_filter,cStats$LogRank);
    #Stats binary
    #cambiar a median antes de subir
    #rcvFilter_reference <- cpFinal$TheCVEvaluations$COX.Reference
    binaryPreds <- rcvFilter_reference$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"Cox with BSWiMS")
    accciTable_filter <- rbind(accciTable_filter,binaryStats$accc)
    errorciTable_filter <- rbind(errorciTable_filter,binaryStats$berror)
    aucTable_filter <- rbind(aucTable_filter,binaryStats$aucs)
    senTable_filter <- rbind(senTable_filter,binaryStats$sensitivity)
    speTable_filter <- rbind(speTable_filter,binaryStats$specificity)
    
    rcvFilter_LASSO <- try(FRESA.CAD::randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvLASSO$selectedFeaturesSet));
    if (!inherits(rcvLASSO, "try-error"))
	  {
      cStats <- predictionStats_survival(rcvFilter_LASSO$survMedianTest,plotname="Cox with LASSO");
      CIFollowUPTable_filter <- rbind(CIFollowUPTable_filter,cStats$CIFollowUp);
      CIRisksTable_filter <- rbind(CIRisksTable_filter,cStats$CIRisk);
      LogRankTable_filter <- rbind(LogRankTable_filter,cStats$LogRank);
      #Stats binary
      binaryPreds <- rcvFilter_LASSO$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
      binaryStats <- predictionStats_binary(binaryPreds,"Cox with Lasso")
      accciTable_filter <- rbind(accciTable_filter,binaryStats$accc)
      errorciTable_filter <- rbind(errorciTable_filter,binaryStats$berror)
      aucTable_filter <- rbind(aucTable_filter,binaryStats$aucs)
      senTable_filter <- rbind(senTable_filter,binaryStats$sensitivity)
      speTable_filter <- rbind(speTable_filter,binaryStats$specificity)
    }  
    
    rcvFilter_BESS <- try(FRESA.CAD::randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvBESS$selectedFeaturesSet));
    if (!inherits(rcvLASSO, "try-error"))
	  {
      cStats <- predictionStats_survival(rcvFilter_BESS$survMedianTest,plotname="Cox with BESS");
      CIFollowUPTable_filter <- rbind(CIFollowUPTable_filter,cStats$CIFollowUp);
      CIRisksTable_filter <- rbind(CIRisksTable_filter,cStats$CIRisk);
      LogRankTable_filter <- rbind(LogRankTable_filter,cStats$LogRank);
      
      #Stats binary
      binaryStats <- predictionStats_binary(rcvFilter_BESS$survMedianTest[,c("Outcome","LinearPredictorsMedian")],"Cox with BeSS")
      accciTable_filter <- rbind(accciTable_filter,binaryStats$accc)
      errorciTable_filter <- rbind(errorciTable_filter,binaryStats$berror)
      aucTable_filter <- rbind(aucTable_filter,binaryStats$aucs)
      senTable_filter <- rbind(senTable_filter,binaryStats$sensitivity)
      speTable_filter <- rbind(speTable_filter,binaryStats$specificity)
    }

    cat("Univariate cox Feature Selection: ");
    rcvFilter_UniCox <- try(FRESA.CAD::randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_cox));
    if (!inherits(rcvLASSO, "try-error"))
	  {
      cStats <- predictionStats_survival(rcvFilter_UniCox$survMedianTest,"Cox with Univariate cox Feature Selection");
      CIFollowUPTable_filter <- rbind(CIFollowUPTable_filter,cStats$CIFollowUp);
      CIRisksTable_filter <- rbind(CIRisksTable_filter,cStats$CIRisk);
      LogRankTable_filter <- rbind(LogRankTable_filter,cStats$LogRank);
      
      #Stats binary
      binaryStats <- predictionStats_binary(rcvFilter_UniCox$survMedianTest[,c("Outcome","LinearPredictorsMedian")],"Unicox")
      accciTable_filter <- rbind(accciTable_filter,binaryStats$accc)
      errorciTable_filter <- rbind(errorciTable_filter,binaryStats$berror)
      aucTable_filter <- rbind(aucTable_filter,binaryStats$aucs)
      senTable_filter <- rbind(senTable_filter,binaryStats$sensitivity)
      speTable_filter <- rbind(speTable_filter,binaryStats$specificity)
    }

    result <- list(CIFollowUPTable_filter = CIFollowUPTable_filter,
                   CIRisksTable_filter = CIRisksTable_filter,
                   LogRankTable_filter = LogRankTable_filter,
                   accciTable_filter = accciTable_filter,
                   errorciTable_filter = errorciTable_filter,
                   aucTable_filter = aucTable_filter,
                   senTable_filter = senTable_filter,
                   speTable_filter = speTable_filter,
                   rcvFilter_reference = rcvFilter_reference,
                   rcvFilter_LASSO = rcvFilter_LASSO,
                   rcvFilter_BESS = rcvFilter_BESS,
                   rcvFilter_UniCox = rcvFilter_UniCox
    )
    
    return(result);
  }
  
  
  ######################Classification Algorithms####################################  
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
    referenceCV <- FRESA.CAD::randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
    referenceName = "BSWiMS";
    referenceFilterName = "Cox.BSWiMS";
    methods <- c(referenceName);
    theFiltersets <- c(referenceName);
  }
  if (class(referenceCV) == "list")
  {
    elapcol <- names(referenceCV[[1]]$theTimes) == "elapsed"
    TheCVEvaluations <- referenceCV;
    for (i in 1:length(referenceCV))
    {
      cStats <- predictionStats_survival(referenceCV[[i]]$survMedianTest,plotname = names(referenceCV)[i]);
      CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
      CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
      LogRankTable <- rbind(LogRankTable,cStats$LogRank);
      #referenceCV <- cpFinal$TheCVEvaluations$Reference
      binaryPreds <- referenceCV[[i]]$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
      binaryStats <- predictionStats_binary(binaryPreds,"BSWiMS")
      accciTable <- rbind(accciTable,binaryStats$accc)
      errorciTable <- rbind(errorciTable,binaryStats$berror)
      aucTable <- rbind(aucTable,binaryStats$aucs);
      senTable <- rbind(senTable,binaryStats$sensitivity);
      speTable <- rbind(speTable,binaryStats$specificity);
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
    cStats <- predictionStats_survival(referenceCV$survMedianTest,plotname = referenceName);
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    binaryPreds <- referenceCV$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"BSWiMS")
    accciTable <- rbind(accciTable,binaryStats$accc)
    errorciTable <- rbind(errorciTable,binaryStats$berror)
    aucTable <- rbind(aucTable,binaryStats$aucs)
    senTable <- rbind(senTable,binaryStats$sensitivity)
    speTable <- rbind(speTable,binaryStats$specificity)
    TheCVEvaluations$Reference <- referenceCV;
    times[[1]] <- referenceCV$theTimes;
    elapcol <- names(referenceCV$theTimes) == "elapsed"
    cputimes[[1]] = mean(referenceCV$theTimes[ elapcol ]);
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(referenceCV$featureFrequency),ncol(selFrequency)] <- referenceCV$featureFrequency;
    jaccard_filter[[1]] <- referenceCV$jaccard;
  }
  reps <- referenceCV$repetitions;
  ######################Predictions union  #################################### 
  test_Predictions <- referenceCV$survMedianTest;
  tnames <- rownames(test_Predictions)

  # 1 - pchisq(cStats$LogRank$chisq, length(cStats$LogRank$n) - 1)
  ######################LASSO#################################### 
  rcvLASSO <- try(FRESA.CAD::randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self"));
  if (!inherits(rcvLASSO, "try-error"))
	{
      methods <- cbind(methods,"LASSO");
		  cStats <- predictionStats_survival(rcvLASSO$survMedianTest,plotname = "LASSO");
      CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
      CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
      LogRankTable <- rbind(LogRankTable,cStats$LogRank);
      
      #rcvLASSO <- cpFinal$TheCVEvaluations$LASSO
      binaryPreds <- rcvLASSO$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
      binaryStats <- predictionStats_binary(binaryPreds,"Lasso")
      accciTable <- rbind(accciTable,binaryStats$accc)
      errorciTable <- rbind(errorciTable,binaryStats$berror)
      aucTable <- rbind(aucTable,binaryStats$aucs)
      senTable <- rbind(senTable,binaryStats$sensitivity)
      speTable <- rbind(speTable,binaryStats$specificity)
      TheCVEvaluations$LASSO <- rcvLASSO;
      times$LASSO <- rcvLASSO$theTimes
      selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
      selFrequency[names(rcvLASSO$featureFrequency),ncol(selFrequency)] <- rcvLASSO$featureFrequency;
      theFiltersets <- c(theFiltersets,"LASSO");
      jaccard_filter$LASSO <- rcvLASSO$jaccard;
      test_Predictions <- cbind(test_Predictions,rcvLASSO$survMedianTest[tnames,3],rcvLASSO$survMedianTest[tnames,4],rcvLASSO$survMedianTest[tnames,5],rcvLASSO$survMedianTest[tnames,6])
      cputimes$LASSO = mean(rcvLASSO$theTimes[ elapcol ])
	}

  ######################GLMNET_RIDGE#################################### 
  rcvGLMNET_RIDGE <- try(FRESA.CAD::randomCV(theData,theOutcome,GLMNET_RIDGE_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self"));
  if (!inherits(rcvGLMNET_RIDGE, "try-error"))
	{
    methods <- cbind(methods,"RIDGE");
    cStats <- predictionStats_survival(rcvGLMNET_RIDGE$survMedianTest,plotname = "GLMNET_RIDGE");
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    #rcvGLMNET_RIDGE <- cpFinal$TheCVEvaluations$GLMNET_RIDGE
    binaryPreds <- rcvGLMNET_RIDGE$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"Ridge")
    accciTable <- rbind(accciTable,binaryStats$accc)
    errorciTable <- rbind(errorciTable,binaryStats$berror)
    aucTable <- rbind(aucTable,binaryStats$aucs)
    senTable <- rbind(senTable,binaryStats$sensitivity)
    speTable <- rbind(speTable,binaryStats$specificity)
    TheCVEvaluations$RIDGE <- rcvGLMNET_RIDGE;
    times$RIDGE <- rcvGLMNET_RIDGE$theTimes
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(rcvGLMNET_RIDGE$featureFrequency),ncol(selFrequency)] <- rcvGLMNET_RIDGE$featureFrequency;
    theFiltersets <- c(theFiltersets,"RIDGE");
    jaccard_filter$RIDGE <- rcvGLMNET_RIDGE$jaccard;
    test_Predictions <- cbind(test_Predictions,rcvGLMNET_RIDGE$survMedianTest[tnames,3],rcvGLMNET_RIDGE$survMedianTest[tnames,4],rcvGLMNET_RIDGE$survMedianTest[tnames,5],rcvGLMNET_RIDGE$survMedianTest[tnames,6])
    cputimes$RIDGE = mean(rcvGLMNET_RIDGE$theTimes[ elapcol ])
  }
  
  ######################GLMNET_ELASTICNET#################################### 
  rcvGLMNET_ELASTICNET <- try(FRESA.CAD::randomCV(theData,theOutcome,GLMNET_ELASTICNET_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self"));
  if (!inherits(rcvGLMNET_ELASTICNET, "try-error"))
	{
    methods <- cbind(methods,"ELASTICNET");
    cStats <- predictionStats_survival(rcvGLMNET_ELASTICNET$survMedianTest,plotname = "GLMNET_ELASTICNET");
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    #rcvGLMNET_ELASTICNET <- cpFinal$TheCVEvaluations$GLMNET_ELASTICNET
    binaryPreds <- rcvGLMNET_ELASTICNET$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"ElasticNet")
    accciTable <- rbind(accciTable,binaryStats$accc)
    errorciTable <- rbind(errorciTable,binaryStats$berror)
    aucTable <- rbind(aucTable,binaryStats$aucs)
    senTable <- rbind(senTable,binaryStats$sensitivity)
    speTable <- rbind(speTable,binaryStats$specificity)
    TheCVEvaluations$ELASTICNET <- rcvGLMNET_ELASTICNET;
    times$ELASTICNET <- rcvGLMNET_ELASTICNET$theTimes
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(rcvGLMNET_ELASTICNET$featureFrequency),ncol(selFrequency)] <- rcvGLMNET_ELASTICNET$featureFrequency;
    theFiltersets <- c(theFiltersets,"ELASTICNET");
    jaccard_filter$ELASTICNET <- rcvGLMNET_ELASTICNET$jaccard;
    test_Predictions <- cbind(test_Predictions,rcvGLMNET_ELASTICNET$survMedianTest[tnames,3],rcvGLMNET_ELASTICNET$survMedianTest[tnames,4],rcvGLMNET_ELASTICNET$survMedianTest[tnames,5],rcvGLMNET_ELASTICNET$survMedianTest[tnames,6])
    cputimes$ELASTICNET = mean(rcvGLMNET_ELASTICNET$theTimes[ elapcol ])
  }
  ######################BESS#################################### 
  rcvBESS <- try(FRESA.CAD::randomCV(theData,theOutcome,BESS,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",method="gsection"));
  if (!inherits(rcvBESS, "try-error"))
	{
    methods <- cbind(methods,"BESS");
    cStats <- predictionStats_survival(rcvBESS$survMedianTest,plotname = "BeSS");
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    #rcvBESS <- cpFinal$TheCVEvaluations$BESS
    #rcvBESS$survMedianTest[,"LinearPredictorsMedian"] <- -rcvBESS$survMedianTest[,"LinearPredictorsMedian"]
    binaryPreds <- rcvBESS$survMedianTest[,c("Outcome","LinearPredictorsMedian")];
    binaryStats <- predictionStats_binary(binaryPreds,"BeSS");
    accciTable <- rbind(accciTable,binaryStats$accc);
    errorciTable <- rbind(errorciTable,binaryStats$berror);
    aucTable <- rbind(aucTable,binaryStats$aucs);
    senTable <- rbind(senTable,binaryStats$sensitivity);
    speTable <- rbind(speTable,binaryStats$specificity);
    TheCVEvaluations$BESS <- rcvBESS;
    times$BESS <- rcvBESS$theTimes
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(rcvBESS$featureFrequency),ncol(selFrequency)] <- rcvBESS$featureFrequency;
    theFiltersets <- c(theFiltersets,"BESS");
    jaccard_filter$BESS <- rcvBESS$jaccard;
    test_Predictions <- cbind(test_Predictions,rcvBESS$survMedianTest[tnames,3],rcvBESS$survMedianTest[tnames,4],rcvBESS$survMedianTest[tnames,5],rcvBESS$survMedianTest[tnames,6])
    cputimes$BESS = mean(rcvBESS$theTimes[ elapcol ])
  }

  ######################BESS SEQUENTIAL#################################### 
  rcvBESSSequential <- try(FRESA.CAD::randomCV(theData,theOutcome,BESS,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self",method="sequential",ic.type="GIC"));
  if (!inherits(rcvBESSSequential, "try-error"))
	{
    methods <- cbind(methods,"BeSS.SEQUENTIAL");
    cStats <- predictionStats_survival(rcvBESSSequential$survMedianTest,plotname = "BeSS.SEQUENTIAL");
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    binaryPreds <- rcvBESSSequential$survMedianTest[,c("Outcome","LinearPredictorsMedian")];
    binaryStats <- predictionStats_binary(binaryPreds,"BeSS.SEQUENTIAL");
    accciTable <- rbind(accciTable,binaryStats$accc);
    errorciTable <- rbind(errorciTable,binaryStats$berror);
    aucTable <- rbind(aucTable,binaryStats$aucs);
    senTable <- rbind(senTable,binaryStats$sensitivity);
    speTable <- rbind(speTable,binaryStats$specificity);
    TheCVEvaluations$BESS.SEQUENTIAL <- rcvBESSSequential;
    times$BESS.SEQUENTIAL <- rcvBESSSequential$theTimes
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(rcvBESSSequential$featureFrequency),ncol(selFrequency)] <- rcvBESSSequential$featureFrequency;
    theFiltersets <- c(theFiltersets,"BESS.SEQUENTIAL");
    jaccard_filter$BESS.SEQUENTIAL <- rcvBESSSequential$jaccard;
    test_Predictions <- cbind(test_Predictions,rcvBESSSequential$survMedianTest[tnames,3],rcvBESSSequential$survMedianTest[tnames,4],rcvBESSSequential$survMedianTest[tnames,5],rcvBESSSequential$survMedianTest[tnames,6])
    cputimes$BESS.SEQUENTIAL = mean(rcvBESSSequential$theTimes[ elapcol ])
  }
  
  ######################BESS SEQUENTIAL BIC#################################### 
  rcvBESSSequentialBIC <- try(FRESA.CAD::randomCV(theData,theOutcome,BESS,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self"));
  if (!inherits(rcvBESSSequentialBIC, "try-error"))
	{
    methods <- cbind(methods,"BeSS.SEQUENTIAL.BIC");
    cStats <- predictionStats_survival(rcvBESSSequentialBIC$survMedianTest,plotname = "BeSS.SEQUENTIAL.BIC");
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    binaryPreds <- rcvBESSSequentialBIC$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"BeSS.SEQUENTIAL.BIC")
    accciTable <- rbind(accciTable,binaryStats$accc);
    errorciTable <- rbind(errorciTable,binaryStats$berror);
    aucTable <- rbind(aucTable,binaryStats$aucs);
    senTable <- rbind(senTable,binaryStats$sensitivity);
    speTable <- rbind(speTable,binaryStats$specificity);
    TheCVEvaluations$BESS.SEQUENTIAL.BIC <- rcvBESSSequentialBIC;
    times$BESS.SEQUENTIAL.BIC <- rcvBESSSequentialBIC$theTimes
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    selFrequency[names(rcvBESSSequentialBIC$featureFrequency),ncol(selFrequency)] <- rcvBESSSequentialBIC$featureFrequency;
    theFiltersets <- c(theFiltersets,"BESS.SEQUENTIAL.BIC");
    jaccard_filter$BESS.SEQUENTIAL.BIC <- rcvBESSSequentialBIC$jaccard;
    test_Predictions <- cbind(test_Predictions,rcvBESSSequentialBIC$survMedianTest[tnames,3],rcvBESSSequentialBIC$survMedianTest[tnames,4],rcvBESSSequentialBIC$survMedianTest[tnames,5],rcvBESSSequentialBIC$survMedianTest[tnames,6])
    cputimes$BESS.SEQUENTIAL.BIC = mean(rcvBESSSequentialBIC$theTimes[ elapcol ])
  }

  ######################Esemble#################################### 
  
 

  predictions <- c("MartinGale","LinearPredictors","FollowUpTimes","Risks");
  columnNamesMethods <- NULL;
  
  for(x in methods)
  {
    for(y in predictions)
    {
      columnNamesMethods <- cbind(columnNamesMethods,paste(x,y,sep=""))
    }
  }  

  colnames(test_Predictions) <- c("Times","Outcome",columnNamesMethods);
  thesets <- c("Survival Algorithm")
  theMethod <- methods;
  rownames(CIFollowUPTable) <- theMethod;
  rownames(CIRisksTable) <- theMethod;
  rownames(LogRankTable) <- theMethod;
  rownames(accciTable) <- theMethod;
  rownames(errorciTable) <- theMethod;
  rownames(aucTable) <- theMethod;
  rownames(senTable) <- theMethod;
  rownames(speTable) <- theMethod;
  

  cputimes <- unlist(cputimes);
  names(cputimes) <- theMethod;
  
  ######################Filters  #################################### 
  
  if (class(referenceCV) != "list")
  {
    classnames <- colnames(test_Predictions);
    cat("Cox\n")
    fmeth <- FilterMethod(survival::coxph,"Cox")
    CIFollowUPTable_filter <- rbind(CIFollowUPTable_filter,fmeth$CIFollowUPTable_filter);
    CIRisksTable_filter <- rbind(CIRisksTable_filter,fmeth$CIFollowUPTable_filter);
    LogRankTable_filter <- rbind(LogRankTable_filter,fmeth$LogRankTable_filter);
    accciTable_filter <- rbind(accciTable_filter,fmeth$accciTable_filter)
    errorciTable_filter <- rbind(errorciTable_filter,fmeth$errorciTable_filter)
    aucTable_filter <- rbind(aucTable_filter,fmeth$aucTable_filter)
    senTable_filter <- rbind(senTable_filter,fmeth$senTable_filter)
    speTable_filter <- rbind(speTable_filter,fmeth$speTable_filter)
    TheCVEvaluations$Cox.Reference <- fmeth$rcvFilter_reference;
    TheCVEvaluations$Cox.LASSO = fmeth$rcvFilter_LASSO;
    TheCVEvaluations$Cox.BESS = fmeth$rcvFilter_BESS;
    TheCVEvaluations$Cox.Unicox = fmeth$rcvFilter_UniCox;
    
    test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_reference$survMedianTest[tnames,3],fmeth$rcvFilter_reference$survMedianTest[tnames,4],fmeth$rcvFilter_reference$survMedianTest[tnames,5],fmeth$rcvFilter_reference$survMedianTest[tnames,6]);
    filters = c("COX.BSWiMS");
     if (!inherits(fmeth$rcvFilter_LASSO, "try-error"))
     {
        test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_LASSO$survMedianTest[tnames,3],fmeth$rcvFilter_LASSO$survMedianTest[tnames,4],fmeth$rcvFilter_LASSO$survMedianTest[tnames,5],fmeth$rcvFilter_LASSO$survMedianTest[tnames,6]);
        filters = c(filters,"COX.LASSO");
     }
      
    if (!inherits(fmeth$rcvFilter_BESS, "try-error"))
    {
      test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_BESS$survMedianTest[tnames,3],fmeth$rcvFilter_BESS$survMedianTest[tnames,4],fmeth$rcvFilter_BESS$survMedianTest[tnames,5],fmeth$rcvFilter_BESS$survMedianTest[tnames,6]);
      filters = c(filters,"COX.BESS");
    }
      
    if (!inherits(fmeth$rcvFilter_UniCox, "try-error"))
    {
      test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_UniCox$survMedianTest[tnames,3],fmeth$rcvFilter_UniCox$survMedianTest[tnames,4],fmeth$rcvFilter_UniCox$survMedianTest[tnames,5],fmeth$rcvFilter_UniCox$survMedianTest[tnames,6]);
      filters = c(filters,"COX.UnivariateCox");
    }
      
    
    columnNamesMethods <- NULL;
    for(x in filters)
    {
      for(y in predictions)
      {
        columnNamesMethods <- cbind(columnNamesMethods,paste(x,y,sep=""))
      }
    }  
    
    colnames(test_Predictions) <-	c(classnames,columnNamesMethods);
    selFrequency <- cbind(selFrequency,numeric(ncol(theData)));
    if (!inherits(fmeth$rcvFilter_UniCox, "try-error"))
    {
      selFrequency[names(fmeth$rcvFilter_UniCox$featureFrequency),ncol(selFrequency)] <- fmeth$rcvFilter_UniCox$featureFrequency;
      theFiltersets <- c(theFiltersets,"Cox.Unicox");
      jaccard_filter$kendall <- fmeth$rcvFilter_UniCox$jaccard;    
    }
    
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
  for (i in 2:ncol(test_Predictions))
  {
    if (test_Predictions[,i] < -1)
    {	
      test_Predictions[,i] <- 1.0/(1.0+exp(-test_Predictions[,i] ));
    }
  }
  
  
  result <- list(errorciTable = errorciTable,accciTable = accciTable,aucTable = aucTable,senTable = senTable,speTable = speTable,
                 errorciTable_filter = errorciTable_filter,accciTable_filter = accciTable_filter,aucTable_filter = aucTable_filter,senTable_filter = senTable_filter,speTable_filter = speTable_filter,
                 CIRisksTable = CIRisksTable,CIFollowUPTable = CIFollowUPTable,LogRankTable = LogRankTable,
                 CIRisksTable_filter = CIRisksTable_filter,CIFollowUPTable_filter = CIFollowUPTable_filter,LogRankTable_filter = LogRankTable_filter,
                 times = list(Reference = referenceCV$theTimes,LASSO = rcvLASSO$theTimes, RIDGE = rcvGLMNET_RIDGE$theTimes, ELASTICNET = rcvGLMNET_ELASTICNET, BESS = rcvBESS$theTimes, BESS.SEQUENTIAL = rcvBESSSequential$theTimes, BESS.SEQUENTIAL.BIC=rcvBESSSequentialBIC$theTimes),
                 jaccard = jaccard,
                 featsize = featsize,
                 TheCVEvaluations = TheCVEvaluations,
                 thesets = thesets,
                 theMethod = theMethod,
                 theFiltersets = theFiltersets,
                 testPredictions = test_Predictions,
                 featureSelectionFrequency = selFrequency,
                 cpuElapsedTimes=cputimes				 
  )
  class(result) <- c("FRESA_benchmark","Survival.COX");
  return(result)
}

