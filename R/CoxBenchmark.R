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
      rcvFilter_reference <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = referenceCV$selectedFeaturesSet);
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
      
      rcvFilter_LASSO <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvLASSO$selectedFeaturesSet);
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
      
      
      rcvFilter_BESS <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = rcvBESS$selectedFeaturesSet);
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
      
      cat("Univariate cox Feature Selection: ");
      rcvFilter_UniCox <- randomCV(theData,theOutcome,clasfun,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = univariate_cox);
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
    if (is.null(referenceCV))
    {
      cat("Modeling BSWiMS: + Model found, - No Model \n"); 
      referenceCV <- randomCV(theData,theOutcome,BSWiMS.model,trainFraction = trainFraction,repetitions = reps,featureSelectionFunction = "Self");
      referenceName = "BSWiMS";
      referenceFilterName = "COX.BSWiMS";
    }
  
    cStats <- predictionStats_survival(referenceCV$survMedianTest,plotname = referenceName);
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    #referenceCV <- cpFinal$TheCVEvaluations$Reference
    binaryPreds <- referenceCV$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"BSWiMS")
    accciTable <- rbind(accciTable,binaryStats$accc)
    errorciTable <- rbind(errorciTable,binaryStats$berror)
    aucTable <- rbind(aucTable,binaryStats$aucs)
    senTable <- rbind(senTable,binaryStats$sensitivity)
    speTable <- rbind(speTable,binaryStats$specificity)
    
    # 1 - pchisq(cStats$LogRank$chisq, length(cStats$LogRank$n) - 1)
    ######################LASSO#################################### 
    rcvLASSO <- randomCV(theData,theOutcome,LASSO_MIN,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self");
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
    
    ######################BESS#################################### 
    rcvBESS <- randomCV(theData,theOutcome,BESS,trainSampleSets = referenceCV$trainSamplesSets,featureSelectionFunction = "Self");
    cStats <- predictionStats_survival(rcvBESS$survMedianTest,plotname = "BeSS", noNegativeRisk = TRUE);
    CIFollowUPTable <- rbind(CIFollowUPTable,cStats$CIFollowUp);
    CIRisksTable <- rbind(CIRisksTable,cStats$CIRisk);
    LogRankTable <- rbind(LogRankTable,cStats$LogRank);
    
    #rcvBESS <- cpFinal$TheCVEvaluations$BESS
    #rcvBESS$survMedianTest[,"LinearPredictorsMedian"] <- -rcvBESS$survMedianTest[,"LinearPredictorsMedian"]
    binaryPreds <- rcvBESS$survMedianTest[,c("Outcome","LinearPredictorsMedian")]
    binaryStats <- predictionStats_binary(binaryPreds,"BeSS")
    accciTable <- rbind(accciTable,binaryStats$accc)
    errorciTable <- rbind(errorciTable,binaryStats$berror)
    aucTable <- rbind(aucTable,binaryStats$aucs)
    senTable <- rbind(senTable,binaryStats$sensitivity)
    speTable <- rbind(speTable,binaryStats$specificity)
    
    ######################Esemble#################################### 
    
    # ens <- cbind(referenceCV$survMedianTest[,1],referenceCV$survMedianTest[,2],rowMeans(cbind(1.0*(referenceCV$survMedianTest[,2] > 0),1.0*(rcvLASSO$survMedianTest[,2] > 0),rcvRF$survMedianTest[,2],rcvKNN$survMedianTest[,2],rcvSVM$survMedianTest[,2])))
    # rowMeans(cbind(1.0*(referenceCV$survMedianTest[,2] > 0),1.0*(rcvLASSO$survMedianTest[,2] > 0),rcvRF$survMedianTest[,2],rcvKNN$survMedianTest[,2],rcvSVM$survMedianTest[,2]))
    # cStats <- predictionStats_binary(ens,plotname = "Ensemble",center = TRUE,cex=0.8);
      
    ######################Filters  #################################### 
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
    
    
    ######################Predictions union  #################################### 
    test_Predictions <- referenceCV$survMedianTest;
    tnames <- rownames(test_Predictions)
    test_Predictions <- cbind(test_Predictions,rcvLASSO$survMedianTest[tnames,3],rcvLASSO$survMedianTest[tnames,4],rcvLASSO$survMedianTest[tnames,5],rcvLASSO$survMedianTest[tnames,6])
    test_Predictions <- cbind(test_Predictions,rcvBESS$survMedianTest[tnames,3],rcvBESS$survMedianTest[tnames,4],rcvBESS$survMedianTest[tnames,5],rcvBESS$survMedianTest[tnames,6])
    test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_reference$survMedianTest[tnames,3],fmeth$rcvFilter_reference$survMedianTest[tnames,4],fmeth$rcvFilter_reference$survMedianTest[tnames,5],fmeth$rcvFilter_reference$survMedianTest[tnames,6]);
    test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_LASSO$survMedianTest[tnames,3],fmeth$rcvFilter_LASSO$survMedianTest[tnames,4],fmeth$rcvFilter_LASSO$survMedianTest[tnames,5],fmeth$rcvFilter_LASSO$survMedianTest[tnames,6]);
    test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_BESS$survMedianTest[tnames,3],fmeth$rcvFilter_BESS$survMedianTest[tnames,4],fmeth$rcvFilter_BESS$survMedianTest[tnames,5],fmeth$rcvFilter_BESS$survMedianTest[tnames,6]);
    test_Predictions <- cbind(test_Predictions,fmeth$rcvFilter_UniCox$survMedianTest[tnames,3],fmeth$rcvFilter_UniCox$survMedianTest[tnames,4],fmeth$rcvFilter_UniCox$survMedianTest[tnames,5],fmeth$rcvFilter_UniCox$survMedianTest[tnames,6]);

    ######################Column names  #################################### 
    predictions <- c("MartinGale","LinearPredictors","FollowUpTimes","Risks");
    methods <- c(referenceName,"LASSO","BESS");
    
    filters <- c("COX.BSWiMS","COX.LASSO","COX.BESS","COX.UnivariateCox");
    # filters <- c(paste("COX.",referenceFilterName,sep=""),"COX.LASSO","COX.BESS","COX.IDI","COX.NRI","COX.tStudent","COX.Wilcoxon","COX.Kendall","COX.mRMR","UnivariateCox")
    columnNamesMethods <- NULL;
    
    for(x in methods)
    {
      for(y in predictions)
      {
        columnNamesMethods <- cbind(columnNamesMethods,paste(x,y,sep=""))
      }
    }
    
    for(x in filters)
    {
      for(y in predictions)
      {
        columnNamesMethods <- cbind(columnNamesMethods,paste(x,y,sep=""))
      }
    }
    
    methodsAndFilters <-c(methods,filters);
    columnNamesMethods <- cbind("Times","Outcome",columnNamesMethods)
    colnames(test_Predictions) <- columnNamesMethods;
    test_Predictions <- as.data.frame(test_Predictions)
    
    thesets <- c("Survival Algorithm")
    
    rownames(CIFollowUPTable) <- methods;
    rownames(CIRisksTable) <- methods;
    rownames(LogRankTable) <- methods;

    rownames(CIFollowUPTable_filter) <- filters;
    rownames(CIRisksTable_filter) <- filters;
    rownames(LogRankTable_filter) <- filters;
    
    rownames(accciTable) <- methods;
    rownames(errorciTable) <- methods;
    rownames(aucTable) <- methods;
    rownames(senTable) <- methods;
    rownames(speTable) <- methods;
    
    rownames(accciTable_filter) <- filters;
    rownames(errorciTable_filter) <- filters;
    rownames(aucTable_filter) <- filters;
    rownames(senTable_filter) <- filters;
    rownames(speTable_filter) <- filters;
    
    ff <- names(referenceCV$featureFrequency)
    ff <- c(ff,names(rcvLASSO$featureFrequency))
    ff <- c(ff,names(rcvBESS$featureFrequency))
    ff <- c(ff,names(fmeth$rcvFilter_reference$featureFrequency))
    ff <- c(ff,names(fmeth$rcvFilter_LASSO$featureFrequency))
    ff <- c(ff,names(fmeth$rcvFilter_BESS$featureFrequency))
    ff <- c(ff,names(fmeth$rcvFilter_UniCox$featureFrequency))
    ff <- unique(ff)
    
    Nvar <- min(c(1000,length(ff)))
    selFrequency <- matrix(0,nrow = Nvar,ncol = length(methodsAndFilters))
    rownames(selFrequency) <- ff[1:Nvar]
    selnames <- rownames(selFrequency)
    colnames(selFrequency) <- methodsAndFilters
    ff <- referenceCV$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,referenceName] <- ff[selnames[fnames]]
    ff <- rcvLASSO$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,"LASSO"] <- ff[selnames[fnames]]
    ff <- rcvBESS$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,"BESS"] <- ff[selnames[fnames]]
    ff <- fmeth$rcvFilter_reference$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,referenceFilterName] <- ff[selnames[fnames]]
    ff <- fmeth$rcvFilter_LASSO$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,"COX.LASSO"] <- ff[selnames[fnames]]
    ff <- fmeth$rcvFilter_BESS$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,"COX.BESS"] <- ff[selnames[fnames]]
    ff <- fmeth$rcvFilter_UniCox$featureFrequency
    fnames <- selnames %in% names(ff)
    selFrequency[fnames,"COX.UnivariateCox"] <- ff[selnames[fnames]]


    elapcol <- names(referenceCV$theTimes) == "elapsed"
    theMethod <- c(referenceName,"LASSO","BESS")
    cputimes <- list(Reference = mean(referenceCV$theTimes[ elapcol ]),LASSO = mean(rcvLASSO$theTimes[ elapcol ]),BESS = mean(rcvBESS$theTimes[ elapcol ]))
    
    
    jaccard_filter = list(Reference = referenceCV$jaccard,
                          LASSO = rcvLASSO$jaccard,
                          BESS = rcvBESS$jaccard,
                          COX.Reference = fmeth$rcvFilter_reference$jaccard,
                          COX.LASSO = fmeth$rcvFilter_LASSO$jaccard,
                          COX.BESS = fmeth$rcvFilter_BESS$jaccard,
                          COX.UniCox = fmeth$rcvFilter_UniCox$jaccard
    );
    
    featsize <- unlist(lapply(jaccard_filter, `[`, c('averageLength')))
    names(featsize) <- methodsAndFilters;
    jaccard <- unlist(lapply(jaccard_filter, `[`, c('Jaccard.SM')))
    names(jaccard) <- methodsAndFilters;
    cputimes <- unlist(cputimes);
    cputimes <- c(cputimes,sum(cputimes));
    names(cputimes) <- methods;
    
    result <- list(errorciTable = errorciTable,accciTable = accciTable,aucTable = aucTable,senTable = senTable,speTable = speTable,
                   errorciTable_filter = errorciTable_filter,accciTable_filter = accciTable_filter,aucTable_filter = aucTable_filter,senTable_filter = senTable_filter,speTable_filter = speTable_filter,
                   CIRisksTable = CIRisksTable,CIFollowUPTable = CIFollowUPTable,LogRankTable = LogRankTable,
                   CIRisksTable_filter = CIRisksTable_filter,CIFollowUPTable_filter = CIFollowUPTable_filter,LogRankTable_filter = LogRankTable_filter,
                   times = list(Reference = referenceCV$theTimes,LASSO = rcvLASSO$theTimes, BESS = rcvBESS$theTimes),
                   jaccard = jaccard,
                   featsize = featsize,
                   TheCVEvaluations = list(Reference = referenceCV,
                                           LASSO = rcvLASSO,
                                           BESS = rcvBESS,
                                           COX.Reference = fmeth$rcvFilter_reference,
                                           COX.LASSO = fmeth$rcvFilter_LASSO,
                                           COX.BESS = fmeth$rcvFilter_BESS,
                                           Cox.UniCox = fmeth$rcvFilter_UniCox
                   ),
                   thesets = thesets,
                   theMethod = methods,
                   theFiltersets = filters,
                   testPredictions = test_Predictions,
                   featureSelectionFrequency = selFrequency,
                   cpuElapsedTimes=cputimes				 
    )
    class(result) <- c("FRESA_benchmark","Survival.COX");
    return(result)
}

