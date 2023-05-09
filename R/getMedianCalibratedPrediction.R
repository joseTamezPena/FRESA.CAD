getMedianSurvCalibratedPrediction <- function(testPredictions)
{
  medianTestCalibrated <- NULL
  allPred <- testPredictions
  ids <- rownames(allPred)
  cnames <-  paste("ID_",rownames(allPred),allPred[,3],sep="_")
  rownames(allPred) <-cnames
  allPred <- as.data.frame(allPred)
  allPred <- as.data.frame(apply(allPred,2,unlist))
  colnames(allPred) <- c("Times","Outcome","Model","LinearPredictors","CalLinearPredictors")
  
  totmodels <- table(allPred$Model)
  posOutcome <- allPred$Outcome==1
  negOutcome <- allPred$Outcome==0
  for (mn in as.numeric(names(totmodels)))
  {
    modelRows <- allPred$Model==mn
    mpredPos <- median(allPred[modelRows & posOutcome,4])
    mpredNeg <- median(allPred[modelRows & negOutcome,4])
    allPred[modelRows,5] <- allPred[modelRows,4] - (mpredPos + mpredNeg)/2
  }
  
  medtest <- boxplot(allPred$LinearPredictors~ids,plot = FALSE);
  rids <- medtest$names
  medtest <- medtest$stats[3,]
  calmedtest <- boxplot(allPred$CalLinearPredictors~ids,plot = FALSE);
  calmedtest <- calmedtest$stats[3,]
  outtest <- boxplot(allPred$Outcome~ids,plot = FALSE);
  outtest <- outtest$stats[3,]
  timestest <- boxplot(allPred$Times~ids,plot = FALSE);
  timestest <- timestest$stats[3,]
  medianTestCalibrated <- cbind(timestest,outtest,medtest,calmedtest)
  medianTestCalibrated <- as.data.frame(medianTestCalibrated)
  colnames(medianTestCalibrated) <- c("Times","Outcome","LinearPredictorsMedian","CalLinearPredictorsMedian")
  rownames(medianTestCalibrated) <- rids
  return (medianTestCalibrated)
}

getMedianLogisticCalibratedPrediction <- function(testPredictions)
{
  medianTestCalibrated <- NULL
  allPred <- testPredictions
  ids <- rownames(allPred)
  cnames <-  paste("ID_",rownames(allPred),allPred[,2],sep="_")
  rownames(allPred) <-cnames
  allPred <- as.data.frame(allPred)
  allPred <- as.data.frame(apply(allPred,2,unlist))
  allPred$CalLinearPredictors <- allPred[,3]
  colnames(allPred) <- c("Outcome","Model","Predictors","CalPredictors")
  minPred <- min(allPred$Predictors)
  maxPred <- max(allPred$Predictors)
  isProbability <- (minPred >= 0) && (maxPred<=1.0)
  if (isProbability) # They are probabilities. Convert to linear predictions
  {
    allPred$Predictors <- log(allPred$Predictors/(1.0-allPred$Predictors))
  }
  
  totmodels <- table(allPred$Model)
  posOutcome <- allPred$Outcome==1
  negOutcome <- allPred$Outcome==0
  for (mn in as.numeric(names(totmodels)))
  {
    modelRows <- allPred$Model==mn
    mpredPos <- median(allPred[modelRows & posOutcome,3])
    mpredNeg <- median(allPred[modelRows & negOutcome,3])
    allPred[modelRows,4] <- allPred[modelRows,3] - (mpredPos + mpredNeg)/2
  }
  
  medtest <- boxplot(allPred$Predictors~ids,plot = FALSE);
  rids <- medtest$names
  medtest <- medtest$stats[3,]
  calmedtest <- boxplot(allPred$CalLinearPredictors~ids,plot = FALSE);
  calmedtest <- calmedtest$stats[3,]
  outtest <- boxplot(allPred$Outcome~ids,plot = FALSE);
  outtest <- outtest$stats[3,]
  medianTestCalibrated <- cbind(outtest,medtest,calmedtest)
  medianTestCalibrated <- as.data.frame(medianTestCalibrated)
  colnames(medianTestCalibrated) <- c("Outcome","PredictorsMedian","CalPredictorsMedian")
  rownames(medianTestCalibrated) <- rids
  
  if (isProbability) # we estimated linear. Convert back to probabilites
  {
    allPred$PredictorsMedian <- 1.0/(1.0+exp(-allPred$PredictorsMedian))
    allPred$CalPredictorsMedian <- 1.0/(1.0+exp(-allPred$CalPredictorsMedian))
  }

  
  return (medianTestCalibrated)
  
}
