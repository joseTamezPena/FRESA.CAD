randomCV <-  function(theData = NULL, 
                      theOutcome = "Class",
                      fittingFunction=NULL, 
                      trainFraction = 0.5, 
                      repetitions = 100,
                      trainSampleSets=NULL,
                      featureSelectionFunction=NULL,
                      featureSelection.control=NULL,
                      asFactor=FALSE,
                      addNoise=FALSE,
                      classSamplingType=c("Proportional","Balanced","Augmented","LOO"),
                      testingSet=NULL,...)
{
  formula.list <- NULL
  classSamplingType <- match.arg(classSamplingType);
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
  
  fullSet <- theData
  if (!is.null(testingSet))
  {
    fullSet <- rbind(fullSet,testingSet);
  }
  
  survpredict <- function(currentModel,Dataset,TestDataset,selectedFeatures)
  {
    theSurvData <-Dataset;
    if (!is.null(currentModel$fit))
    {
      predi_t <- predict(currentModel,TestDataset);
      tpred <- attr(predi_t,"dataAdjusted")
      if (!is.null(tpred))
      {
        TestDataset <- tpred
      }
      currentModel <- currentModel$fit;
    }
    fclass <- class(currentModel)
    if (length(fclass)>1) fclass <- fclass[1];
    
    
    if(length(selectedFeatures)>=nrow(theSurvData))
    {
      selectedFeatures <- head(selectedFeatures, nrow(theSurvData)-1)
      warning("The number of selected features is bigger than the number of observations, the top will be used.")
    }

    numberCoeficients <- length(selectedFeatures)
    
    if(length(selectedFeatures)==0)
    {
      warning("Method did not select any features")
      return(list(
                  linearPredictors=rep(NA,nrow(TestDataset)),
                  risks = list(fit=rep(NA,nrow(TestDataset)),se.fit=rep(NA,nrow(TestDataset))),
                  hr = rep(NA,nrow(TestDataset))));
    }
    
    baseformula <- as.character(theformula);
    formulaCox <- as.formula(paste(paste(baseformula[2],"~"), paste(selectedFeatures, collapse='+')));
    cox <- try(survival::coxph(formula=formulaCox, data=theSurvData));

    if (fclass == "FRESA_GLMNET")
    {
      cox$coefficients <- currentModel$coef[1:numberCoeficients];
    }
    if (fclass == "fitFRESA")
    {
        cox$coefficients <- currentModel$bagging$bagged.model$coefficients[-c(1)][1:numberCoeficients];
    }
    if(fclass=="FRESA_BESS")
    {
      names(currentModel$fit$bestmodel$coefficients) <- selectedFeatures;
      cox$coefficients <- currentModel$fit$bestmodel$coefficients[1:numberCoeficients];
    }
    
    if (!inherits(cox,"try-error") && !inherits(cox,"list") )
    {
      #risk
      risks <- try(predict(cox,newdata=TestDataset,type="risk"))
      if (!inherits(risks,"try-error"))
      {
#        
#
        #linear predictos
        linearPredictors <- predict(cox,newdata=TestDataset,type="lp")
        
        hr <- round(coef(summary(cox))[,2],3)
        survPreds <- list(linearPredictors=linearPredictors,
                          risks = risks,
                          hr = hr)
        return (survPreds)
      }
      else{
        warning("Cox Fit Follow-up Times Error");
        return(list(
                    linearPredictors=rep(0,nrow(TestDataset)),
                    risks = list(fit=rep(0,nrow(TestDataset)),se.fit=rep(0,nrow(TestDataset))),
                    hr = rep(0,nrow(TestDataset))));
      }
    }
    else{
      warning("Cox Fit Error");
        return(list(
                    linearPredictors=rep(0,nrow(TestDataset)),
                    risks = list(fit=rep(0,nrow(TestDataset)),se.fit=rep(0,nrow(TestDataset))),
                    hr = rep(0,nrow(TestDataset))));
    }
  }
  

  jaccard <-  function(featureSet)
  {
    Jaccard.SM <- 0;
    averageLength <- 0;
    tota <- 0;
    loops <- length(featureSet);
    for (n in 1:loops)
    {
      feat <- featureSet[[n]];
      #		print(feat);
      lft <- length(feat);
      averageLength <- averageLength+lft
      if (lft>0)
      {
        if (n<loops)
        {
          for (i in (n+1):loops)
          {
            feat2 <- featureSet[[i]];
            if (length(feat2) > 0)
            {
              Jaccard.SM = Jaccard.SM+sum(duplicated(c(feat2,feat)))/length(unique(c(feat2,feat)));
              tota = tota + 1;
            }
          }
        }
      }
    }
    averageLength <- averageLength/loops;
    if (tota>1) Jaccard.SM = Jaccard.SM/tota;
    result <- list(Jaccard.SM=Jaccard.SM,averageLength=averageLength);
    return(result);
  }
  
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    install.packages("glmnet", dependencies = TRUE)
  } 
  if (!requireNamespace("BeSS", quietly = TRUE)) {
    install.packages("BeSS", dependencies = TRUE)
  } 
  
  theformula <- NULL;
  theTime <- NULL;
  varsmod <- NULL;
  isSurv <- FALSE;
  if (inherits(theOutcome,"formula"))
  {
    theformula <- theOutcome;
    varsmod <- all.vars(theformula);
    theOutcome <- varsmod[1];
    varsmod <- varsmod[varsmod!="."]
    baseformula <- as.character(theformula);
    if (sum(str_count(baseformula,"Surv")) > 0)
    {
      theOutcomeSurv <- theformula;
      theOutcome <- varsmod[2];
      theTime <- varsmod[1];
      isSurv <- TRUE;
    }
  }
  else
  {
    varsmod <- theOutcome;
    theformula <- formula(paste(theOutcome,"~ ."));
  }
  
  dataTable <- table(theData[,theOutcome]);
  theClasses <- as.numeric(names(dataTable));
  classLen <- length(theClasses);
  
  selectedFeaturesSet <- list();
  testClases <- ((classLen < 10) && (min(dataTable) >= 5));
  BootReplace = FALSE;
  if (inherits(trainFraction,"character"))
  {
	 BootReplace <- (trainFraction == "Bootstrap");
	 trainFraction <- 0.5 + 0.5*BootReplace;
  }
  samplePerClass <- as.integer((nrow(theData)/classLen)*trainFraction);
  if (testClases)
  {
    ssubsets <- list();
    samplePerClass <- as.integer((nrow(theData)/classLen)*trainFraction);
    if (classSamplingType == "Balanced")
    {
      samplePerClass <- as.integer(min(dataTable)*trainFraction);
    }
    if (classSamplingType == "Augmented")
    {
      samplePerClass <- as.integer(max(dataTable)*trainFraction);
    }
    jind <- 1;
    for (s in theClasses)
    {
      ssubsets[[jind]] <- subset(theData,get(theOutcome) == s);
      jind <- jind + 1;
    }
  }
  trainSamplesSets <- list();
  tset <- 1;
  testPredictions <- NULL;
  trainPredictions <- NULL;
  theTimes <- NULL;
  survTestPredictions <- NULL;
  survTrainPredictions <- NULL;
  avgsel <- 0;

  theVars <-colnames(theData)[!colnames(theData) %in% c(as.character(theTime),as.character(theOutcome))];
  survHR <- data.frame(matrix(ncol = length(theVars), nrow = 0));
  colnames(survHR) <- theVars;
  
  if (!is.null(trainSampleSets))
  {
    repetitions <- trainSampleSets$repetitions
  }
  MADERROR <- numeric(repetitions);
  for (rept in 1:repetitions)
  {
    #		cat(ncol(theData),":",nrow(theData),"\n")
    cat(".");
    nfet <- TRUE;
    selectedFeaturesSet[[rept]] <- character();
    #		cat(length(selectedFeaturesSet),"\n");
    if (testClases)
    {
      jind <- 1;
      trainSet <- NULL;
      testSet <- NULL;
      for (s in theClasses)
      {
        sampleTrain <- NULL;
        if (is.null(trainSampleSets))
        {
          if (classSamplingType == "LOO")
          {
            smp <- seq(1,nrow(ssubsets[[jind]]));
            sampleTrain <- smp[-(1 + (rept %% nrow(ssubsets[[jind]])))];
          }
          else
          {
            minfrac <- min(c(trainFraction,0.99));
            if (classSamplingType == "Proportional")
            {
              sampleTrain <- sample(nrow(ssubsets[[jind]]),as.integer(nrow(ssubsets[[jind]])*minfrac),replace=BootReplace);
            }
            else
            {
              ssize <- nrow(ssubsets[[jind]])*minfrac;
              if (samplePerClass > ssize)
              {
                sampleTrain <- sample(nrow(ssubsets[[jind]]),ssize,replace=BootReplace);
                if (classSamplingType == "Augmented")
                {
                  therest <- samplePerClass-ssize;
                  nsample <- sample(ssize,therest,replace=TRUE);
                  sampleTrain <- append(sampleTrain,sampleTrain[nsample]);
                }
              }
              else
              {
                sampleTrain <- sample(nrow(ssubsets[[jind]]),samplePerClass,replace=BootReplace);
              }
            }
          }
		  usample <- unique(sampleTrain);
		  if (length(usample) == nrow(ssubsets[[jind]]))
		  {	
			 sampleTrain <- sampleTrain[sampleTrain != usample[1]]; # remove at least one sample
		  }
          trainSamplesSets[[tset]] <- sampleTrain;
        }
        else
        {
          sampleTrain <- trainSampleSets[[tset]];
        }
        tset <- tset + 1;
        
        samdup <- duplicated(sampleTrain)
        if (sum(samdup) > 0 )
        {
          noPerturbedTrain <- ssubsets[[jind]][sampleTrain[!samdup],];
          tobePerturbed <- as.data.frame(ssubsets[[jind]][sampleTrain[samdup],])
          rownames(tobePerturbed) <- paste(rownames(tobePerturbed),"D",sep="");
          fnames <- !(colnames(ssubsets[[jind]]) %in% c(varsmod));
          cols <- sum(fnames);
          if (cols > 1)
          {
            fnames <- colnames(ssubsets[[jind]])[fnames];
            for (nf in fnames)
            {
              dto <- tobePerturbed[,nf]
              tb <- table(noPerturbedTrain[,nf])
              
              if (length(tb) > 5)
              {
#                if (nf == fnames[3])
#                {
#                  hist(noPerturbedTrain[,nf])
#                  hist(tobePerturbed[,nf])
#               }
                tbnames <- names(tb);
                indx <- 1:length(tb);
                names(indx) <- tbnames;
                ralea <- runif(length(dto),-1.0,1.0);
                peturblen <- min(1,0.025*length(tb));
                nidx <- indx[as.character(dto)] + floor(peturblen*ralea + 0.5);
                nidx[nidx < 1] <- 1;
                nidx[nidx > length(tb)] <- length(tb);
                tobePerturbed[,nf] <- 0.5*(dto + as.numeric(names(tb[nidx])));
#                if (nf == fnames[3])
#                {
#                  hist(tobePerturbed[,nf])
#                }
              }
#              else
#              {
#                tobechanged <- (runif(length(dto)) < 0.1);
#                tobePerturbed[tobechanged,nf] <- noPerturbedTrain[sample(nrow(noPerturbedTrain),sum(tobechanged)),nf];
#              }
            }
          }
          trainSet <- rbind(trainSet,noPerturbedTrain,tobePerturbed);
        }
        else
        {
          trainSet <- rbind(trainSet,ssubsets[[jind]][sampleTrain,]);
        }
        
        testSet <- rbind(testSet,ssubsets[[jind]][-unique(sampleTrain),]);
        
        jind <- jind + 1;
#        cat("Train Class: ",s," rows:",nrow(trainSet),"\n");
#        cat("Test Class: ",s," rows:",nrow(testSet),"\n");
      }
    }
    else
    {
      if (is.null(trainSampleSets))
      {
        sampleTrain <- sample(nrow(theData),nrow(theData)*trainFraction,replace=BootReplace)
		usample <- unique(sampleTrain);
		if (length(usample) == nrow(theData))
		{	
			 sampleTrain <- sampleTrain[sampleTrain != usample[1]]; # remove at least one sample
		}
        trainSamplesSets[[tset]] <- sampleTrain;
      }
      else
      {
        sampleTrain <- trainSampleSets[[tset]];
      }
      tset <- tset + 1;
      trainSet <- theData[sampleTrain,];
      testSet <- theData[-unique(sampleTrain),];
    }
    
    if (!is.null(testingSet))
    {
       toremove <- (rownames(testingSet) %in% rownames(trainSet)) | (rownames(testingSet) %in% rownames(testSet))
       testSet <- rbind(testSet,testingSet[!toremove,]);
    }
    
    selnames <- character();
    if (inherits(featureSelectionFunction,"list"))
    {
      if (!is.null(featureSelectionFunction[[rept]]))
      {
        trainSet <- trainSet[,c(varsmod,featureSelectionFunction[[rept]])];
        nfet <- (length(featureSelectionFunction[[rept]]) > 0)
        selnames <- featureSelectionFunction[[rept]];
        if (!is.null(selnames))
        {
          selectedFeaturesSet[[rept]] <- selnames;
        }
        #				cat("List: ",length(selectedFeaturesSet),":",ncol(trainSet),"\n");
      }
      else
      {
        nfet <- FALSE;
        selnames <- character();	
      }
    }
    else
    {
      if (inherits(featureSelectionFunction,"function"))
      {
        #				print(tracemem(trainSet))
        if (!is.null(featureSelection.control))
        {
          frank <- do.call(featureSelectionFunction,c(list(trainSet,theOutcome),featureSelection.control));
          if (length(frank)>0)
          {
            selectedFeaturesSet[[rept]] <- names(frank);
          }
        }
        else
        {
          if(isSurv)
          {
            frank <- featureSelectionFunction(trainSet,theformula)
          }
          else
          {
            frank <- featureSelectionFunction(trainSet,theOutcome)
          }
          
          if (length(frank)>0)
          {
            selectedFeaturesSet[[rept]] <- names(frank);
          }
        }
        
        if(isSurv){
          selectedFeaturesSet[[rept]] <- selectedFeaturesSet[[rept]][!selectedFeaturesSet[[rept]] %in% theTime]
        }
        
        nfet <- (length(selectedFeaturesSet[[rept]]) > 0)
        if (nfet) trainSet <- trainSet[,c(varsmod,selectedFeaturesSet[[rept]])];
        selnames <- selectedFeaturesSet[[rept]];
        #				cat("Function:",length(selectedFeaturesSet),":",ncol(trainSet),"\n");
      }
    }
    #		cat(length(selectedFeaturesSet),":",ncol(trainSet),"\n");
    if (nfet)
    {
      #			print(selectedFeaturesSet[[rept]]);
      if (addNoise)
      {
        fnames <- !(colnames(trainSet) %in% c(varsmod));
        cols <- sum(fnames);
        if (cols > 1)
        {
          fnames <- colnames(trainSet)[fnames];
          nlevel = 1.0*addNoise;
          nlthr <- min(nlevel/3,0.5);
          sdg <- apply(trainSet[,fnames],2,sd,na.rm = TRUE);
          iqrg <- apply(trainSet[,fnames],2,IQR,na.rm = TRUE)/(2*abs(qnorm(0.25)));
          iqrg[iqrg == 0] <- sdg[iqrg == 0];
          iqrg[iqrg == 0] <- 1.0e-10;
          rows <- nrow(trainSet);
          for (nf in fnames)
          {
            dto <- trainSet[,nf];
            tb <- table(dto)
            if (length(tb) > 5)
            {
              noise <- as.numeric(rnorm(rows,0,iqrg[nf]/length(tb)));

              tbnames <- names(tb);
              indx <- 1:length(tb);
              names(indx) <- tbnames;
              ralea <- runif(length(dto),-1.0,1.0);
              peturblen <- min(1,nlthr*0.05*length(tb))

              nidx <- indx[as.character(dto)] + floor(peturblen*ralea + 0.5);
              nidx[nidx < 1] <- 1;
              nidx[nidx > length(tb)] <- length(tb);
              dto <- 0.5*(dto + as.numeric(names(tb[nidx])));

              trainSet[,nf] <- dto + nlevel*noise;
            }
#            else
#            {
#              tobechanged <- (runif(length(dto)) < nlthr/5.0);
#              trainSet[tobechanged,nf] <- dto[sample(length(dto),sum(tobechanged))];
#            }
          }
         }
      }
      if ((testClases) && (asFactor))
      {
        trainSet[,theOutcome] <- as.factor(trainSet[,theOutcome]); 
        #					print(selnames)
      }
      #			cat(ncol(trainSet),":",nrow(trainSet),"\n")
      
      theTimes <- append(theTimes,system.time(currentModel <- try(fittingFunction(theformula,trainSet,...))));
#	  cat("Here\n")
      if ( inherits(currentModel, "try-error"))
      {
        if ((testClases) && (!asFactor))
        {
          warning("Fit Error. I will try outcome as Factor");
          trainSet[,theOutcome] <- as.factor(trainSet[,theOutcome]); 
          currentModel <- try(fittingFunction(theformula,trainSet,...));
        }
        else
        {
          olength <- length(selnames);
          if (olength>2) 
          {
            selnames <- correlated_Remove(trainSet,selnames);
            attr(selnames,"CorrMatrix") <- NULL;

            if ((length(selnames)>1) && (olength>length(selnames)))
            {
              warning("Fit Error. I will try with less features. Original Number of features:",olength," New number of features:",length(selnames),"\n");
              trainSet <- trainSet[,c(theOutcome,selnames)];
              
              currentModel <- try(fittingFunction(theformula,trainSet,...));
            }
          }
        }
        if ( inherits(currentModel, "try-error"))
        {
          cat("Fit Error. Number of features:",length(selnames),"\n");
        }
      }
      if ( !inherits(currentModel, "try-error"))
      {
        #				if ((class(featureSelectionFunction) == "character") || (class(featureSelectionFunction) == "function"))
        {
          fclass <- class(currentModel)
          if (length(fclass)>1) fclass <- fclass[1];
          #				print(fclass);
          if (fclass == "fitFRESA")
          {
            ffet <- currentModel$bagging$frequencyTable;
            if (!is.null(ffet))
            {
              selectedFeaturesSet[[rept]] <- names(ffet);
            }
          }
          if (!is.null(currentModel$selectedfeatures))
          {
            ffet <- currentModel$selectedfeatures;
#			print(ffet);
            selectedFeaturesSet[[rept]] <- ffet;
          }
          if (fclass == "FRESA_BESS")
          {
            if (length(selectedFeaturesSet[[rept]])>0){
              if (head(selectedFeaturesSet[[rept]],1)==""){
                selectedFeaturesSet[[rept]] <- list()
              }
            }
          }
          
          vs <- NULL;
          if (!is.null(currentModel$importance))
          {
            vs <- currentModel$importance[,1];
          }
          if (!is.null(currentModel$variable.importance))
          {
            vs <- currentModel$variable.importance;
          }
          if (!is.null(vs))
          {
            #						print(vs)
            if (length(vs)>1)
            {
              if (length(vs) > 2) vs <- vs[order(-vs)];
              if (sum(vs > 0.0) > 1)
              {
                vs <- vs[vs > 0.0];							
              }
              selectedFeaturesSet[[rept]] <- names(vs);
            }
          }
        }
        if (!is.null(currentModel$formula.list))
        {
          formula.list <- c(formula.list,currentModel$formula.list)
        }
        else
        {
          if (!is.null(currentModel$fit))
          {
            if (!is.null(currentModel$fit$formula.list))
            {
              formula.list <- c(formula.list,currentModel$fit$formula.list)
            }
          }
        }
        #				print(selectedFeaturesSet);
        if ((length(selectedFeaturesSet[[rept]])>0) || is.null(featureSelectionFunction))
        {
          pred <- rpredict(currentModel,testSet,asFactor,classLen,...);
          ctestPredictions <- cbind(testSet[,theOutcome],rep(rept,nrow(testSet)),pred);
          pred <- rpredict(currentModel,trainSet,asFactor,classLen,...);
          ctrainPredictions <- cbind(trainSet[,theOutcome],rep(rept,nrow(trainSet)),pred);
          rownames(ctestPredictions) <- rownames(testSet);
          rownames(ctrainPredictions) <- rownames(trainSet);
          testPredictions <- rbind(testPredictions,ctestPredictions);
          trainPredictions <- rbind(trainPredictions,ctrainPredictions);
          if (isSurv)
          {
            #SURVPREDICT
            survPreds <- survpredict(currentModel,trainSet,testSet,selectedFeaturesSet[[rept]]);
            csurvTestPredictions <- cbind(testSet[,theTime],testSet[,theOutcome],rep(rept,nrow(testSet)),survPreds$linearPredictors,as.vector(survPreds$risks));
            survPreds <- survpredict(currentModel,trainSet,trainSet,selectedFeaturesSet[[rept]]);
            csurvTrainPredictions <- cbind(trainSet[,theTime],trainSet[,theOutcome],rep(rept,nrow(trainSet)),survPreds$linearPredictors,as.vector(survPreds$risks));
            rownames(csurvTestPredictions) <- rownames(testSet)
            rownames(csurvTrainPredictions) <- rownames(trainSet)
            survTestPredictions <- rbind(survTestPredictions, csurvTestPredictions)
            survTrainPredictions <- rbind(survTrainPredictions, csurvTrainPredictions)
            # HR
            hr<-vector(mode="numeric", length=length(theVars))
            names(hr)<-theVars
            hr[names(survPreds$hr)]<-survPreds$hr
            survHR <- rbind(survHR,hr)
            names(survHR)<-theVars
            
          }
        }
        else
        {
          outx <- theData[sampleTrain,theOutcome];
          pred <- rep(NA,nrow(testSet));
          ctestPredictions <- cbind(testSet[,theOutcome],rep(rept,nrow(testSet)),pred);
          pred <- rep(NA,length(outx));
          ctrainPredictions <- cbind(outx,rep(rept,length(outx)),pred);
          rownames(ctestPredictions) <- rownames(testSet);
          rownames(ctrainPredictions) <- rownames(theData[sampleTrain,]);
          testPredictions <- rbind(testPredictions,ctestPredictions);
          trainPredictions <- rbind(trainPredictions,ctrainPredictions);
        }
      }
      else
      {
        outx <- theData[sampleTrain,theOutcome];
        pred <- rep(NA,nrow(testSet));
        ctestPredictions <- cbind(testSet[,theOutcome],rep(rept,nrow(testSet)),pred);
        pred <- rep(NA,length(outx));
        ctrainPredictions <- cbind(outx,rep(rept,length(outx)),pred);
        rownames(ctestPredictions) <- rownames(testSet);
        rownames(ctrainPredictions) <- rownames(theData[sampleTrain,]);
        testPredictions <- rbind(testPredictions,ctestPredictions);
        trainPredictions <- rbind(trainPredictions,ctrainPredictions);
      }
    }
    else
    {
      outx <- theData[sampleTrain,theOutcome];
      pred <- rep(NA,nrow(testSet));
      ctestPredictions <- cbind(testSet[,theOutcome],rep(rept,nrow(testSet)),pred);
      pred <- rep(NA,length(outx));
      ctrainPredictions <- cbind(outx,rep(rept,length(outx)),pred);
      rownames(ctestPredictions) <- rownames(testSet);
      rownames(ctrainPredictions) <- rownames(theData[sampleTrain,]);
      testPredictions <- rbind(testPredictions,ctestPredictions);
      trainPredictions <- rbind(trainPredictions,ctrainPredictions);
    }
    MADERROR[rept] = 0;
    avgsel <- avgsel + length(selectedFeaturesSet[[rept]])
    if ((!is.null(testPredictions) && length(rownames(testPredictions)) > 3 ))
    {
      boxstaTest <- try(boxplot(as.numeric(as.character(testPredictions[,3]))~rownames(testPredictions),plot = FALSE));
      if (!inherits(boxstaTest, "try-error"))
      {
        medianTest <- cbind(fullSet[boxstaTest$names,theOutcome],boxstaTest$stats[3,])
        tb <- table(rownames(testPredictions));
        MADERROR[rept] = mean(abs(medianTest[,1]-medianTest[,2]));
        if ((rept %% 10) == 0)
        {
          cat(rept," Tested:",nrow(medianTest),"Avg. Selected:",avgsel/rept,"Min Tests:",min(tb),"Max Tests:",max(tb),"Mean Tests:",mean(tb),". MAD:",MADERROR[rept],"\n");
        }
      }
    }
  }
  #	cat("done ",nrow(testPredictions),":",ncol(testPredictions),"\n");
  medianTest <- NULL;
  medianTrain <- NULL;
  boxstaTest <- NULL;
  boxstaTrain <- NULL;
  #Surv
  medianLinearPredictorsSurvTest <- NULL;
  medianLinearPredictorsSurvTrain <- NULL;
  boxstaLinearPredictorsSurvTest <- NULL;
  boxstaLinearPredictorsSurvTrain <- NULL;
  medianRisksSurvTest <- NULL;
  medianRisksSurvTrain <- NULL;
  boxstaRisksSurvTest <- NULL;
  boxstaRisksSurvTrain <- NULL;
  boxstaHRSurvTest <- NULL;
  boxstaHRSurvTrain <- NULL;
  
  medianSurvTest <- NULL;
  medianSurvTrain <- NULL;
  jaccard.sm <- NULL;
  featureFrequency <- NULL;
  if ((!is.null(testPredictions) && length(rownames(testPredictions)) > 3 ))
  {
    if (ncol(testPredictions) == 3)
    {
      colnames(testPredictions) <- c("Outcome","Model","Prediction");
      colnames(trainPredictions) <- c("Outcome","Model","Prediction");
    }
    boxstaTest <- try(boxplot(as.numeric(as.character(testPredictions[,3]))~rownames(testPredictions),plot = FALSE));
    if (!inherits(boxstaTest, "try-error"))
    {
      medianTest <- cbind(fullSet[boxstaTest$names,theOutcome],boxstaTest$stats[3,])
      rownames(medianTest) <- boxstaTest$names
    }
    else
    {
      warning("boxplot test failed");
      medianTest <- cbind(fullSet[,theOutcome],rep(0,nrow(fullSet)));
      rownames(medianTest) <- rownames(fullSet);
    }
    colnames(medianTest) <- c("Outcome","Median");
    
    boxstaTrain <- try(boxplot(as.numeric(as.character(trainPredictions[,3]))~rownames(trainPredictions),plot = FALSE));
    if (!inherits(boxstaTrain, "try-error"))
    {
      medianTrain  <- cbind(trainPredictions[boxstaTrain$names,1],boxstaTrain$stats[3,])
      rownames(medianTrain) <- boxstaTrain$names
    }
    else
    {
      warning("boxplot train failed");
      medianTrain <- cbind(theData[,theOutcome],rep(0,nrow(theData)));
      rownames(medianTrain) <- rownames(theData);
    }
    
    colnames(medianTrain) <- c("Outcome","Median");
    trainSamplesSets$repetitions <- repetitions;
    if (length(selectedFeaturesSet)>1) 
    {
      jaccard.sm <- jaccard(selectedFeaturesSet);
      featureFrequency <- table(unlist(selectedFeaturesSet));
      featureFrequency <- featureFrequency[order(-featureFrequency)];
    }
  }
  
  
  #[,theTime],testSet[,theOutcome]
  medianSurvTest <- data.frame(matrix(0, ncol = 6, nrow = nrow(theData)))
  colnames(medianSurvTest) <- c("Times","Outcome","LinearPredictorsMedian","RisksMedian");
  rownames(medianSurvTest) <- rownames(theData)
  medianSurvTest[,1] = theData[,theTime]
  medianSurvTest[,2] = theData[,theOutcome]
  
  medianSurvTrain <- data.frame(matrix(0, ncol = 6, nrow = nrow(theData)))
  colnames(medianSurvTrain) <- c("Times","Outcome","LinearPredictorsMedian","RisksMedian");
  rownames(medianSurvTrain) <- rownames(theData)
  medianSurvTrain[,1] = theData[,theTime]
  medianSurvTrain[,2] = theData[,theOutcome]
  
  # #Surv medians and boxsta
  if (!is.null(survTestPredictions) && length(rownames(survTestPredictions))>3)
  {
    if (ncol(survTestPredictions) == 7)
    {
      colnames(survTestPredictions) <- c("Times","Outcome","Model","LinearPredictors","Risks");
      colnames(survTrainPredictions) <- c("Times","Outcome","Model","LinearPredictors","Risks");
    }
    
    #   ######################Linear Predictors####################################  
    boxstaLinearPredictorsSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,5]))~rownames(survTestPredictions),plot = FALSE));
    if (!inherits(boxstaLinearPredictorsSurvTest, "try-error"))
    {
      medianSurvTest[boxstaLinearPredictorsSurvTest$names,4] <- boxstaLinearPredictorsSurvTest$stats[3,]
    }
    
    boxstaLinearPredictorsSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,5]))~rownames(survTrainPredictions),plot = FALSE));
    if (!inherits(boxstaLinearPredictorsSurvTrain, "try-error"))
    {
      medianSurvTrain[boxstaLinearPredictorsSurvTrain$names,4] <- boxstaLinearPredictorsSurvTrain$stats[3,]
    }
    
    #   ######################Risks####################################  
    boxstaRisksSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,7]))~rownames(survTestPredictions),plot = FALSE));
    if (!inherits(boxstaRisksSurvTest, "try-error"))
    {
      medianSurvTest[boxstaRisksSurvTest$names,6] <- boxstaRisksSurvTest$stats[3,]
    }
    
    boxstaRisksSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,7]))~rownames(survTrainPredictions),plot = FALSE));
    if (!inherits(boxstaRisksSurvTrain, "try-error"))
    {
      medianSurvTrain[boxstaRisksSurvTrain$names,6] <- boxstaRisksSurvTrain$stats[3,]
    }
  }
  
  results <- list(testPredictions = testPredictions,
                  trainPredictions = trainPredictions,
                  survTestPredictions = survTestPredictions,
                  survTrainPredictions = survTrainPredictions,
                  medianTest = medianTest,
                  medianTrain = medianTrain,
                  boxstaTest = boxstaTest,
                  boxstaTrain = boxstaTrain,
                  survMedianTest = medianSurvTest,
                  survMedianTrain = medianSurvTrain,
                  survBoxstaTest = list(
                                        linearPredictors=boxstaLinearPredictorsSurvTest,
                                        risks=boxstaRisksSurvTest),
                  survBoxstaTrain = list(
                                         linearPredictors=boxstaLinearPredictorsSurvTrain,
                                         risks=boxstaRisksSurvTrain),
                  survHR = survHR,
                  trainSamplesSets = trainSamplesSets,
                  selectedFeaturesSet = selectedFeaturesSet,
                  featureFrequency = featureFrequency,
                  formula.list = formula.list,
                  jaccard = jaccard.sm,
                  theTimes = theTimes,
                  MADERROR = MADERROR,
                  repetitions=repetitions
  );
  class(results) <- "RandomHOCV"
  return (results);
}
