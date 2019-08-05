randomCV <-  function(theData = NULL, theOutcome = "Class",fittingFunction=NULL, trainFraction = 0.5, repetitions = 100,trainSampleSets=NULL,featureSelectionFunction=NULL,featureSelection.control=NULL,asFactor=FALSE,addNoise=FALSE,classSamplingType=c("Augmented","NoAugmented","Proportional"),...)
{


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
	
	survpredict <- function(currentModel,Dataset,TestDataset,selectedFeatures)
	{
	  theSurvData <-Dataset;
	  fclass <- class(currentModel)
	  if (length(fclass)>1) fclass <- fclass[1];
	  
	  if(length(selectedFeatures)>=nrow(theSurvData))
	  {
	    selectedFeatures <- head(selectedFeatures, nrow(theSurvData)-1)
	    warning("The number of selected features is bigger than the number of observations, the top will be used.")
	  }
	  
	  if (fclass == "FRESA_LASSO")
	  {
	    #Creating lasso object
	    baseformula <- as.character(theformula);
	    formulaCox <- as.formula(paste(paste(baseformula[2],"~"), paste(selectedFeatures, collapse='+')));
	    cox <- try(survival::coxph(formula=formulaCox, data=theSurvData));
	    #changing the coef to the ones with lasso
	    cox$coefficients<-lft
	  }
	  if (fclass == "fitFRESA")
	  {
	    #Creating lasso object
	    # theSurvData <- testSet
	    # selectedFeatures <- selectedFeaturesSet[[1]];
	    baseformula <- as.character(theformula);
	    formulaCox <- as.formula(paste(paste(baseformula[2],"~"), paste(selectedFeatures, collapse='+')));
	    cox <- try(survival::coxph(formula=formulaCox,data=theSurvData));
	    #changing the coef to the ones with lasso
	    cox$coefficients<-currentModel$bagging$bagged.model$coefficients[-c(1)];
	    #cox$coefficients<-currentModel$BSWiMS.model$back.model$coefficients[-c(1)]
	  }
	  if(fclass=="bess")
	  {
	    baseformula <- as.character(theformula);
	    formulaCox <- as.formula(paste(paste(baseformula[2],"~"), paste(selectedFeatures, collapse='+')));
	    cox <- try(survival::coxph(formula=formulaCox,data=theSurvData));
	    #changing the coef to the ones with lasso
	    names(currentModel$fit$bestmodel$coefficients)<-selectedFeatures;
	    cox$coefficients<-currentModel$fit$bestmodel$coefficients;
	  }
	  if(fclass=="coxph.null")
	  {
	    
	    baseformula <- as.character(theformula);
	    formulaCox <- as.formula(paste(paste(baseformula[2],"~"), paste(selectedFeatures, collapse='+')));
	    cox <- try(survival::coxph(formula=formulaCox,data=theSurvData));
	  }
	  
	  if (!inherits(cox,"try-error") && class(cox)!="list")
	  {
	    followUpTimes <- try(predict(cox,newdata=TestDataset,type="expected"))
	    if (!inherits(followUpTimes,"try-error"))
	    {
	      #Martingale resid 
	      martingaleResid <- (as.integer(as.matrix(TestDataset[theOutcome]))-1) - followUpTimes
	      #linear predictos
	      linearPredictors <- predict(cox,newdata=TestDataset,type="lp")
	      #risk
	      #risks <- predict(cox,type="risk",se.fit=TRUE)
	      risks <- predict(cox,newdata=TestDataset,type="risk")
	      
	      hr <- round(coef(summary(cox))[,2],3)
	      survPreds <- list(martingaleResid=martingaleResid,
	                        linearPredictors=linearPredictors,
	                        followUpTimes=followUpTimes,
	                        risks = risks,
	                        hr = hr)
	      return (survPreds)
	    }
	    else{
	      warning("Cox Fit Error");
	      return(list(martingaleResid=rep(NA,nrow(theSurvData)),
	                  linearPredictors=rep(NA,nrow(theSurvData)),
	                  followUpTimes=rep(NA,nrow(theSurvData)),
	                  risks = list(fit=rep(NA,nrow(theSurvData)),se.fit=rep(NA,nrow(theSurvData))),
	                  hr = rep(NA,nrow(theSurvData))));
	    }
	  }
	  else{
	    warning("Cox Fit Error");
	    return(list(martingaleResid=rep(NA,nrow(theSurvData)),
	                linearPredictors=rep(NA,nrow(theSurvData)),
	                followUpTimes=rep(NA,nrow(theSurvData)),
	                risks = list(fit=rep(NA,nrow(theSurvData)),se.fit=rep(NA,nrow(theSurvData))),
	                hr = rep(NA,nrow(theSurvData))));
	  }
	}
rpredict <-  function(currentModel,DataSet)
{
	 fclass <- class(currentModel)
    if(fclass=="bess"){
      pred <- try(predict(currentModel$fit,DataSet))
    }
    else{
      pred <- try(predict(currentModel,DataSet))
    }

	if (inherits(pred, "try-error"))
	{
		pred <- numeric(nrow(DataSet));
	}
	else
	{
		if (class(pred) == "list")
		{
			if (is.null(pred$posterior))
			{
				pred <-as.numeric(as.character(pred[[1]]));
			}
			else
			{
				pred <-as.numeric(pred$posterior[,2]);
			}
		}
		if (class(pred) == "factor")
		{
			pred <- as.numeric(as.character(pred));
		}
		if (class(pred) == "array")
		{
			pnames <- colnames(pred);
			pred <- pnames[apply(pred[,,1],1,which.max)];
			pred <- as.numeric(pred);
		}
		if (class(pred) == "matrix") 
		{
			if (ncol(pred)>1)
			{
				pnames <- colnames(pred);
				pred <- pnames[apply(pred,1,which.max)];
				pred <- as.numeric(pred);
			}
			else
			{
				pred <- as.vector(pred);
			}
		}
	}
	return (pred)
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
	if (class(theOutcome)=="formula")
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
	samplePerClass <- as.integer((nrow(theData)/classLen)*trainFraction);
	if (testClases)
	{
		ssubsets <- list();
		samplePerClass <- as.integer((nrow(theData)/classLen)*trainFraction);
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
  	
	theVars <-colnames(theData)[!colnames(theData) %in% c(as.character(theTime),as.character(theOutcome))];
	survHR <- data.frame(matrix(ncol = length(theVars), nrow = 0));
  	colnames(survHR) <- theVars;
	  
	if (!is.null(trainSampleSets))
	{
		repetitions <- trainSampleSets$repetitions
	}

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
					if (classSamplingType == "Proportional")
					{
						sampleTrain <- sample(nrow(ssubsets[[jind]]),as.integer(nrow(ssubsets[[jind]])*trainFraction));
					}
					else
					{
						maxfrac <- max(c(trainFraction,0.95));
						ssize <- min(c(nrow(ssubsets[[jind]])-1,as.integer(nrow(ssubsets[[jind]])*maxfrac))); # minimum training size
						if (samplePerClass > ssize)
						{
							sampleTrain <- sample(nrow(ssubsets[[jind]]),ssize);
							if (classSamplingType == "Augmented")
							{
								therest <- samplePerClass-ssize;
								nsample <- sample(ssize,therest,replace=TRUE);
								sampleTrain <- append(sampleTrain,sampleTrain[nsample]);
							}
						}
						else
						{
							sampleTrain <- sample(nrow(ssubsets[[jind]]),samplePerClass);
						}
					}
					trainSamplesSets[[tset]] <- sampleTrain;
				}
				else
				{
					sampleTrain <- trainSampleSets[[tset]];
				}
				tset <- tset + 1;
				
				trainSet <- rbind(trainSet,ssubsets[[jind]][sampleTrain,]);
				testSet <- rbind(testSet,ssubsets[[jind]][-unique(sampleTrain),]);

				jind <- jind + 1;
#				cat("Class: ",s," rows:",nrow(trainSet),"\n");
			}
		}
		else
		{
			if (is.null(trainSampleSets))
			{
				sampleTrain <- sample(nrow(theData),nrow(theData)*trainFraction)
				trainSamplesSets[[tset]] <- sampleTrain;
			}
			else
			{
				sampleTrain <- trainSampleSets[[tset]];
			}
			tset <- tset + 1;
			trainSet <- theData[sampleTrain,];
			testSet <- theData[-sampleTrain,];
		}

		selnames <- character();
		if (class(featureSelectionFunction) == "list")
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
			if (class(featureSelectionFunction) == "function")
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
				if (cols>0)
				{
					stdg <- apply(trainSet[,fnames],2,sd,na.rm = TRUE);
					stdg <- min(stdg);
					if (stdg == 0)
					{
						stdg <- 1.0e-5;
					}
					noiselevel <- 0.05*stdg;
					rows <- nrow(trainSet);
					trainSet[,fnames] <- trainSet[,fnames]+as.data.frame(matrix(noiselevel*rnorm(rows*cols),nrow=rows,ncol=cols));
				}
			}
			if ((testClases) && (asFactor))
			{
					trainSet[,theOutcome] <- as.factor(trainSet[,theOutcome]); 
#					print(selnames)
			}
#			cat(ncol(trainSet),":",nrow(trainSet),"\n")

			theTimes <- append(theTimes,system.time(currentModel <- try(fittingFunction(theformula,trainSet,...))));
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
					if (fclass == "FRESA_LASSO")
					{
					  cf <- coef(currentModel$fit, s = currentModel$s)
					  if (class(cf) == "list")
					  {
					    cenet <- as.matrix(cf[[1]]);
					    if (!is.null(cenet))
					    {
					      lft <- cenet[as.vector(cenet[,1] != 0),,drop=FALSE]
					      numbers<-as.numeric(lft)
					      names(numbers)<-rownames(lft)
					      if(isSurv)
					      {
					        if (length(lft)>0)
					        {
					          selectedFeaturesSet[[rept]] <- names(lft);
					        }
					      }
					      else{
					        if (length(lft)>1)
					        {
					          if(isSurv)
					          {
					            if (length(lft)>0)
					            {
					              selectedFeaturesSet[[rept]] <- names(lft);
					            }
					          }
					          else{
					            selectedFeaturesSet[[rept]] <- names(lft)[-1];
					          }
					          for (cl in 2:length(cf))
					          {
					            cenet <- as.matrix(cf[[cl]]);
					            lft <- cenet[as.vector(cenet[,1] != 0),,drop=FALSE]
					            numbers<-as.numeric(lft)
					            names(numbers)<-rownames(lft)
					            lft<-numbers
					            if(isSurv)
					            {
					              if (length(lft)>0)
					              {
					                selectedFeaturesSet[[rept]] <- append(selectedFeaturesSet[[rept]],names(lft));
					              }
					            }
					            else{
					              if (length(lft)>1)
					              {
					                selectedFeaturesSet[[rept]] <- append(selectedFeaturesSet[[rept]],names(lft)[-1]);
					              }
					            }
					          }
					          selectedFeaturesSet[[rept]] <- unique(selectedFeaturesSet[[rept]]);
					        }
					      }
					      
					    }
					  }else{
					    cenet <- as.matrix(cf);
					    lft <- cenet[as.vector(cenet[,1]!=0),,drop=FALSE]
					    numbers<-as.numeric(lft)
					    names(numbers)<-rownames(lft)
					    lft<-numbers
					    
					    if(isSurv)
					    {
					      if (length(lft)>0)
					      {
					        selectedFeaturesSet[[rept]] <- names(lft);
					      }
					    }
					    else{
					      if (length(lft)>1)
					      {
					        selectedFeaturesSet[[rept]] <- names(lft)[-1];
					      }
					    }
					  }
					}
					if (fclass == "bess")
					{
					  bessCoefficients <- currentModel$fit$bestmodel$coefficients
					  if (!is.null(bessCoefficients))
					  {
					    selectedFeaturesSet[[rept]] <- gsub("xbest","",names(bessCoefficients));
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
#				print(selectedFeaturesSet);
				if ((length(selectedFeaturesSet[[rept]])>0) || is.null(featureSelectionFunction))
				{
					pred <- rpredict(currentModel,testSet);
					ctestPredictions <- cbind(testSet[,theOutcome],rep(rept,nrow(testSet)),pred);
					pred <- rpredict(currentModel,trainSet);
					ctrainPredictions <- cbind(trainSet[,theOutcome],rep(rept,nrow(trainSet)),pred);
					rownames(ctestPredictions) <- rownames(testSet);
					rownames(ctrainPredictions) <- rownames(trainSet);
					testPredictions <- rbind(testPredictions,ctestPredictions);
					trainPredictions <- rbind(trainPredictions,ctrainPredictions);
					if (isSurv)
					{
						#SURVPREDICT
						survPreds <- survpredict(currentModel,trainSet,testSet,selectedFeaturesSet[[rept]]);
						csurvTestPredictions <- cbind(testSet[,theTime],testSet[,theOutcome],rep(rept,nrow(testSet)),as.vector(survPreds$martingaleResid),survPreds$linearPredictors,as.vector(survPreds$followUpTimes),as.vector(survPreds$risks));
						survPreds <- survpredict(currentModel,trainSet,trainSet,selectedFeaturesSet[[rept]]);
						csurvTrainPredictions <- cbind(trainSet[,theTime],trainSet[,theOutcome],rep(rept,nrow(trainSet)),as.vector(survPreds$martingaleResid),survPreds$linearPredictors,as.vector(survPreds$followUpTimes),as.vector(survPreds$risks));
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
	}
#	cat("done ",nrow(testPredictions),":",ncol(testPredictions),"\n");
	medianTest <- NULL;
	medianTrain <- NULL;
	boxstaTest <- NULL;
	boxstaTrain <- NULL;
	  #Surv
	medianMartingaleResidSurvTest <- NULL;
	medianMartingaleResidSurvTrain <- NULL;
	boxstaMartingaleResidSurvTest <- NULL;
	boxstaMartingaleResidSurvTrain <- NULL;
	medianLinearPredictorsSurvTest <- NULL;
	medianLinearPredictorsSurvTrain <- NULL;
	boxstaLinearPredictorsSurvTest <- NULL;
	boxstaLinearPredictorsSurvTrain <- NULL;
	medianFollowUpTimesSurvTest <- NULL;
	medianFollowUpTimesSurvTrain <- NULL;
	boxstaFollowUpTimesSurvTest <- NULL;
	boxstaFollowUpTimesSurvTrain <- NULL;
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
	if (!is.null(testPredictions))
	{
		if (ncol(testPredictions) == 3)
		{
			colnames(testPredictions) <- c("Outcome","Model","Prediction");
			colnames(trainPredictions) <- c("Outcome","Model","Prediction");
		}
		boxstaTest <- try(boxplot(as.numeric(as.character(testPredictions[,3]))~rownames(testPredictions),plot = FALSE));
		if (!inherits(boxstaTest, "try-error"))
		{
			medianTest <- cbind(theData[boxstaTest$names,theOutcome],boxstaTest$stats[3,])
			rownames(medianTest) <- boxstaTest$names
		}
		else
		{
			warning("boxplot test failed");
			medianTest <- cbind(theData[,theOutcome],rep(0,nrow(theData)));
			rownames(medianTest) <- rownames(theData);
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

	#Surv medians and boxsta
	if (!is.null(survTestPredictions))
	{
	  if (ncol(survTestPredictions) == 7)
	  {
	    colnames(survTestPredictions) <- c("Times","Outcome","Model","MartinGale","LinearPredictors","FollowUpTimes","Risks");
	    colnames(survTrainPredictions) <- c("Times","Outcome","Model","MartinGale","LinearPredictors","FollowUpTimes","Risks");
	  }
	  
	  ######################Martin Gale#################################### 
	  boxstaMartingaleResidSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,4]))~rownames(survTestPredictions),plot = FALSE));
	  if (!inherits(boxstaMartingaleResidSurvTest, "try-error"))
	  {
	    medianSurvTest <- cbind(theData[boxstaMartingaleResidSurvTest$names,theTime],theData[boxstaMartingaleResidSurvTest$names,theOutcome],boxstaMartingaleResidSurvTest$stats[3,])
	    rownames(medianSurvTest) <- boxstaMartingaleResidSurvTest$names
	  }else
	  {
	    warning("boxplot test failed");
	    medianSurvTest <- cbind(theData[,theTime],theData[,theOutcome],rep(0,nrow(theData)));
	    rownames(medianSurvTest) <- rownames(theData);
	  }
	  
	  boxstaMartingaleResidSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,4]))~rownames(survTrainPredictions),plot = FALSE));
	  if (!inherits(boxstaMartingaleResidSurvTrain, "try-error"))
	  {
	    medianSurvTrain  <- cbind(theData[boxstaMartingaleResidSurvTrain$names,theTime],theData[boxstaMartingaleResidSurvTrain$names,theOutcome],boxstaMartingaleResidSurvTrain$stats[3,])
	    rownames(medianSurvTrain) <- boxstaMartingaleResidSurvTrain$names
	  }else
	  {
	    warning("boxplot train failed");
	    medianSurvTrain <- cbind(theData[,theTime],theData[,theOutcome],rep(0,nrow(theData)));
	    rownames(medianSurvTrain) <- rownames(theData);
	  }
	  
	  ######################Linear Predictors####################################  
	  boxstaLinearPredictorsSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,5]))~rownames(survTestPredictions),plot = FALSE));
	  if (!inherits(boxstaLinearPredictorsSurvTest, "try-error"))
	  {
	    medianSurvTest <- cbind(medianSurvTest,boxstaLinearPredictorsSurvTest$stats[3,])
	  }else
	  {
	    warning("boxplot test failed");
	    medianSurvTest <- cbind(medianSurvTest,rep(0,nrow(theData)));
	  }
	  
	  
	  boxstaLinearPredictorsSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,5]))~rownames(survTrainPredictions),plot = FALSE));
	  if (!inherits(boxstaLinearPredictorsSurvTrain, "try-error"))
	  {
	    medianSurvTrain  <- cbind(medianSurvTrain,boxstaLinearPredictorsSurvTrain$stats[3,])
	  }else
	  {
	    warning("boxplot train failed");
	    medianSurvTrain  <- cbind(medianSurvTrain,rep(0,nrow(theData)));
	  }
	  
	  ######################Follow Up Times####################################  
	  boxstaFollowUpTimesSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,6]))~rownames(survTestPredictions),plot = FALSE));
	  if (!inherits(boxstaFollowUpTimesSurvTest, "try-error"))
	  {
	    medianSurvTest <- cbind(medianSurvTest,boxstaFollowUpTimesSurvTest$stats[3,])
	  }else
	  {
	    warning("boxplot test failed");
	    medianSurvTest <- cbind(medianSurvTest,rep(0,nrow(theData)));
	  }
	  
	  boxstaFollowUpTimesSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,6]))~rownames(survTrainPredictions),plot = FALSE));
	  if (!inherits(boxstaFollowUpTimesSurvTrain, "try-error"))
	  {
	    medianSurvTrain  <- cbind(medianSurvTrain,boxstaFollowUpTimesSurvTrain$stats[3,])
	  }else
	  {
	    warning("boxplot train failed");
	    medianSurvTrain  <- cbind(medianSurvTrain,rep(0,nrow(theData)));
	  }
	  
	  ######################Risks####################################  
	  boxstaRisksSurvTest <- try(boxplot(as.numeric(as.character(survTestPredictions[,7]))~rownames(survTestPredictions),plot = FALSE));
	  if (!inherits(boxstaRisksSurvTest, "try-error"))
	  {
	    medianSurvTest <- cbind(medianSurvTest,boxstaRisksSurvTest$stats[3,])
	  }else
	  {
	    warning("boxplot test failed");
	    medianSurvTest <- cbind(medianSurvTest,rep(0,nrow(theData)));
	  }
	  
	  boxstaRisksSurvTrain <- try(boxplot(as.numeric(as.character(survTrainPredictions[,7]))~rownames(survTrainPredictions),plot = FALSE));
	  if (!inherits(boxstaRisksSurvTrain, "try-error"))
	  {
	    medianSurvTrain  <- cbind(medianSurvTrain,boxstaRisksSurvTrain$stats[3,])
	  }else
	  {
	    warning("boxplot train failed");
	    medianSurvTrain  <- cbind(medianSurvTrain,rep(0,nrow(theData)));
	  }
	  
	  
	  colnames(medianSurvTest) <- c("Times","Outcome","MartinGaleMedian","LinearPredictorsMedian","FollowUpTimesMedian","RisksMedian");
	  colnames(medianSurvTrain) <- c("Times","Outcome","MartinGaleMedian","LinearPredictorsMedian","FollowUpTimesMedian","RisksMedian");
	}
  
  

	#	cat("done ",nrow(medianTest),":",ncol(medianTest),"\n");
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
					survBoxstaTest = list(martingaleResid=boxstaMartingaleResidSurvTest,
										linearPredictors=boxstaLinearPredictorsSurvTest,
										followUpTimes=boxstaFollowUpTimesSurvTest,
										risks=boxstaRisksSurvTest),
					survBoxstaTrain = list(martingaleResid=boxstaMartingaleResidSurvTrain,
											linearPredictors=boxstaLinearPredictorsSurvTrain,
											followUpTimes=boxstaFollowUpTimesSurvTrain,
											risks=boxstaRisksSurvTrain),
					survHR = survHR,
					trainSamplesSets = trainSamplesSets,
					selectedFeaturesSet = selectedFeaturesSet,
					featureFrequency = featureFrequency,
					jaccard = jaccard.sm,
					theTimes = theTimes,
					repetitions=repetitions
		);
	return (results);
}