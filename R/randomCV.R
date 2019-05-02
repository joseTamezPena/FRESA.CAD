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

rpredict <-  function(currentModel,DataSet)
{
#	if (!is.null(currentModel$coefficients)) print(currentModel$coefficients)
	pred <- try(predict(currentModel,DataSet))
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

	theformula <- NULL;
	varsmod <- NULL;
	if (class(theOutcome)=="formula")
	{
		theformula <- theOutcome;
		varsmod <- all.vars(theformula);
		theOutcome <- varsmod[1];
		if (sum(str_count(theformula,"Surv")) > 0)
		{
			theOutcome <- varsmod[2];
		}
	}
	else
	{
		varsmod <- theOutcome;
		theformula <- formula(paste(theOutcome,"~ ."));
	}



	theClasses <- as.numeric(names(table(theData[,theOutcome])));
	classLen <- length(theClasses);
	selectedFeaturesSet <- list();
	if(classLen < 10)
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
		if(classLen < 10)
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
					frank <- featureSelectionFunction(trainSet,theOutcome)
					if (length(frank)>0)
					{
						selectedFeaturesSet[[rept]] <- names(frank);
					}
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
			if ((classLen < 10) && (asFactor))
			{
					trainSet[,theOutcome] <- as.factor(trainSet[,theOutcome]); 
#					print(selnames)
			}
#			cat(ncol(trainSet),":",nrow(trainSet),"\n")

			theTimes <- append(theTimes,system.time(currentModel <- try(fittingFunction(theformula,trainSet,...))));
			if ( inherits(currentModel, "try-error"))
			{
				if ((classLen < 10) && (!asFactor))
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
								lft <- cenet[as.vector(cenet[,1] != 0),]
								if (length(lft)>1)
								{
									selectedFeaturesSet[[rept]] <- names(lft)[-1];
									for (cl in 2:length(cf))
									{
										cenet <- as.matrix(cf[[cl]]);
										lft <- cenet[as.vector(cenet[,1] != 0),]
										if (length(lft)>1)
										{
											selectedFeaturesSet[[rept]] <- append(selectedFeaturesSet[[rept]],names(lft)[-1]);
										}
									}
									selectedFeaturesSet[[rept]] <- unique(selectedFeaturesSet[[rept]]);
								}
							}
						}
						else
						{
							cenet <- as.matrix(cf);
							lft <- cenet[as.vector(cenet[,1] != 0),]
							if (length(lft)>1)
							{
								selectedFeaturesSet[[rept]] <- names(lft)[-1];
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
#	cat("done ",nrow(medianTest),":",ncol(medianTest),"\n");
	results <- list(testPredictions = testPredictions,
					trainPredictions = trainPredictions,
					medianTest = medianTest,
					medianTrain = medianTrain,
					boxstaTest = boxstaTest,
					boxstaTrain = boxstaTrain,
					trainSamplesSets = trainSamplesSets,
					selectedFeaturesSet = selectedFeaturesSet,
					featureFrequency = featureFrequency,
					jaccard = jaccard.sm,
					theTimes = theTimes,
					repetitions=repetitions
		);
	return (results);
}
