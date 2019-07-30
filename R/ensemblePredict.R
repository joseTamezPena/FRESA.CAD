ensemblePredict <-
function (formulaList,trainData,testData=NULL, predictType = c("prob", "linear"),type = c("LOGIT", "LM","COX","SVM"),Outcome=NULL,nk = 0) 
{

#	cat("Median\n")
#		print(formulaList);
	if (length(formulaList)==0)
	{
		cat("No formulas\n");
		warning("No formulas");
		result <- list(ensemblePredict=NULL,
		medianKNNPredict=NULL,predictions=NULL,KNNpredictions=NULL,wPredict=NULL)
		return (result)		
	}
	if (is.null(Outcome)) 
	{
		nk = -1;
#		for (i in 1:length(formulaList))
#		{
#			formulaList[i] <- paste(formulaList[i],"+1")
#		}
		varlist <- attr(terms(formula(formulaList[1])),"variables")
		dependent <- as.character(varlist[[2]])
		if (length(dependent)==3)
		{
			vOutcome = dependent[3];
		}
		else
		{
			vOutcome = dependent[1];
		}		
	}
	else
	{
		vOutcome=Outcome;
		for (i in 1:length(formulaList))
		{
			if (formulaList[i] != "=-=End=-=")
			{
				if (gregexpr(pattern ='~',as.character(formulaList[i]))[1]>0)
				{
					feat <- unlist(strsplit(as.character(formulaList[i]),"[~]"));
					if (nchar(feat[2])>0)
					{
						formulaList[i] <- feat[2];
					}
					else
					{
						formulaList[i] <- "1";
					}
				}
			}
		}
	}
	EquTrainSet <- trainData;
	minTrainSamples <- nrow(trainData);
	maxTrainSamples = minTrainSamples;
	casesample  <- NULL;
	controlsample <- NULL;
	noequalSets <- FALSE;
	nrowcases <- minTrainSamples
	nrowcontrols <- minTrainSamples
	if ((type == "LOGIT") || (type == "COX"))
	{
		casesample = subset(trainData,get(vOutcome)  == 1);
		controlsample = subset(trainData,get(vOutcome) == 0);
		trainCaseses <- casesample;
		trainControls <- controlsample;
		nrowcases <- nrow(casesample);
		nrowcontrols <- nrow(controlsample);
		
		minTrainSamples <- min(c(nrow(casesample),nrow(controlsample)));
		maxTrainSamples <- max(c(nrow(casesample),nrow(controlsample)));
		noequalSets <- (minTrainSamples < 0.80*maxTrainSamples);
#		if (noequalSets) warning("Case and Control sets no equal. Setting the size of the control and cases equal\n");
	}
	
	if (nk==0)
	{
		nk = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
	}
	if (is.null(testData))
	{
		medianKNN=NULL;
		out=NULL;
		KNNpredictions=NULL;
		outKNN=NULL;
		medianout <- vector(mode="numeric",length = nrow(trainData));
		wmedpredict <- medianout;
		cat("\n");
		for ( n in 1:nrow(trainData))
		{
			if ((n %% 10)==0) cat(".");
			mp <- ensemblePredict(formulaList,trainData[-n,],trainData[n,],predictType,type,Outcome,nk = -1)
			medianout[n] = mp$ensemblePredict[1]
			wmedpredict[n] = mp$wPredict[1]
			out <- rbind(out,cbind(trainData[n,vOutcome],mp$predictions[1,-1]));
		}
		rownames(out) <- rownames(trainData);
		cat("\n");
	}
	else
	{
		wmedpred <- numeric(nrow(testData))
		swts <- 0;
		theoutcome <- testData[,vOutcome];
		out <- NULL;
		if (!is.null(Outcome)) 
		{
			if (nchar(formulaList[[1]])>0)
			{
				ftmp <- formula(paste(Outcome,"~ 1+",formulaList[[1]]))
			}
			else
			{
				ftmp <- formula(paste(Outcome,"~ 1"));
			}
		}
		else
		{
			ftmp <- formula(formulaList[1])
		}
		if (noequalSets)
		{
			if (maxTrainSamples > nrowcases)  trainCaseses <- casesample[sample(1:nrowcases,maxTrainSamples,replace=TRUE),]
			if (maxTrainSamples > nrowcontrols)  trainControls <- controlsample[sample(1:nrowcontrols,maxTrainSamples,replace=TRUE),]
			EquTrainSet <- rbind(trainCaseses,trainControls)
		}
		if (nk>0)
		{
			outKNN <- cbind(theoutcome,getKNNpredictionFromFormula(ftmp,EquTrainSet,testData,Outcome=vOutcome,nk)$binProb)
			rownames(outKNN) <- rownames(testData);
		}
		else
		{
			outKNN <- NULL;
			medianKNN <- NULL;
		}

		bestmodel <- modelFitting(ftmp,EquTrainSet,type,fitFRESA=TRUE)
		if (inherits(bestmodel, "try-error"))
		{
			curprediction <- numeric(nrow(testData));
			wmedpred <- curprediction;
			swts <- 0;
		}
		else
		{		
			curprediction <- predict.fitFRESA(bestmodel,testData,predictType);
			trainOutcome <- EquTrainSet[,vOutcome];
			varOutcome <- var(trainOutcome);
			trainPrediction <- predict.fitFRESA(bestmodel,EquTrainSet,predictType);
			residual <- as.vector(abs(trainPrediction-trainOutcome));
			observations <- nrow(EquTrainSet);
			if (predictType[1]=="linear")
			{
				Rwts <- (varOutcome-sum(residual^2)/observations)/varOutcome; #Correlation
			}
			else
			{
				Rwts <- 2.0*(sum(1.0*(residual<0.5))/observations - 0.5); # 2*(ACC-0.5)
			}
			if (Rwts<=0) Rwts <- 1.0e-4;
			wmedpred <- Rwts*curprediction;
			swts <- Rwts;
		}




		totSamples <- cbind(theoutcome,curprediction);
		rownames(totSamples) <- rownames(testData);
		out <- totSamples;
		if (length(formulaList)>1)
		{
			for (i in 2:length(formulaList))
			{
				if (formulaList[i] != "=-=End=-=")
				{

					if (!is.null(Outcome))
					{		
						if (nchar(formulaList[[i]])>0)
						{
							ftmp <- formula(paste(Outcome,"~1+",formulaList[[i]]));
						}
						else
						{
							ftmp <- formula(paste(Outcome,"~1"));
						}
					}
					else
					{
						ftmp <- formula(formulaList[i]);
					}
					if (noequalSets && (runif(1)<0.20))
					{
						if (maxTrainSamples > nrowcases)  trainCaseses <- casesample[sample(1:nrowcases,maxTrainSamples,replace=TRUE),]
						if (maxTrainSamples > nrowcontrols)  trainControls <- controlsample[sample(1:nrowcontrols,maxTrainSamples,replace=TRUE),]
						EquTrainSet <- rbind(trainCaseses,trainControls)
					}
					if (nk>0) 
					{
						outKNN <- cbind(outKNN,getKNNpredictionFromFormula(ftmp,EquTrainSet,testData,Outcome=vOutcome,nk)$binProb);
					}
					if (length(all.vars(ftmp))>1)
					{
						fm <- modelFitting(ftmp,EquTrainSet,type,fitFRESA=TRUE)
						if (inherits(fm, "try-error"))
						{
							if (noequalSets)
							{
								if (maxTrainSamples > nrowcases)  trainCaseses <- casesample[sample(1:nrowcases,maxTrainSamples,replace=TRUE),]
								if (maxTrainSamples > nrowcontrols)  trainControls <- controlsample[sample(1:nrowcontrols,maxTrainSamples,replace=TRUE),]
								EquTrainSet <- rbind(trainCaseses,trainControls)
							}
							fm <- modelFitting(ftmp,EquTrainSet,type,fitFRESA=TRUE)
						}
						curprediction <- predict.fitFRESA(fm,testData,predictType);
						trainOutcome <- EquTrainSet[,vOutcome];
						varOutcome <- var(trainOutcome);
						trainPrediction <- predict.fitFRESA(fm,EquTrainSet,predictType);
						residual <- as.vector(abs(trainPrediction-trainOutcome));
						observations <- nrow(EquTrainSet);
						if (predictType[1]=="linear")
						{
							Rwts <- (varOutcome-sum(residual^2)/observations)/varOutcome; #Correlation
						}
						else
						{
							Rwts <- 2.0*(sum(1.0*(residual<0.5))/observations - 0.5); # 2*(ACC-0.5)
						}
						if (Rwts<=0) Rwts <- 1.0e-4;
						wmedpred <- wmedpred + Rwts*curprediction;
						swts <- swts + Rwts;
						out <- cbind(out,curprediction);
					}
					else
					{
						warning("No formula");
					}
				}
			}
			if (ncol(out)<3)
			{
				out <- cbind(out,curprediction);
			}
			out <- as.data.frame(out);
			medianout <- rowMedians(out[,-1],na.rm = TRUE);
			wmedpredict <- wmedpred;
			if (swts>0) wmedpredict <- wmedpred/swts;
			if (nk>0) 
			{
				outKNN <- as.data.frame(outKNN);
				medianKNN <- rowMedians(outKNN[,-1],na.rm = TRUE);
			}
		}
		else
		{
			out <- as.data.frame(out);
			medianout <- out[,-1];
			wmedpredict <- medianout;
			if (nk>0) 
			{
				outKNN <- as.data.frame(outKNN);
				medianKNN <- outKNN[,-1];
			}
		}
	}
	result <- list(ensemblePredict=medianout,
	medianKNNPredict=medianKNN,predictions=out,KNNpredictions=outKNN,wPredict=wmedpredict)
    return (result)
}