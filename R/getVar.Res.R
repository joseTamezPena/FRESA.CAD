getVar.Res <-
function (object,data,Outcome="Class", type=c("LM","LOGIT","COX"),testData=NULL,callCpp=TRUE) 
{
#	print(colnames(data));
#	print(all.vars(object$formula));
	data <- data[,all.vars(object$formula)];
	ttestData <- !is.null(testData);
  
	samples <- nrow(data);
	varsList <- unlist(as.list(attr(terms(object),"variables")))
	termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
	sizevec <- length(termList);
	
	outCome = paste(varsList[2]," ~ 1");
	frm1 = outCome;
	if (sizevec>0)
	{
		for ( i in 1:sizevec)
		{
			frm1 <- paste(frm1,paste("+",termList[i]));
		}
	}
	ftmp <- formula(frm1);
	FullModel <- modelFitting(ftmp,data,type,fitFRESA=callCpp);

	FullResiduals <- residualForFRESA(FullModel,data,Outcome);
	FullTrainMSE <- mean(FullResiduals^2);
	dataOutcome <- data[,Outcome];

	if (ttestData)
	{
		testResiduals <- residualForFRESA(FullModel,testData,Outcome);
		FullTestMSE <- mean(testResiduals^2);
		testdataOutcome <- testData[,Outcome];
	}
	else
	{
		FullTestMSE <- FullTrainMSE;
	}

	

	model_tpvalue <- numeric(sizevec);
	model_bpvalue <- numeric(sizevec);
	model_neri <- numeric(sizevec);
	model_wpvalue <- numeric(sizevec);
	model_fpvalue <- numeric(sizevec);
	testmodel_tpvalue <- numeric(sizevec);
	testmodel_bpvalue <- numeric(sizevec);
	testmodel_neri <- numeric(sizevec);
	testmodel_wpvalue <- numeric(sizevec);
	testmodel_fpvalue <- numeric(sizevec);
	unitest.MSE <- numeric(sizevec);
	unitrain.MSE <- numeric(sizevec);
	redtest.MSE <- numeric(sizevec);
	redtrain.MSE <- numeric(sizevec);
	obsser <- (nrow(data)>2) # at least 3 samples for improvement analysis

	if (sizevec>0)
	{
		for ( i in 1:sizevec)
		{
			iwhere <- i;
			frm1 = outCome;
			for ( j in 1:sizevec)
			{
				if (i!=j)
				{
					frm1 <- paste(frm1,"+",termList[j]);
				}
			}
			ftmp <- formula(frm1);
			redModel <- modelFitting(ftmp,data,type,fitFRESA=callCpp)

			if ( inherits(redModel, "try-error"))
			{
				redModel <- FullModel;
			}
			
			redResiduals <- residualForFRESA(redModel,data,Outcome);
			if (ttestData) 
			{
				redTestResiduals <- residualForFRESA(redModel,testData,Outcome);
			}
			else
			{
				redTestResiduals <- redResiduals;
			}

			if (obsser)
			{
				if (!callCpp)
				{
					iprob <- improvedResiduals(redResiduals,FullResiduals);
					if (ttestData) testiprob <- improvedResiduals(redTestResiduals,testResiduals);
				}
				else
				{
					iprob <- .Call("improvedResidualsCpp",redResiduals,FullResiduals," ",0);
					if (ttestData) 
					{
						testiprob <- .Call("improvedResidualsCpp",redTestResiduals,testResiduals," ",samples);
					}
					else 
					{
						testiprob <- iprob;
					}
				}
				model_tpvalue[iwhere] <- iprob$tP.value;
				model_bpvalue[iwhere] <- iprob$BinP.value;
				model_wpvalue[iwhere] <- iprob$WilcoxP.value;
				model_fpvalue[iwhere] <- iprob$FP.value;
				model_neri[iwhere] <- iprob$NeRI;
				testmodel_tpvalue[iwhere] <- testiprob$tP.value;
				testmodel_bpvalue[iwhere] <- testiprob$BinP.value;
				testmodel_wpvalue[iwhere] <- testiprob$WilcoxP.value;
				testmodel_fpvalue[iwhere] <- testiprob$FP.value;
				testmodel_neri[iwhere] <- testiprob$NeRI;
			}
			redtrain.MSE[iwhere] = mean(redTestResiduals^2);
			
			uniModel <- modelFitting(formula(paste(outCome,"+",termList[i])),data,type,fitFRESA=callCpp);
			uniPredict_train <- predict.fitFRESA(uniModel,data,'linear');
			unitrain.MSE[iwhere] <- mean((uniPredict_train-dataOutcome)^2);
			if (ttestData) 
			{
				uniPredict <- predict.fitFRESA(uniModel,testData,'linear');
				unitest.MSE[iwhere] <- mean((uniPredict-testdataOutcome)^2);
				redtest.MSE[iwhere] = mean(redResiduals^2);
			}
			else
			{
				unitest.MSE[iwhere] <- unitrain.MSE[iwhere];
				redtest.MSE[iwhere] = redtrain.MSE[iwhere];
			}			
		}
	}
#		cat("End\n");

	 result <- list(
	 tP.value=model_tpvalue,
	 BinP.value=model_bpvalue,
	 WilcoxP.value=model_wpvalue,
	 FP.value=model_fpvalue,
	 NeRIs=model_neri,
	 testData.tP.value=testmodel_tpvalue,
	 testData.BinP.value=testmodel_bpvalue,
	 testData.WilcoxP.value=testmodel_wpvalue,
	 testData.FP.value=testmodel_fpvalue,
	 testData.NeRIs=testmodel_neri,
	 unitestMSE = unitest.MSE,
	 unitrainMSE = unitrain.MSE,
	 redtestMSE = redtest.MSE,
	 redtrainMSE = redtrain.MSE,
	 FullTrainMSE = FullTrainMSE,
	 FullTestMSE = FullTestMSE
	 );

    return (result)
}
