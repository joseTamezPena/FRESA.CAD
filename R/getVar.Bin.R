getVar.Bin <-
function (object,data,Outcome="Class", type = c("LOGIT", "LM","COX"),testData=NULL,callCpp=TRUE) 
{

#	print(colnames(data));
#	print(all.vars(object$formula));
#	print(Outcome);
	data <- data[,all.vars(object$formula)];
	ttestdata <- !is.null(testData);
#	if (is.null(testData))
#	{
#		testData <- data;
#	}
#	print(colnames(testData));
	
  
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
#		cat(frm1,"\n")
	}


	model_zidi <- numeric(sizevec);
	model_znri <- numeric(sizevec);
	model_idi <- numeric(sizevec);
	model_nri <- numeric(sizevec);

	t.model_zidi <- numeric(sizevec);
	t.model_znri <- numeric(sizevec);
	t.model_idi <- numeric(sizevec);
	t.model_nri <- numeric(sizevec);
	test.accuracy <- numeric(sizevec);
	train.accuracy <- numeric(sizevec);
	redtest.accuracy <- numeric(sizevec);
	redtrain.accuracy <- numeric(sizevec);
	test.auc <- numeric(sizevec);
	train.auc <- numeric(sizevec);
	redtest.auc <- numeric(sizevec);
	redtrain.auc <- numeric(sizevec);
	fullTestAccuracy <- NA;
	fullTrainAccuracy <-NA;
	fullTestAUC <- NA;
	fullTrainAUC <-NA;

	ftmp <- formula(frm1);
	sizetrain <- nrow(data);
	traindataOutcome <- as.vector(data[,Outcome]);
	trainsizecases <- sum(traindataOutcome);
	
	if (ttestdata)
	{
		sizetest <- nrow(testData);
		testdataOutcome <- as.vector(testData[,Outcome]);
		testsizecases <- sum(testdataOutcome);
	}
	
	FullModel <- modelFitting(ftmp,data,type,fitFRESA=callCpp)
	if ((sizevec>0)&&( !inherits(FullModel, "try-error")))
	{
	
		FullPredict_train <- predict.fitFRESA(FullModel,data,'prob');
		sametrain <- 1.0*(traindataOutcome == 1.0*(FullPredict_train>=0.5));
		fullTrainAccuracy <- sum(sametrain)/sizetrain;
		fullTrainAUC <- 0.5*(sum(traindataOutcome*sametrain)/trainsizecases+sum((traindataOutcome==0)*sametrain)/(sizetrain-trainsizecases));

		if (ttestdata)
		{
			FullPredict <- predict.fitFRESA(FullModel,testData, 'prob');
			sametest <- 1.0*(testdataOutcome == 1.0*(FullPredict>=0.5));
			fullTestAccuracy <- sum(sametest)/sizetest;
			fullTestAUC  <- 0.5*(sum(testdataOutcome*sametest)/testsizecases+sum((testdataOutcome==0)*sametest)/(sizetest-testsizecases));
		}
		else
		{
			fullTestAccuracy <- fullTrainAccuracy;
			fullTestAUC  <- fullTrainAUC;
		}
		
		for ( i in 1:sizevec)
		{
			iwhere <- i;
		
			frm1 = outCome;
				for ( j in 1:sizevec)
				{
					if (i!=j)
					{
						frm1 <- paste(frm1,paste("+",termList[j]));
					}
				}
			ftmp <- formula(frm1);
#			cat("Red Formula:",termList[i],"->",frm1,"\n")
			redModel <- modelFitting(ftmp,data,type,fitFRESA=callCpp)
			if (inherits(redModel, "try-error"))
			{
				redModel <- FullModel;
			}

			redPredict_train <- predict.fitFRESA(redModel,data,'prob');
			if (ttestdata) {redPredict <- predict.fitFRESA(redModel,testData,'prob');}
			else {redPredict <- redPredict_train;}
			
			if (!callCpp)
			{
				iprob_t <- improveProb(redPredict_train,FullPredict_train,traindataOutcome);
				if (ttestdata) {iprob <- improveProb(redPredict,FullPredict,testdataOutcome);}
				else {iprob <- iprob_t;}
			}
			else
			{
				iprob_t <-.Call("improveProbCpp",redPredict_train,FullPredict_train,traindataOutcome);
				if (ttestdata) {iprob <- .Call("improveProbCpp",redPredict,FullPredict,testdataOutcome);}
				else {iprob <- iprob_t;}
			}

			uniModel <- modelFitting(formula(paste(outCome,"+",termList[i])),data,type,fitFRESA=callCpp);
			
			redPredict_train <- 1.0*(redPredict_train>=0.5);
			sametrain <- 1.0*(traindataOutcome == redPredict_train);
			redtrain.accuracy[iwhere] <- sum(sametrain)/sizetrain;
			redtrain.auc[iwhere] <- 0.5*(sum(traindataOutcome*sametrain)/trainsizecases+sum((traindataOutcome==0)*sametrain)/(sizetrain-trainsizecases));
			uniPredict_train <- 1.0*(predict.fitFRESA(uniModel,data,'prob')>=0.5);
			sametrain <- 1.0*(traindataOutcome == uniPredict_train);
			train.accuracy[iwhere] <- sum(sametrain)/sizetrain;
			train.auc[iwhere] <- 0.5*(sum(traindataOutcome*sametrain)/trainsizecases+sum((traindataOutcome==0)*sametrain)/(sizetrain-trainsizecases));
			
			if (ttestdata) 
			{
				redPredict <- 1.0*(redPredict>=0.5);
				sametest <- 1.0*(testdataOutcome == redPredict);
				redtest.accuracy[iwhere] <- sum(sametest)/sizetest;
				redtest.auc[iwhere] <- 0.5*(sum(testdataOutcome*sametest)/testsizecases+sum((testdataOutcome==0)*sametest)/(sizetest-testsizecases));
				uniPredict <- 1.0*(predict.fitFRESA(uniModel,testData,'prob')>=0.5);
				sametest <- 1.0*(testdataOutcome == uniPredict);
				test.accuracy[iwhere] <- sum(sametest)/sizetest;
				test.auc[iwhere] <- 0.5*(sum(testdataOutcome*sametest)/testsizecases+sum((testdataOutcome==0)*sametest)/(sizetest-testsizecases));
			}
			else
			{
				redtest.accuracy[iwhere] <- redtrain.accuracy[iwhere];
				redtest.auc[iwhere] <- redtrain.auc[iwhere] 
				test.accuracy[iwhere] <- train.accuracy[iwhere];
				test.auc[iwhere] <- train.auc[iwhere];
			}
		
			model_zidi[iwhere] <- iprob$z.idi;
			model_idi[iwhere] <- iprob$idi;
			model_nri[iwhere] <- iprob$nri;
			model_znri[iwhere] <- iprob$z.nri;

			t.model_zidi[iwhere] <- iprob_t$z.idi;
			t.model_idi[iwhere] <- iprob_t$idi;
			t.model_nri[iwhere] <- iprob_t$nri;
			t.model_znri[iwhere] <- iprob_t$z.nri;
		}

	}
	else
	{
#		cat("Error:",frm1,"\n")
		if (sizevec>0)
		{
			for ( i in 1:sizevec)
			{
				iwhere <- i;

				uniModel <- modelFitting(formula(paste(outCome,"+",termList[i])),data,type,fitFRESA=callCpp);
				uniPredict_train <- 1.0*(predict.fitFRESA(uniModel,data,'prob')>=0.5);
				sametrain <- 1.0*(traindataOutcome == uniPredict_train);
				train.accuracy[iwhere] <- sum(sametrain)/sizetrain;
				train.auc[iwhere] <- 0.5*(sum(traindataOutcome*sametrain)/trainsizecases+sum((traindataOutcome==0)*sametrain)/(sizetrain-trainsizecases));

				if (ttestdata)
				{
					uniPredict <- 1.0*(predict.fitFRESA(uniModel,testData,'prob')>=0.5);
					sametest <- 1.0*(testdataOutcome == uniPredict);
					test.accuracy[iwhere] <- sum(sametest)/sizetest;
					test.auc[iwhere] <- 0.5*(sum(testdataOutcome*sametest)/testsizecases+sum((testdataOutcome==0)*sametest)/(sizetest-testsizecases));
				}
				else			
				{
					test.accuracy[iwhere] <- train.accuracy[iwhere];
					test.auc[iwhere]  <- train.auc[iwhere];
				}
			}
		}
		model_zidi <- rep(0,sizevec);
		model_idi <- model_zidi;
		model_nri <- model_zidi;
		model_znri <- model_zidi;

		t.model_zidi <- model_zidi;
		t.model_idi <- model_zidi;
		t.model_nri <- model_zidi;
		t.model_znri <- model_zidi;
	}
#	cat("end getvar.bin\n")
	 result <- list(
		 z.IDIs=t.model_zidi,
		 z.NRIs=t.model_znri,
		 IDIs=t.model_idi,
		 NRIs=t.model_nri,
		 testData.z.IDIs=model_zidi,
		 testData.z.NRIs=model_znri,
		 testData.IDIs=model_idi,
		 testData.NRIs=model_nri,
		 uniTrainAccuracy = train.accuracy,
		 uniTestAccuracy = test.accuracy,
		 redtestAccuracy = redtest.accuracy,
		 redtrainAccuracy = redtrain.accuracy,
		 uniTrainAUC = train.auc,
		 uniTestAUC = test.auc,
		 redtestAUC = redtest.auc,
		 redtrainAUC = redtrain.auc,
		 fullTestAccuracy = fullTestAccuracy,
		 fullTrainAccuracy = fullTrainAccuracy,
		 fullTestAUC = fullTestAUC,
		 fullTrainAUC = fullTrainAUC
	 );

    return (result)
}
