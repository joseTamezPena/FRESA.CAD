crossValidationFeatureSelection_Res <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),startOffset=0,elimination.bootstrap.steps=100,trainFraction=0.67,trainRepetition=9,setIntersect=1,unirank=NULL,print=TRUE,plots=TRUE,lambda="lambda.1se",equivalent=FALSE,bswimsCycles=10,usrFitFun=NULL,featureSize=0)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 


	enetSamples <- NULL;
	enetTrainSamples <- NULL;
	totSamples <- NULL;
	totTrainSamples <- NULL;
	Full.totSamples <- NULL;
	Full.totTrainSamples <- NULL;
	uniTrainMSS <- NULL;
	uniTestMSS <- NULL;
	
	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);

	filter.p.value = 2.0*pvalue;
	if (length(pvalue)>1)
	{
		filter.p.value = pvalue[2];
		pvalue=pvalue[1];
	}

#	cat(type,"\n")

	Fullsammples <- nrow(data);
	if ( K > Fullsammples) K=Fullsammples

	if (!is.null(unirank))
	{
		uprank <- update.uniRankVar(unirank,data=data,FullAnalysis=FALSE)
		variableList <- uprank$orderframe;
	}

	if (size>nrow(variableList)) size <- nrow(variableList);
	if (size<5) size=min(c(5,nrow(variableList)));

	shortVarList <- as.vector(variableList[1:size,1]);

	varlist <- paste(Outcome,"~1");
	for (i in 1:length(shortVarList))
	{
		varlist <- paste(varlist,shortVarList[i],sep="+");
	}
	varlist <-  all.vars(formula(varlist))[-1];
	shortVarList <- varlist[varlist %in% colnames(data)];
	enetshortVarList <- shortVarList;
	
    LASSOVariables <- NULL;
	if (type=="LM")
	{
#		print(enetshortVarList);
		Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="gaussian"));
	}
	else
	{
		Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
	}
	if (inherits(Fullenet, "try-error"))
	{
		cat("enet Error")
		Fullenet <- NULL;
	}
	else
	{
		cenet <- as.matrix(coef(Fullenet,s=lambda))
#		print(LASSOVariables <- list(names(cenet[as.vector(cenet[,1] !=0 ),])))
		lanames <- names(cenet[as.vector(cenet[,1] != 0),])
		print(LASSOVariables <- paste(lanames[lanames !=  "(Intercept)" ],collapse=" + "))
	}

	if (type!="COX")
	{
		baseformula <- paste(Outcome,"~",covariates);
		abaseformula <- paste(Outcome,"~ .");
		extvar <- Outcome;
	}
	else
	{
		baseformula <- paste("Surv(",timeOutcome,",",Outcome,") ~",covariates);
		abaseformula <- paste("Surv(",timeOutcome,",",Outcome,") ~ .");
		extvar <- c(timeOutcome,covariates);
	}


	CVselection.pValue <- pvalue;

	FULLBSWiMS.models <- BSWiMS.model(formula=baseformula,data=data,type=type,testType=testType,pvalue=pvalue,variableList=variableList,size=size,loops=loops,elimination.bootstrap.steps=elimination.bootstrap.steps,fraction=fraction,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize);
	Full_CurModel_S <- FULLBSWiMS.models$forward.model;
	Full_UCurModel_S <- FULLBSWiMS.models$update.model;
	Full_redCurmodel_S <- FULLBSWiMS.models$BSWiMS.model;

	cat ("Update   :",Full_UCurModel_S$formula,"\n")
	cat ("At  RMSE :",Full_redCurmodel_S$at.RMSE.formula,"\n")
	cat ("B:SWiMS  :",Full_redCurmodel_S$back.formula,"\n")
	
		
	
	formulas <- vector();
	AtOptFormulas <- vector();
	ForwardFormulas <- vector();
	baggFormulas <- vector();
	equiFormulas <- vector();
	allBSWIMSFormulas <- character();
	vtrainRMS <- vector();
	vblindRMS <- vector();

	vtrainSpearman <- vector();
	vtrainPearson <- vector();

	FullvtrainRMS <- vector();
	FullvblindRMS <- vector();
	FullvtrainSpearman <- vector();
	FullvtrainPearson <- vector();

	blindFoldPearson <-  vector();
	blindFoldSpearman <- vector();
	blindFoldCstat <- vector();
	blindFoldMS <- vector();
	CVBlindPearson  <- vector();
	CVBlindSpearman  <- vector();
	CVBlindRMS <- vector();

	inserted = 0;
	lastmodel = 0;
	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
			sampleFolds <- cvTools::cvFolds(nrow(data), K,1, "random");
		}
				

		TrainSet <- data[sampleFolds$subsets[sampleFolds$which != j,],];
		BlindSet <- data[sampleFolds$subsets[sampleFolds$which == j,],];
		
		outindex <- grep(Outcome, colnames(TrainSet))-1;
		tsize <- nrow(TrainSet);
		ETraningSet <-as.data.frame(.Call("equalizedSampling", as.matrix(TrainSet),outindex,5));
		ETraningSet <- ETraningSet[sample(nrow(ETraningSet),tsize),]
		colnames(ETraningSet) <- colnames(TrainSet);
		if (plots)
		{
			hist(TrainSet[,Outcome],breaks=5,main="Before Equalization");
			hist(ETraningSet[,Outcome],breaks=5,main="After Equalization");
		}

		blindsampleidx <- as.vector(rownames(BlindSet));
		sampleidx <- as.vector(rownames(TrainSet));

		if (!is.null(unirank))
		{
			variableList <- update.uniRankVar(unirank,data=TrainSet,FullAnalysis=FALSE)$orderframe;

			unitPvalues <- (1.0-pnorm(variableList$ZUni));
			names(unitPvalues) <- variableList$Name;
			if (unirank$categorizationType == "Raw")
			{
				adjPvalues <- p.adjust(unitPvalues,"BH")
				gadjPvalues <- adjPvalues[adjPvalues < 2*filter.p.value]			
				noncornames <- correlated_Remove(data,names(gadjPvalues),thr=0.99);
				if (length(noncornames) > 1) featureSize <- featureSize*length(noncornames)/length(gadjPvalues);
#				cat(length(noncornames),":",length(gadjPvalues),":",length(noncornames)/length(gadjPvalues),"\n");
			}
			filter.z.value <- abs(qnorm(filter.p.value))
			varMax <- sum(variableList$ZUni >= filter.z.value);
			pvarMax <- sum(p.adjust(unitPvalues,"BH") < 2*filter.p.value);
			cat(ncol(TrainSet),":",nrow(variableList),": Unadjusted size:",varMax," Adjusted Size:",pvarMax,"\n")
			size= min(c(pvarMax,varMax));
			if (size<5) size=min(c(5,nrow(variableList)));

			shortVarList <- as.vector(variableList[1:size,1]);

			varlist <- paste(Outcome,"~1");
			for (nn in 1:length(shortVarList))
			{
				varlist <- paste(varlist,shortVarList[nn],sep="+");
			}
			varlist <-  all.vars(formula(varlist))[-1];
			shortVarList <- varlist[varlist %in% colnames(TrainSet)];
			
#			print(shortVarList)
		}




		cat("Samples Train :",nrow(TrainSet),"Equalized Samples Train :",nrow(ETraningSet),"Samples Test :",nrow(BlindSet),"\n");
		cat ("Loop :",i,"\n")

		if (!is.null(Fullenet))
		{
			if (type=="LM")
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(ETraningSet[,shortVarList]),as.vector(ETraningSet[,Outcome]),family="gaussian"));
			}
			else
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(ETraningSet[,shortVarList]),as.vector(ETraningSet[,Outcome]),family="binomial"));
			}
			cenet <- as.matrix(coef(foldenet,s=lambda))
#			LASSOVariables[[i+1]] <- names(cenet[as.vector(cenet[,1] != 0),])
			lanames <- names(cenet[as.vector(cenet[,1] != 0),])
			LASSOVariables <- append(LASSOVariables,paste(lanames[lanames !=  "(Intercept)" ],collapse=" + "))
			if (i == 1)
			{
				enetSamples <- cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList]),s=lambda),i);
				enetTrainSamples <- cbind(TrainSet[,Outcome],predict(foldenet,as.matrix(TrainSet[,shortVarList]),s=lambda),i);
			}
			else
			{
				enetSamples <- rbind(enetSamples,cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList]),s=lambda),i));
				enetTrainSamples <- rbind(enetTrainSamples,cbind(TrainSet[,Outcome],predict(foldenet,as.matrix(TrainSet[,shortVarList]),s=lambda),i));
			}
##			print(LASSOVariables)
		}

		
		
		par(mfrow=c(1,1))

#		cat(type,"\n")

#		cat(Full_redCurmodel_S$back.formula," <-Back formula\n");


		BSWiMS.models <- BSWiMS.model(formula=baseformula,data=TrainSet,type=type,testType=testType,pvalue=CVselection.pValue,variableList=variableList,size=size,loops=loops,elimination.bootstrap.steps=elimination.bootstrap.steps,fraction=fraction,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize);

		CurModel_S <- BSWiMS.models$forward.model;
		UCurModel_S <- BSWiMS.models$update.model;
		redCurmodel_S <- BSWiMS.models$BSWiMS.model;


		if (length(CurModel_S$var.names)>0)
		{

			thebaglist <- CurModel_S$formula.list
			if (length(thebaglist)>0)
			{
				bagg <- baggedModel(thebaglist,ETraningSet,type,Outcome,timeOutcome,univariate=variableList,useFreq=loops); 
			}
			else
			{
				bagg <- list(bagged.model= redCurmodel_S$back.model,formula=redCurmodel_S$back.formula);
			}
			baggedfoldmodel <- BSWiMS.models$bagging$bagged.model;
			if (is.null(baggedfoldmodel))
			{
				baggedfoldmodel <- redCurmodel_S$back.model;
			}
			
			Full_model <- modelFitting(Full_redCurmodel_S$back.formula,ETraningSet,type,fitFRESA=TRUE)
			redfoldmodel.AtOpt <- redCurmodel_S$at.opt.model;
			redfoldmodel <- redCurmodel_S$back.model;
			
			cat ("\nUpdate   :",UCurModel_S$formula,"\n")
			cat ("At RMSE  :",redCurmodel_S$at.RMSE.formula,"\n")
			cat ("B:SWiMS  :",redCurmodel_S$back.formula,"\n")
			cat ("BB:SWiMS  :",BSWiMS.models$bagging$formula,"\n")
			
			{
				inserted = inserted +1;

				predictTest <- predict.fitFRESA(baggedfoldmodel,BlindSet,"linear");
				predictTest.AtOpt <- predict.fitFRESA(redfoldmodel.AtOpt,BlindSet,"linear");
				predictTest.ForwardModel <- predict.fitFRESA(UCurModel_S$final.model,BlindSet,"linear");
				Full_predictTest <- predict.fitFRESA(Full_model,BlindSet,"linear");
				predictTrain <- predict.fitFRESA(baggedfoldmodel,TrainSet,"linear");
				Full_predictTrain <- predict.fitFRESA(Full_model,TrainSet,"linear");
				
				baggedForwardPredict <- predict.fitFRESA(bagg$bagged.model,BlindSet,"linear");
				
				
				medianPred <- ensemblePredict(BSWiMS.models$forward.selection.list,ETraningSet,BlindSet, predictType = "linear",type = type)$ensemblePredict
				eq=NULL;
				if ((length(redfoldmodel$coefficients) > 1) && equivalent)
				{
					collectFormulas <- BSWiMS.models$forward.selection.list;
					bagg2 <- baggedModel(collectFormulas,data,type,Outcome,timeOutcome,univariate=variableList,useFreq=loops);
					shortcan <- bagg2$frequencyTable;
					modelterms <- attr(terms(redfoldmodel),"term.labels");
					eshortlist <- unique(c(names(shortcan),str_replace_all(modelterms,":","\\*")));
					eshortlist <- eshortlist[!is.na(eshortlist)];
					if (length(eshortlist)>0)
					{
						nameslist <- c(all.vars(baggedfoldmodel$formula),as.character(variableList[eshortlist,2]));
						nameslist <- unique(nameslist[!is.na(nameslist)]);
						if (!is.null(unirank) && (unirank$categorizationType != "RawRaw"))
						{
							eqdata <- ETraningSet[,nameslist];
						}
						else
						{
							eqdata <- ETraningSet;
						}
						eq <- reportEquivalentVariables(redfoldmodel,pvalue = 0.25*pvalue,
									  data=eqdata,
									  variableList=cbind(eshortlist,eshortlist),
									  Outcome = Outcome,
									  timeOutcome=timeOutcome,
									  type = type,osize=featureSize,
									  method="BH",fitFRESA=TRUE);
									  
						eqpredict <- ensemblePredict(eq$formula.list,ETraningSet,BlindSet, predictType = "linear",type = type)$ensemblePredict;
						equiFormulas <- append(equiFormulas,eq$formula.list);
					}
					else
					{
						eqpredict <- predictTest;
					}
				}
				else
				{
					eqpredict <- predictTest;
				}
				
				palt <- NULL;
				if (!is.null(usrFitFun))
				{
					fit <- usrFitFun(formula(abaseformula),TrainSet[,c(extvar,shortVarList)]);
					palt <- predict(fit,BlindSet);
					vset <- all.vars(baggedfoldmodel$formula);
					if (!is.null(eq))
					{
						vset <- unique(append(vset,all.vars(eq$equivalentModel$formula)));
					}
					if (length(vset) == (1+1*(type=="COX")))
					{
						vset <- all.vars(UCurModel_S$final.model$formula);
					}
					if (length(vset) > (1+1*(type=="COX")))
					{
						fit <- usrFitFun(formula(abaseformula),TrainSet[,vset]);
						palt <- cbind(palt,predict(fit,BlindSet));
					}
					else
					{
						palt <- cbind(palt,numeric(nrow(BlindSet)));
					}
				}
				

				trainResiduals <- residualForFRESA(baggedfoldmodel,TrainSet,Outcome);
				blindResiduals <- residualForFRESA(baggedfoldmodel,BlindSet,Outcome);
				FulltrainResiduals <- residualForFRESA(Full_model,TrainSet,Outcome);
				FullblindResiduals <- residualForFRESA(Full_model,BlindSet,Outcome);


				px <- cbind(BlindSet[,Outcome],predictTest,i,blindResiduals,medianPred,baggedForwardPredict,predictTest.ForwardModel,predictTest.AtOpt,eqpredict);
				if (!is.null(usrFitFun)) {px <- cbind(px,palt);}
				rownames(px) <- blindsampleidx;
				totSamples <- rbind(totSamples,px);
				px <- cbind(BlindSet[,Outcome],Full_predictTest,i,FullblindResiduals);
				rownames(px) <- blindsampleidx;
				Full.totSamples <- rbind(Full.totSamples,px);
				px <- cbind(TrainSet[,Outcome],predictTrain,i,trainResiduals);
				rownames(px) <- sampleidx;
				totTrainSamples <- rbind(totTrainSamples,px);
				px <- cbind(TrainSet[,Outcome],Full_predictTrain,i,FulltrainResiduals);
				rownames(px) <- sampleidx;
				Full.totTrainSamples <- rbind(Full.totTrainSamples,px);

				formulas <- append(formulas,BSWiMS.models$bagging$formula);
				AtOptFormulas <- append(AtOptFormulas,redCurmodel_S$at.RMSE.formula);
				ForwardFormulas <- append(ForwardFormulas,UCurModel_S$formula);
				baggFormulas <- append(baggFormulas,bagg$formula);
				allBSWIMSFormulas <- append(allBSWIMSFormulas,BSWiMS.models$formula.list);

				
				trainRMS <- sqrt(sum(trainResiduals^2)/nrow(TrainSet));
				trainPearson <- cor.test(TrainSet[,Outcome], predictTrain, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				trainSpearman <- cor.test(TrainSet[,Outcome], predictTrain, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

				FulltrainRMS <- sqrt(sum(FulltrainResiduals^2)/nrow(TrainSet));
				FulltrainPearson <- cor.test(TrainSet[,Outcome], Full_predictTrain, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				FulltrainSpearman <- cor.test(TrainSet[,Outcome], Full_predictTrain, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

				blindRMS <- sqrt(sum(blindResiduals^2)/nrow(BlindSet));
				FullblindRMS <- sqrt(sum(FullblindResiduals^2)/nrow(BlindSet));

				if (nrow(BlindSet)>5)
				{
					foldPearson <- cor.test(BlindSet[,Outcome], predictTest, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					foldSpearman <- cor.test(BlindSet[,Outcome], predictTest, method = "spearman",na.action=na.omit,exact=FALSE)$estimate				
					cstat <- rcorr.cens(predictTest,BlindSet[,Outcome], outx=FALSE)[1];
					foldRMS <- sum(blindResiduals^2)/(nrow(BlindSet)-1);
					cat("Fold RMS: ",sqrt(foldRMS),"Fold Test Pearson: ", foldPearson, "Fold Test Spearman: ",foldSpearman,"Fold Cstat:",cstat,"\n");

					blindFoldMS <- append(blindFoldMS,foldRMS);
					blindFoldPearson <- append(blindFoldPearson,foldPearson);
					blindFoldSpearman <- append(blindFoldSpearman,foldSpearman);
					blindFoldCstat <- append(blindFoldCstat,cstat);
	#				cat("Accu RMS: ",sqrt(mean(blindFoldMS)),"Accu Test Pearson: ", mean(blindFoldPearson), "Accu Test Spearman: ",mean(blindFoldSpearman),"Accu Cstat:",mean(blindFoldCstat),"\n");
				}

	# univariate analysis of top model residuals
				uniEval <- getVar.Res(Full_UCurModel_S$final.model,TrainSet,Outcome,type = type,testData=BlindSet);
				if (i==1)
				{
					uniTrainMSS <- rbind(uniEval$unitrainMSE);
					uniTestMSS <- rbind(uniEval$unitestMSE);
				}
				else
				{
					uniTrainMSS <- rbind(uniTrainMSS,uniEval$unitrainMSE);
					uniTestMSS <- rbind(uniTestMSS,uniEval$unitestMSE);
				}


				if (nrow(totSamples)>5)
				{
					blindPearson <- cor.test(totSamples[,1], totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					blindSpearman <- cor.test(totSamples[,1], totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)
					blindForwardSpearman <- cor.test(totSamples[,1], totSamples[,6], method = "spearman",na.action=na.omit,exact=FALSE)$estimate

					cstat <- rcorr.cens(totSamples[,2],totSamples[,1], outx=FALSE)[1];

					FullblindPearson <- cor.test(Full.totSamples[,1], Full.totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
					FullblindSpearman <- cor.test(Full.totSamples[,1], Full.totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)

					AcumRMS <- sqrt(sum(totSamples[,4]^2)/nrow(totSamples));
					AcumFullRMS <- sqrt(sum(Full.totSamples[,4]^2)/nrow(totSamples));


					cat("Samples: ", nrow(totSamples),"Full RMS:",AcumFullRMS," Accumulated Blind RMS: ", AcumRMS," c-index : ",cstat,"\n");
					cat("Full Blind RMS: ", FullblindRMS, " Full Train RMS: ",FulltrainRMS,"\n");
					cat("Blind Pearson: ", blindPearson, " Train Pearson: ",trainPearson,"\n");
					cat("Blind Spearman: ", blindSpearman$estimate,"(", blindSpearman$p.value,")  Train Spearman: ",trainSpearman,"Forward Bagged Spearman:",blindForwardSpearman,"\n");
					cat("Full Blind Pearson: ", FullblindPearson , " Full Train Pearson: ",FulltrainPearson,"\n");
					cat("Full Blind Spearman: ", FullblindSpearman$estimate, "(",FullblindSpearman$p.value,") Full Train Spearman: ",FulltrainSpearman,"\n");

					

				}
				
				vblindRMS <- append(vblindRMS,blindRMS);
				FullvblindRMS <- append(FullvblindRMS,FullblindRMS);
				
				vtrainRMS <- append(vtrainRMS,trainRMS);
				vtrainSpearman <- append(vtrainSpearman,trainSpearman);
				vtrainPearson <- append(vtrainPearson,trainPearson);

				FullvtrainRMS <- append(FullvtrainRMS,FulltrainRMS);
				FullvtrainSpearman <- append(FullvtrainSpearman,FulltrainSpearman);
				FullvtrainPearson <- append(FullvtrainPearson,FulltrainPearson);
			}


		}
		if ( (i %% K) == 0)
		{
			foldtest <- totSamples[totSamples[,3]>lastmodel,];
			CVBlindPearson <- append(CVBlindPearson,cor.test(foldtest[,1], foldtest[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindSpearman <- append(CVBlindSpearman,cor.test(foldtest[,1], foldtest[,2], method = "spearman",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindRMS <- append(CVBlindRMS,sqrt(sum((foldtest[,1]-foldtest[,2])^2)/nrow(foldtest)));
			lastmodel = i;
		}
	}
			
#	print(LASSOVariables)

	if (!is.null(usrFitFun)) 
	{
		colnames(totSamples) <- c("Outcome","Prediction","Model","Residuals","Median","Bagged","Forward","Backwards","Equivalent","usrFitFunction","usrFitFunction_Sel");
	}
	else
	{
		colnames(totSamples) <- c("Outcome","Prediction","Model","Residuals","Median","Bagged","Forward","Backwards","Equivalent");
	}

	rnames <- rownames(totSamples);
#	print(rnames)
#	totSamples <- as.data.frame(totSamples);
#	rnames <- rownames(totSamples);
#	print(rnames)
#	rownames(totSamples) <- rnames;
#	cat("Here 0\n")

	colnames(Full.totSamples) <- c("Outcome","Prediction","Model","Residuals");
#	Full.totSamples <- as.data.frame(Full.totSamples);

	colnames(totTrainSamples) <- c("Outcome","Prediction","Model","Residuals");
#	totTrainSamples <- as.data.frame(totTrainSamples);

	colnames(Full.totTrainSamples) <- c("Outcome","Prediction","Model","Residuals");
#	Full.totTrainSamples <- as.data.frame(Full.totTrainSamples);

#	print(totSamples)
#	cat("Here 1\n")
	BSWiMS.ensemble.prediction <- NULL
	bsta <- boxplot(totSamples[,"Prediction"]~rownames(totSamples),plot=FALSE)
	sta <- cbind(bsta$stats[3,])
	rownames(sta) <- bsta$names
#	BSWiMS.ensemble.prediction <- cbind(data[,Outcome],sta[rownames(data),])
	BSWiMS.ensemble.prediction <- cbind(data[rownames(sta),Outcome],sta)
	colnames(BSWiMS.ensemble.prediction) <- c("Outcome","Prediction");
	BSWiMS.ensemble.prediction <- as.data.frame(BSWiMS.ensemble.prediction);
#	print(BSWiMS.ensemble.prediction);
#	cat("Here 2\n")

	if (!is.null(enetSamples))
	{
		colnames(enetSamples) <- c("Outcome","Prediction","Model");
#		enetSamples <- as.data.frame(enetSamples);
		colnames(enetTrainSamples) <- c("Outcome","Prediction","Model");
#		enetTrainSamples <- as.data.frame(enetTrainSamples);
	}

	if (!is.null(uniTrainMSS))
	{
		uniTrainMSS <- as.data.frame(uniTrainMSS);
		uniTestMSS <- as.data.frame(uniTestMSS);
		colnames(uniTrainMSS) <-  attr(terms(formula(Full_UCurModel_S$formula)),'term.labels');
		colnames(uniTestMSS) <-  attr(terms(formula(Full_UCurModel_S$formula)),'term.labels');
	}

	
	blindRMS <- sqrt(sum((totSamples[,"Residuals"])^2)/nrow(totSamples));
	blindPearson <- cor.test(totSamples[,"Outcome"], totSamples[,"Prediction"], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	blindSpearman <- cor.test(totSamples[,"Outcome"], totSamples[,"Prediction"], method = "spearman",na.action=na.omit,exact=FALSE)$estimate

#	cat("Here 3\n")

	FullblindRMS <- sqrt(sum((Full.totSamples[,"Residuals"])^2)/nrow(totSamples));
	FullblindPearson <- cor.test(Full.totSamples[,"Outcome"], Full.totSamples[,"Prediction"], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	FullblindSpearman <- cor.test(Full.totSamples[,"Outcome"], Full.totSamples[,"Prediction"], method = "spearman",na.action=na.omit,exact=FALSE)$estimate

#	cat("Here 4\n")

	cstat <- rcorr.cens(totSamples[,"Prediction"],totSamples[,"Outcome"], outx=FALSE)[1];

	cat("##### CV of Initial Model (Biased) ###### \n");
	cat("Full Blind RMS: ", FullblindRMS,"\n")
	cat("Full Blind Spearman: ", FullblindSpearman,"\n")
	cat("Full Blind Pearson: ", FullblindPearson,"\n")
	if (length(blindFoldMS)>0)
	{
		cat("##### By Fold Analysis ###### \n Samples: ",length(blindFoldMS),"\n");
		cat("Sampled RMS: ", sqrt(mean(blindFoldMS,na.rm=TRUE)),"\n")
		cat("Mean Blind Spearman: ", mean(blindFoldSpearman,na.rm=TRUE)," (",sd(blindFoldSpearman,na.rm=TRUE),")\n")
		cat("Mean Blind Pearson: ", mean(blindFoldPearson,na.rm=TRUE)," (",sd(blindFoldPearson,na.rm=TRUE),")\n")
		cat("Mean Blind cstat: ",mean(blindFoldCstat,na.rm=TRUE)," (",sd(blindFoldCstat,na.rm=TRUE),")\n")
	}
	cat("##### Full Set Coherence Analysis ###### \n");
	cat("Blind  RMS: ", blindRMS," c-index : ",cstat,"\n");
	cat("Blind  Spearman: ", blindSpearman,"\n")
	cat("Blind  Pearson: ", blindPearson,"\n")
	if (length(CVBlindPearson)>0)
	{
		cat("##### Models Coherence Analysis ###### \n Samples: ",length(CVBlindPearson),"\n");
		cat("Mean Blind Spearman: ", mean(CVBlindSpearman,na.rm=TRUE)," (",sd(CVBlindSpearman,na.rm=TRUE),")\n")
		cat("Mean Blind Pearson: ", mean(CVBlindPearson,na.rm=TRUE)," (",sd(CVBlindPearson,na.rm=TRUE),")\n")
	}

	
	
	result <- list(		formula.list=formulas, 
						Models.testPrediction=totSamples,
						FullBSWiMS.testPrediction=Full.totSamples,
						BSWiMS=Full_redCurmodel_S,
						forwardSelection=Full_CurModel_S,
						updatedforwardModel=Full_UCurModel_S,
						testRMSE = blindRMS,
						testPearson = blindPearson,
						testSpearman = blindSpearman,
						FullTestRMSE = FullblindRMS,
						FullTestPearson = FullblindPearson,
						FullTestSpearman = FullblindSpearman,
						trainRMSE = vtrainRMS,
						trainPearson = vtrainPearson,
						trainSpearman = vtrainSpearman,
						FullTrainRMS = FullvtrainRMS,
						FullTrainPearson = FullvtrainPearson,
						FullTrainSpearman = FullvtrainSpearman,
						byFoldTestMS = blindFoldMS,
						byFoldTestSpearman = blindFoldSpearman,
						byFoldTestPearson = blindFoldPearson,
						byFoldCstat = blindFoldCstat,
						testRMSEAtFold = vblindRMS,
						FullTestRMSEAtFold = FullvblindRMS,
						Fullenet=Fullenet,
						LASSO.testPredictions=enetSamples,
						LASSOVariables=LASSOVariables,
						CVBlindPearson=CVBlindPearson,
						CVBlindSpearman=CVBlindSpearman,
						CVBlindRMS=CVBlindRMS,
						Models.trainPrediction=totTrainSamples,
						FullBSWiMS.trainPrediction=Full.totTrainSamples,
						LASSO.trainPredictions=enetTrainSamples,
						uniTrainMSS = uniTrainMSS,
						uniTestMSS = uniTestMSS,
						BSWiMS.ensemble.prediction = BSWiMS.ensemble.prediction,
						AtOptFormulas.list = AtOptFormulas,
						ForwardFormulas.list = ForwardFormulas,
						baggFormulas.list = baggFormulas,
						equiFormulas.list = equiFormulas,
						allBSWiMSFormulas.list = allBSWIMSFormulas,
						LassoFilterVarList = enetshortVarList,
						BSWiMS.models=FULLBSWiMS.models

					);
	return (result)
}
