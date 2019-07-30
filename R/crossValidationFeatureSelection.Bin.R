crossValidationFeatureSelection_Bin <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),selectionType=c("zIDI","zNRI"),startOffset=0,elimination.bootstrap.steps=100,trainFraction=0.67,trainRepetition=9,bootstrap.steps=100,nk=0,unirank=NULL,print=TRUE,plots=TRUE,lambda="lambda.1se",equivalent=FALSE,bswimsCycles=10,usrFitFun=NULL,featureSize=0)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 

	nlenght <- ncol(data)-1;

	enetSamples <- NULL;
	enetTrainSamples <- NULL;

	casesample = subset(data,get(Outcome)  == 1);
	controlsample = subset(data,get(Outcome) == 0);

	casesamplesize <- nrow(casesample);

	controlsamplesize <- nrow(controlsample);
	
	filter.p.value = 2.0*pvalue;
	if (length(pvalue)>1)
	{
		filter.p.value = pvalue[2];
		pvalue=pvalue[1];
	}

	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);


	acc = 0.0;
	sen = 0.0;
	spe = 0.0;
	sizecases = 0;
	sizecontrol = 0;
	totsize = 0;
	paracc = 0;
	psen = 0;
	pspe = 0;

	Full.acc = 0.0;
	Full.sen = 0.0;
	Full.spe = 0.0;
	Full.paracc = 0;
	Full.psen = 0;
	Full.pspe = 0;



	formulas <- character();
	AtOptFormulas <-character();
	ForwardFormulas <- character();
	baggFormulas <- character();
	equiFormulas <- character();
	allBSWIMSFormulas <- character();
	trainCorrelations <- vector();
	trainAccuracy <- vector();
	trainSensitivity <- vector();
	trainSpecificity <- vector();
	trainAUC <- vector();
	testAccuracy <- vector();
	testSensitivity <- vector();
	testSpecificity <- vector();
	testAUC <- vector();

	blindCorrelations <- vector();
	WholeFoldBlindAccuracy <- vector();
	WholeFoldBlindSpecificity <- vector();
	WholeFoldBlindSensitivity <- vector();
	WholeFoldBlindAUC <- vector();
	FoldBlindAccuracy <- vector();
	FoldBlindSpecificity <- vector();
	FoldBlindSensitivity <- vector();
	TopUniCoherenceTest <- vector();
	selection.pValue <- pvalue;

	CVselection.pValue <- pvalue;

	par(mfrow=c(1,1))
	
	if (!is.null(unirank))
	{
		uprank <- update.uniRankVar(unirank,data=data,FullAnalysis=FALSE)
		variableList <- uprank$orderframe;
#		print(variableList[1:20,])
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
	Fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
	LASSOVariables <- NULL;
	if (inherits(Fullenet, "try-error"))
	{
		cat("enet Error")
		Fullenet <- NULL;
	}
	else
	{
		cenet <- as.matrix(coef(Fullenet,s=lambda))
		lanames <- names(cenet[as.vector(cenet[,1] != 0),])
		print(LASSOVariables <- paste(lanames[lanames !=  "(Intercept)" ],collapse=" + "))
	}
	
#	cat ("END LASSO\n")

	
	
	selType = selectionType;

	ex_covariates=covariates;
	
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
		extvar <- c(timeOutcome,Outcome);
	}
		
#	cat ("Data:",ncol(data),"Features:",nrow(variableList),"Size:",size,"\n")
	FULLBSWiMS.models <- BSWiMS.model(formula=baseformula,data=data,type=type,testType=selectionType,pvalue=pvalue,variableList=variableList,size=size,loops=loops,elimination.bootstrap.steps=elimination.bootstrap.steps,fraction=fraction,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize);
	CurModel_Full <- FULLBSWiMS.models$forward.model;
	UCurModel_Full <- FULLBSWiMS.models$update.model;
	redCurmodel_Full <- FULLBSWiMS.models$BSWiMS.model;
	
	Full_formula <- redCurmodel_Full$back.model;
	FullBootCross <- bootstrapValidation_Bin(1.0000,bootstrap.steps,redCurmodel_Full$back.formula,Outcome,data,type,plots=plots)
	redBootCross <- FullBootCross;

	cat ("CV pvalue    :",CVselection.pValue,"\n")
	cat ("Update    :",UCurModel_Full$formula,"\n")
	cat ("At Accuray:",redCurmodel_Full$at.Accuracy.formula,"\n")
	cat ("B:SWiMS   :",redCurmodel_Full$back.formula,"\n")

	if (is.null(FullBootCross))
	{
		stop("no initial model found\n");
	}
	if (print) summary(FullBootCross,2)

	
	inserted = 0;
	rocadded = 0;
	split.blindSen <- NULL;
	blindreboot <- NULL;
	KNNSamples <- NULL;
	Full.KNNSamples <- NULL;
	totSamples <- NULL;
	Full.totSamples <- NULL;
	totTrainSamples <- NULL;
	Full.totTrainSamples <- NULL;
	uniTrainAccuracy <- NULL;
	uniTestAccuracy <- NULL;
	totSamples <- NULL;

	Fullsammples <- min(casesamplesize,controlsamplesize);
	if ( K > Fullsammples) K=Fullsammples
#	cat("Number of folds: ",K,"\n");
	specificities <- c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05);


	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
			casefolds <- cvTools::cvFolds(casesamplesize, K,1,  "random");
			controlfolds <- cvTools::cvFolds(controlsamplesize, K,1,  "random");
			cycleinsert=0;
			totalUniCor=0;
		}

		CaseTrainSet <- casesample[casefolds$subsets[casefolds$which != j,],];
		CaseBlindSet <- casesample[casefolds$subsets[casefolds$which == j,],];
		ControlTrainSet <- controlsample[controlfolds$subsets[controlfolds$which != j,],];
		ControlBlindSet <- controlsample[controlfolds$subsets[controlfolds$which == j,],];

		TrainSet <- rbind(CaseTrainSet,ControlTrainSet);
		BlindSet <- rbind(CaseBlindSet,ControlBlindSet);
		framesize <- nrow(BlindSet);
		minTrainSamples <- min(nrow(CaseTrainSet),nrow(ControlTrainSet));

		if (nk==0)
		{
			nk = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
		}
			

		KnnTrainSet <- rbind(CaseTrainSet[sample(1:nrow(CaseTrainSet),minTrainSamples,replace=FALSE),],ControlTrainSet[sample(1:nrow(ControlTrainSet),minTrainSamples,replace=FALSE),])
		
		par(mfrow=c(1,1))

		redBootCross <- bootstrapValidation_Bin(1.0000,bootstrap.steps,redCurmodel_Full$back.formula,Outcome,TrainSet,type,plots=plots)

		Full.p <- predict.fitFRESA(redBootCross$boot.model,BlindSet, 'linear');
		
		Fullknnclass <- getKNNpredictionFromFormula(FULLBSWiMS.models$bagging$formula,KnnTrainSet,BlindSet,Outcome,nk)
		


		if (!is.null(unirank))
		{
#			cat("Unirank\n")
			uprank <- update.uniRankVar(unirank,data=TrainSet,FullAnalysis=FALSE)
			variableList <- uprank$orderframe;

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
			cat(ncol(TrainSet),": Unadjusted size:",varMax," Adjusted Size:",pvarMax,"\n")
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

			
#			print(variableList[,1]);
		}
		
		
		
		if (!is.null(Fullenet))
		{
#			cat(length(shortVarList)," :In elastic Net\n")
			foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="binomial"));
			cenet <- as.matrix(coef(foldenet,s=lambda))
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
#			print(LASSOVariables)
		}




		cat ("Loop :",i,"Input Cases =",sum(data[,Outcome] > 0 ),"Input Control =",sum(data[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Train Cases =",sum(TrainSet[,Outcome] > 0 ),"Train Control =",sum(TrainSet[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Blind Cases =",sum(BlindSet[,Outcome] > 0 ),"Blind Control =",sum(BlindSet[,Outcome] == 0),"\n")
		cat ("K   :",nk,"KNN T Cases =",sum(KnnTrainSet[,Outcome] > 0 ),"KNN T Control =",sum(KnnTrainSet[,Outcome] == 0),"\n")

		lastinserted = inserted;


		BSWiMS.models <- BSWiMS.model(formula=baseformula,data=TrainSet,type=type,testType=selectionType,pvalue=CVselection.pValue,variableList=variableList,size=size,loops=loops,elimination.bootstrap.steps=elimination.bootstrap.steps,fraction=fraction,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize);
		CurModel_S <- BSWiMS.models$forward.model;
		UCurModel_S <- BSWiMS.models$update.model;
		redCurmodel_S <- BSWiMS.models$BSWiMS.model;
		redBootCross_S <- redCurmodel_S$bootCV;


		if (length(CurModel_S$var.names)>0)
		{
			 

# lets use beforeFSC model for prediction 
			redfoldmodel.AtOpt <- redCurmodel_S$at.opt.model;
			forwardmodel <- UCurModel_S$final.model;
			baggedfoldmodel <- BSWiMS.models$bagging$bagged.model;
			redfoldmodel <- redCurmodel_S$back.model;



			cat ("Update     :",UCurModel_S$formula,"\n")
			cat ("At Accuracy:",redCurmodel_S$at.Accuracy.formula,"\n")
			cat ("B:SWiMS    :",redCurmodel_S$back.formula,"\n")

			if (!is.null(forwardmodel))
			{

				if (print) 
				{
					cat ("\n The last CV bootstrapped model")
					s <- summary(redBootCross_S,2)
				}
			
				thebaglist <- CurModel_S$formula.list;
				bagg <- baggedModel(thebaglist,TrainSet,type,Outcome,timeOutcome,univariate=variableList,useFreq=loops); 


				
			
#					cat ("Predict 1\n")
				
				p.AtOpt <- predict.fitFRESA(redfoldmodel.AtOpt,BlindSet, 'linear');
				p.forward <- predict.fitFRESA(forwardmodel,BlindSet, 'linear');
				baggedForwardPredict <- predict.fitFRESA(bagg$bagged.model,BlindSet,'linear');
				firstBSWIMSPredict <- predict.fitFRESA(redfoldmodel,BlindSet, 'linear');
				if (length(BSWiMS.models$formula.list)>1)
				{
					medianBSWIMSPredict <- ensemblePredict(BSWiMS.models$formula.list,TrainSet,BlindSet, predictType = "linear",type = type)$ensemblePredict;
				}
				else
				{
					medianBSWIMSPredict <- firstBSWIMSPredict;
				}
#					cat ("Predict 2\n")
				medianPred <- ensemblePredict(BSWiMS.models$forward.selection.list,TrainSet,BlindSet, predictType = "linear",type = type)$ensemblePredict;
				eq <- NULL;
#					cat ("End Predict\n")

				if (length(redfoldmodel$coefficients)> 1)
				{
					p <- predict.fitFRESA(baggedfoldmodel,BlindSet,'linear');
					if (equivalent)
					{
						collectFormulas <- BSWiMS.models$forward.selection.list;
						bagg2 <- baggedModel(collectFormulas,data,type,Outcome,timeOutcome,univariate=variableList,useFreq=loops);
						modelterms <- attr(terms(redfoldmodel),"term.labels");
						shortcan <- bagg2$frequencyTable;
						eshortlist <- unique(c(names(shortcan),str_replace_all(modelterms,":","\\*")));
						eshortlist <- eshortlist[!is.na(eshortlist)];
						if (length(eshortlist)>0)
						{
							nameslist <- c(all.vars(baggedfoldmodel$formula),as.character(variableList[eshortlist,2]));
							nameslist <- unique(nameslist[!is.na(nameslist)]);
							if (!is.null(unirank) && (unirank$categorizationType != "RawRaw"))
							{
								eqdata <- TrainSet[,nameslist];
							}
							else
							{
								eqdata <- TrainSet;
							}
							eq <- reportEquivalentVariables(redfoldmodel,pvalue = 0.25*pvalue,
										  data=eqdata,
										  variableList=cbind(eshortlist,eshortlist),
										  Outcome = Outcome,
										  timeOutcome=timeOutcome,								  
										  type = type,osize=featureSize,
										  method="BH",fitFRESA=TRUE);
										  
							eqpredict <- ensemblePredict(eq$formula.list,TrainSet,BlindSet, predictType = "linear",type = type)$ensemblePredict;
							equiFormulas <- append(equiFormulas,eq$formula.list);
						}
					}
					else
					{
						eqpredict <- p;
					}

				}
				else
				{
					p <- firstBSWIMSPredict;
					eqpredict <- firstBSWIMSPredict;
				}

#				cat("B KNN\n")
				knnclass <- getKNNpredictionFromFormula(baggedfoldmodel$formula,KnnTrainSet,BlindSet,Outcome,nk);
#				cat("A KNN \n")

				palt <- NULL;
				if (!is.null(usrFitFun))
				{
					fit <- usrFitFun(formula(abaseformula),TrainSet[,c(extvar,shortVarList)]);
					palt <- predict(fit,BlindSet);
					if (is.null(baggedfoldmodel))
					{
						vset <- all.vars(bagg$formula);
					}
					else
					{
						vset <- all.vars(baggedfoldmodel$formula);
					}
					if (!is.null(eq))
					{
						vset <- unique(append(vset,all.vars(eq$equivalentModel$formula)));
					}
					if (length(vset) == (1+1*(type=="COX")))
					{
						vset <- all.vars(forwardmodel$formula);
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

				
				inserted = inserted + 1
				cycleinsert = cycleinsert + 1

				tcor <- cor.test(predict.fitFRESA(redfoldmodel,TrainSet, 'linear'),predict.fitFRESA(redBootCross$boot.model,TrainSet, 'linear'), method = "spearman",na.action=na.omit,exact=FALSE)$estimate
				trainCorrelations <- append(trainCorrelations,tcor);
				trainAccuracy <- append(trainAccuracy,redBootCross_S$base.Accuracy);
				trainSensitivity <- append(trainSensitivity,redBootCross_S$base.Sensitivity);
				trainSpecificity <- append(trainSpecificity,redBootCross_S$base.Specificity);
				trainAUC <- append(trainAUC,mean(redBootCross_S$train.ROCAUC));


				bcor <- 0;
				if (framesize>5)
				{
					bcor <- cor.test(p, Full.p, method = "spearman",na.action=na.omit,exact=FALSE)$estimate;
					blindCorrelations <- append(blindCorrelations,bcor);
				}

				if (((sumca <- sum(BlindSet[,Outcome]>0)) > 1) && ((sumco <- sum(BlindSet[,Outcome]==0)) > 1))
				{

					atRoc <- pROC::roc(as.vector(BlindSet[,Outcome]), p,plot=FALSE,ci=TRUE,auc=TRUE,of='se',specificities=specificities,boot.n=100,smooth=FALSE,progress= 'none',quiet = TRUE)
					splitRoc <- atRoc$ci[,2];
					FullRocBlindAUC <- pROC::roc(as.vector(BlindSet[,Outcome]), Full.p,plot=FALSE,auc=TRUE,ci=FALSE,quiet = TRUE)$auc
					WholeFoldBlindAUC <- append(WholeFoldBlindAUC,FullRocBlindAUC);
					if (rocadded == 0)
					{
						split.blindSen <- splitRoc;
					}
					else
					{
						split.blindSen <- rbind(split.blindSen,splitRoc);
					}
					rocadded = rocadded + 1;
				}
				
				totsize <- totsize + framesize;
				scase <- sum(BlindSet[,Outcome] == 1);
				scotr <- sum(BlindSet[,Outcome] == 0);
				sizecases <- sizecases + scase;
				sizecontrol <- sizecontrol + scotr;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( p >= 0.0 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( p < 0.0 )) , na.rm = TRUE)
				acc <- acc + psen + pspe;
				sen <- sen + psen;
				spe <- spe + pspe;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( Full.p >= 0.0 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( Full.p < 0.0 )) , na.rm = TRUE)
				Full.acc <- Full.acc + psen + pspe;
				Full.sen <- Full.sen + psen;
				Full.spe <- Full.spe + pspe;
				paracc = acc/totsize;
				psen = 0;
				pspe = 0;
				if (sizecases>0) 
				{
					psen = sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					pspe = spe/sizecontrol;
				}

				Full.paracc = Full.acc/totsize;
				Full.psen = 0;
				Full.pspe = 0;
				if (sizecases>0) 
				{
					Full.psen = Full.sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					Full.pspe = Full.spe/sizecontrol;
				}

				WholeFoldBlindAccuracy <- append(WholeFoldBlindAccuracy,redBootCross$blind.accuracy);
				WholeFoldBlindSpecificity <- append(WholeFoldBlindSpecificity,redBootCross$blind.specificity);
				WholeFoldBlindSensitivity <- append(WholeFoldBlindSensitivity,redBootCross$blind.sensitivity);

				FoldBlindAccuracy <- append(FoldBlindAccuracy,redBootCross_S$blind.accuracy);
				FoldBlindSpecificity <- append(FoldBlindSpecificity,redBootCross_S$blind.specificty);
				FoldBlindSensitivity <- append(FoldBlindSensitivity,redBootCross_S$blind.sensitivity);

				Full.ptrain <- predict.fitFRESA(redBootCross$boot.model,TrainSet, 'linear');
				ptrain <- predict.fitFRESA(redfoldmodel,TrainSet, 'linear');

				if ( cycleinsert == 1)
				{
					cvcycle.predictions <- cbind(BlindSet[,Outcome],p.AtOpt,i);
				}

				px <- cbind(BlindSet[,Outcome],p,i,medianPred,baggedForwardPredict,p.forward,p.AtOpt,eqpredict,firstBSWIMSPredict,medianBSWIMSPredict);
				if (!is.null(usrFitFun)) {px <- cbind(px,palt);}
				rownames(px) <- rownames(BlindSet);
				totSamples <- rbind(totSamples,px);
				px <- cbind(BlindSet[,Outcome],p.AtOpt,i);
				rownames(px) <- rownames(BlindSet);
				cvcycle.predictions <- rbind(cvcycle.predictions,px);
				px <- cbind(BlindSet[,Outcome],Full.p,i);
				rownames(px) <- rownames(BlindSet);
				Full.totSamples <- rbind(Full.totSamples,px);
				px <- cbind(BlindSet[,Outcome],abs(knnclass$prob$prob-1*(knnclass$prediction=="0")),i);
				rownames(px) <- rownames(BlindSet);
				KNNSamples <- rbind(KNNSamples,px);
				px <- cbind(BlindSet[,Outcome],abs(Fullknnclass$prob$prob-1*(Fullknnclass$prediction=="0")),i);
				rownames(px) <- rownames(BlindSet);
				Full.KNNSamples <- rbind(Full.KNNSamples,px);
				px <- cbind(TrainSet[,Outcome],ptrain,i);
				rownames(px) <- rownames(TrainSet);
				totTrainSamples <- rbind(totTrainSamples,px);
				px <- cbind(TrainSet[,Outcome],Full.ptrain,i);
				rownames(px) <- rownames(TrainSet);
				Full.totTrainSamples <- rbind(Full.totTrainSamples,px);
				
				
				formulas <- append(formulas,BSWiMS.models$bagging$formula);
				AtOptFormulas <- append(AtOptFormulas,redCurmodel_S$at.Accuracy.formula);
				ForwardFormulas <- append(ForwardFormulas,UCurModel_S$formula);
				baggFormulas <- append(baggFormulas,bagg$formula);
				allBSWIMSFormulas <- append(allBSWIMSFormulas,BSWiMS.models$formula.list);

				knnACC <- sum(KNNSamples[,1] == (KNNSamples[,2]>0.5))/totsize;
				knnSEN <- sum((KNNSamples[,1]>0.5) & (KNNSamples[,2]>0.5))/sizecases;
				knnSPE <- sum((KNNSamples[,1]<0.5) & (KNNSamples[,2]<0.5))/sizecontrol;

				Full.knnACC <- sum(Full.KNNSamples[,1] == (Full.KNNSamples[,2]>0.5))/totsize;
				Full.knnSEN <- sum((Full.KNNSamples[,1]>0.5) & (Full.KNNSamples[,2]>0.5))/sizecases;
				Full.knnSPE <- sum((Full.KNNSamples[,1]<0.5) & (Full.KNNSamples[,2]<0.5))/sizecontrol;


				cat ("Loop :",i,"Blind Cases =",scase,"Blind Control =",scotr,"Total =",totsize, "Size Cases =",sizecases,"Size Control =",sizecontrol,"\n")
				cat ("Accumulated Models CV Accuracy        =",paracc,"Sensitivity =",psen,"Specificity =",pspe,"Forw. Ensemble Accuracy=",mean(1.0*((totSamples[,4] > 0) == (totSamples[,1] > 0))),"\n")
				cat ("Initial Model Accumulated CV Accuracy =",Full.paracc,"Sensitivity =",Full.psen,"Specificity =",Full.pspe,"\n");
				cat ("Initial Model Bootstrapped Accuracy   =",redBootCross$blind.accuracy,"Sensitivity =",redBootCross$blind.sensitivity,"Specificity =",redBootCross$blind.specificity,"\n")
				cat ("Current Model Bootstrapped Accuracy   =",redBootCross_S$blind.accuracy,"Sensitivity =",redBootCross_S$blind.sensitivity,"Specificity =",redBootCross_S$blind.specificity,"\n")
				cat ("Current KNN Accuracy   =",knnACC,"Sensitivity =",knnSEN,"Specificity =",knnSPE,"\n")
				cat ("Initial KNN Accuracy   =",Full.knnACC,"Sensitivity =",Full.knnSEN,"Specificity =",Full.knnSPE,"\n")
				cat ("Train Correlation: ",tcor," Blind Correlation :",bcor,"\n KNN to Model Confusion Matrix: \n")
				print(table(KNNSamples[,2]>0.5,totSamples[,2]>0.0))
			}
			
		}
		else
		{
			cat ("Loop :",i,"No Model.\n")
		}

		uniEval <- getVar.Bin(UCurModel_Full$final.model,TrainSet,Outcome,type = type,testData=BlindSet);
		if (i==1)
		{
			uniTrainAccuracy <- rbind(uniEval$uniTrainAccuracy);
			TopUniTrainCor <- vector();
		}
		else
		{
			uniTrainAccuracy <- rbind(uniTrainAccuracy,uniEval$uniTrainAccuracy);
		}
		if ( j == 1)
		{
			cvcycle.uniAccuracies <- uniEval$uniTestAccuracy * framesize;
			totblindadded = framesize;
			topUniTestCor <- vector();
			totalUniCor = 0; 
		}
		else
		{
			cvcycle.uniAccuracies <- rbind(cvcycle.uniAccuracies,uniEval$uniTestAccuracy * framesize);
			totblindadded = totblindadded + framesize;
		}

		if ((lastinserted<inserted)&&(length(redCurmodel_S$back.model$coefficients)>1))
		{		
			uniEvalCor <- getVar.Bin(redCurmodel_S$back.model,TrainSet,Outcome,type = type,testData=BlindSet);
			TopUniTrainCor <- append(TopUniTrainCor,uniEvalCor$uniTrainAccuracy[1]);
			topUniTestCor <- append(topUniTestCor,uniEvalCor$uniTestAccuracy[1] * framesize);
			totalUniCor <- totalUniCor + framesize
		}
		
				
		if ( j == K)
		{
			if (totalUniCor>0) TopUniCoherenceTest <- append(TopUniCoherenceTest,sum(topUniTestCor)/totalUniCor)
			if (i == K)
			{
				uniTestAccuracy <- rbind(colSums(cvcycle.uniAccuracies)/totblindadded);
			}
			else
			{
				uniTestAccuracy <- rbind(uniTestAccuracy,colSums(cvcycle.uniAccuracies)/totblindadded);
			}
		}



		if ( j == K)
		{
			nsamp <- nrow(cvcycle.predictions)
			if (nsamp>0)
			{
				atRocAUC <- pROC::roc(as.vector(cvcycle.predictions[,1]), cvcycle.predictions[,2],plot=FALSE,auc=TRUE,smooth=FALSE,quiet = TRUE)$auc;
				testAccuracy <- append(testAccuracy,sum(cvcycle.predictions[,1] == 1.0*(cvcycle.predictions[,2]>=0.0))/nsamp);
				testSensitivity <- append(testSensitivity,sum((cvcycle.predictions[,1] == 1) & (cvcycle.predictions[,2]>=0.0))/sum(cvcycle.predictions[,1] == 1));
				testSpecificity <- append(testSpecificity,sum((cvcycle.predictions[,1] == 0) & (cvcycle.predictions[,2] <0.0))/sum(cvcycle.predictions[,1] == 0));
				testAUC <- append(testAUC,atRocAUC);
			}
#			print(testAccuracy)
#			print(testAUC)
		}
		
	}
	if (length(formulas)==0)
	{
		stop("No Significant Models Found\n");
	}
	if (!is.null(usrFitFun)) 
	{
		colnames(totSamples) <- c("Outcome","Prediction","Model","Ensemble.Forward","Forward.Selection.Bagged","Forward","Backwards","eB.SWiMS","first.B.SWiMS","Ensemble.B.SWiMS","usrFitFunction","usrFitFunction_Sel");
	}
	else
	{
		colnames(totSamples) <- c("Outcome","Prediction","Model","Ensemble.Forward","Forward.Selection.Bagged","Forward","Backwards","eB.SWiMS","first.B.SWiMS","Ensemble.B.SWiMS");
	}
#	totSamples <- as.data.frame(totSamples);

	colnames(Full.totSamples) <- c("Outcome","Prediction","Model");
#	Full.totSamples <- as.data.frame(Full.totSamples);

	colnames(totTrainSamples) <- c("Outcome","Prediction","Model");
#	totTrainSamples <- as.data.frame(totTrainSamples);
	colnames(Full.totTrainSamples) <- c("Outcome","Prediction","Model");
#	Full.totTrainSamples <- as.data.frame(Full.totTrainSamples);

	colnames(KNNSamples) <- c("Outcome","Prediction","Model");
#	KNNSamples <- as.data.frame(KNNSamples);
	
	colnames(Full.KNNSamples) <- c("Outcome","Prediction","Model");
#	Full.KNNSamples <- as.data.frame(Full.KNNSamples);

	
	
	BSWiMS.ensemble.prediction <- NULL

	bsta <- boxplot(totSamples[,"Prediction"]~rownames(totSamples),plot=FALSE)
	sta <- cbind(bsta$stats[3,])
	rownames(sta) <- bsta$names
	BSWiMS.ensemble.prediction <- cbind(data[rownames(sta),Outcome],sta)
	colnames(BSWiMS.ensemble.prediction) <- c("Outcome","Prediction");
	BSWiMS.ensemble.prediction <- as.data.frame(BSWiMS.ensemble.prediction);
	
	if (!is.null(enetSamples))
	{
		colnames(enetSamples) <- c("Outcome","Prediction","Model");
#		enetSamples <- as.data.frame(enetSamples);
		colnames(enetTrainSamples) <- c("Outcome","Prediction","Model");
#		enetTrainSamples <- as.data.frame(enetTrainSamples);
	}

	sumSen = NA;
	
	if (plots)
	{
		plotModels.ROC(totSamples,theCVfolds=K,predictor="Prediction",main="B:SWiMS");
		par(mfrow=c(1,1))
		incBsen=0
		aucBlindTest <- pROC::roc(as.vector(totSamples[,1]),totSamples[,2],col="red",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3,quiet = TRUE)$auc
		par(new=TRUE)
		aucCVBlind <- pROC::roc(as.vector(Full.totSamples[,1]),Full.totSamples[,2],col="blue",auc=TRUE,plot=TRUE,ci=FALSE,smooth=FALSE,quiet = TRUE)$auc
		aucBoot=0;
		aucTrain=0;
		if (!is.null(FullBootCross$testPrediction))
		{
			par(new=TRUE)
			aucTrain <- pROC::roc( as.vector(FullBootCross$outcome), FullBootCross$boot.model$linear.predictors,col="green",plot=TRUE,auc=TRUE,smooth=FALSE,quiet = TRUE)$auc;        
			par(new=TRUE)
			aucBoot <- pROC::roc( as.vector(FullBootCross$testOutcome), FullBootCross$testPrediction,col="black",auc=TRUE,plot=TRUE,smooth=FALSE,quiet = TRUE)$auc;
		}
		ley.names <- c(paste("Bootstrapped: Train Model ROC (",sprintf("%.3f",aucTrain),")"),paste("Bootstrapped: Blind ROC (",sprintf("%.3f",aucBoot),")"),
		paste("CV: Blind ROC (",sprintf("%.3f",aucCVBlind),")"),paste("CV: Blind Fold Models Coherence (",sprintf("%.3f",aucBlindTest),")"))
		ley.colors <- c("green","black","blue","red")
		ley.lty <- c(1,1,1,3)
		if (rocadded>0)
		{
			boxplot(split.blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
			sumSen <- colMeans(split.blindSen,na.rm = TRUE);
			sennames <- names(sumSen);
			sumSen <- append(0,sumSen);
			sumSen <- append(sumSen,1);
			sennames <- append("1",sennames);
			sennames <- append(sennames,"0");
			names(sumSen) <- sennames;
			spevalues <- as.numeric(names(sumSen));
			lines(spevalues,sumSen,col="red",lwd=2.0);
			auc = 0;
			for (i in 2:length(spevalues))
			{
				auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
			}
			ley.names <- append(ley.names,paste("CV Blind: Mean ROC of Models (",sprintf("%.3f",auc),")"));
			ley.colors <- append(ley.colors,"red");
			ley.lty  <- append(ley.lty,1);
		}
		else
		{
			sumSen = NA;
		}
		
		legend(0.6,0.30, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n")
	}


	if (!is.null(uniTrainAccuracy))
	{
		uniTrainAccuracy <- as.data.frame(uniTrainAccuracy);
		uniTestAccuracy <- as.data.frame(uniTestAccuracy);
		colnames(uniTrainAccuracy) <-  attr(terms(formula(UCurModel_Full$formula)),'term.labels');
		colnames(uniTestAccuracy) <-  attr(terms(formula(UCurModel_Full$formula)),'term.labels');
	}
	
	result <- list(formula.list=formulas,
	Models.testPrediction=totSamples,
	FullBSWiMS.testPrediction=Full.totSamples,
	TestRetrained.blindPredictions=blindreboot,
	LastTrainBSWiMS.bootstrapped=redCurmodel_S$bootCV,
	Test.accuracy=paracc,
	Test.sensitivity=psen,
	Test.specificity=pspe,
	Train.correlationsToFull=trainCorrelations,
	Blind.correlationsToFull=blindCorrelations,
	FullModelAtFoldAccuracies=WholeFoldBlindAccuracy,
	FullModelAtFoldSpecificties=WholeFoldBlindSpecificity,
	FullModelAtFoldSensitivities=WholeFoldBlindSensitivity,
	FullModelAtFoldAUC=WholeFoldBlindAUC,
	CVTrain.Accuracies=trainAccuracy,
	CVTrain.Sensitivity=trainSensitivity,
	CVTrain.Specificity=trainSpecificity,
	CVTrain.AUCs=trainAUC,
	CVTest.Accuracies=testAccuracy,
	CVTest.Sensitivity=testSensitivity,
	CVTest.Specificity=testSpecificity,
	CVTest.AUCs=testAUC,
	AtCVFoldModelBlindAccuracies=FoldBlindAccuracy,
	AtCVFoldModelBlindSpecificities=FoldBlindSpecificity,
	AtCVFoldModelBlindSensitivities=FoldBlindSensitivity,
	forwardSelection = CurModel_Full,
	updateforwardSelection = UCurModel_Full,
	BSWiMS = redCurmodel_Full,
	FullBSWiMS.bootstrapped=FullBootCross,
	Models.testSensitivities = split.blindSen,
	FullKNN.testPrediction=Full.KNNSamples,
	KNN.testPrediction=KNNSamples,
	Fullenet=Fullenet,
	LASSO.testPredictions=enetSamples,
	LASSOVariables=LASSOVariables,
	uniTrain.Accuracies=uniTrainAccuracy,
	uniTest.Accuracies=uniTestAccuracy,
	uniTest.TopCoherence=TopUniCoherenceTest,
	uniTrain.TopCoherence=TopUniTrainCor,
	Models.trainPrediction=totTrainSamples,
	FullBSWiMS.trainPrediction=Full.totTrainSamples,
	LASSO.trainPredictions=enetTrainSamples,
	BSWiMS.ensemble.prediction = BSWiMS.ensemble.prediction,
	ForwardFormulas.list = ForwardFormulas,
	AtOptFormulas.list = AtOptFormulas,
	baggFormulas.list = baggFormulas,
	equiFormulas.list = equiFormulas,
	allBSWiMSFormulas.list = allBSWIMSFormulas,
	LassoFilterVarList = enetshortVarList,
	BSWiMS.models=FULLBSWiMS.models
	);
	return (result)
}
