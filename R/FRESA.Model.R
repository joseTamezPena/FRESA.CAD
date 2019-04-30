FRESA.Model <-
function(formula,data,OptType=c("Binary","Residual"),pvalue=0.05,filter.p.value=0.10,loops=32,maxTrainModelSize=20,elimination.bootstrap.steps=100,bootstrap.steps=100,print=FALSE,plots=FALSE,CVfolds=1,repeats=1,nk=0,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail","Tail","RawRaw"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,var.description=NULL,testType=c("zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"),lambda="lambda.1se",equivalent=FALSE,bswimsCycles=20,usrFitFun=NULL)
{

	a = as.numeric(Sys.time());
	set.seed(a);

	categorizationType <- match.arg(categorizationType);
	cl <- match.call();
	
	cvObject <-  NULL;
	univariate <- NULL;
	eq=NULL;
	bagg=NULL;

	type = "LM";
	if (class(formula)=="character")
	{
		formula <- formula(formula);
	}
	if (class(formula)=="formula")
	{
		featureSize = ncol(data)-1;
		OptType <- match.arg(OptType)
		
		varlist <- attr(terms(formula),"variables")
		
		dependent <- as.character(varlist[[2]])
		
		timeOutcome = NA;
		Outcome = NA;
		
		type = "LM";
		if (length(dependent)==3)
		{
			type = "COX"
			timeOutcome = dependent[2];
			Outcome = dependent[3];
			dependentout = paste(dependent[1],"(",dependent[2],",",dependent[3],")");
		}
		else
		{
			Outcome = dependent[1];
			dependentout = Outcome;	
		}
		
		setIntersect <- attr(terms(formula),"intercept")
		if (setIntersect == 0)
		{
			covariates = "0";
		}
		else
		{
			covariates = "1";
		}
		
		termslist <- attr(terms(formula),"term.labels");
		acovariates <- covariates[1];
		if (length(termslist)>0)
		{
			for (i in 1:length(termslist))
			{
				covariates <- paste(covariates,"+",termslist[i]);
				acovariates <- append(acovariates,termslist[i]);
			}
		}

		startOffset = length(termslist);
		variables <- vector();
		descrip <- vector();
		pnames <- as.vector(colnames(data));
		for (i in 1:length(pnames))
		{
			detected = 0;
			if (length(termslist)>0)
			{
				for (j in 1:length(termslist))
				{
					if (termslist[j] == pnames[i]) detected = 1;
				}
			}
			if (Outcome == pnames[i]) detected = 1;
			if (!is.na(timeOutcome) )
			{
				if (timeOutcome == pnames[i]) detected = 1;
			}
			if (detected == 0)
			{
				variables <- append(variables,pnames[i]);
				if (!is.null(var.description))
				{
					descrip <- append(descrip,var.description[i]);
				}
			}
		}
		if (!is.null(var.description))
		{
			variables <- cbind(variables,descrip);
		}
		else
		{
			variables <- cbind(variables,variables);
		}
		
		colnames(variables) <- c("Var","Description");
	
		if (CVfolds>nrow(data))
		{
			cat("Setting to LOO CV\n");
			CVfolds=nrow(data);
		}
		
		trainFraction <- 1.0-1.0/CVfolds;
		trainRepetition <- repeats*CVfolds;
		fraction = 1.0000;					# will be working with 1.0000 fraction of the samples for bootstrap training 
		varMax = nrow(variables);
		baseModel <- paste(dependentout,"~",covariates);
		cvObject = NULL;
		reducedModel = NULL;
		bootstrappedModel = NULL;
		UpdatedModel = NULL;
		filter.z.value <- abs(qnorm(filter.p.value))
		cutpvalue <- 3.0*filter.p.value
		if (cutpvalue > 0.45) cutpvalue=0.45;
		selectionType = match.arg(testType);
		testType = match.arg(testType);
		theScores <- names(table(data[,Outcome]))
		if (((length(theScores)>2)||(min(data[,Outcome])<0))&&(OptType == "Binary"))
		{
			OptType = "Residual";
		}
		
		if (categorizationType=="RawRaw")
		{
			rownames(variables) <- variables[,1];
			unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType="Raw",type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Regression",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
			univariate <- unirank$orderframe;
			featureSize <- nrow(univariate);
			unitPvalues <- (1.0-pnorm(univariate$ZUni));
			names(unitPvalues) <- univariate$Name;
			adjPvalues <- p.adjust(unitPvalues,"BH");
			variables <- variables[names(adjPvalues[adjPvalues <= 2*filter.p.value]),];
		}

		if (OptType == "Binary")
		{


			if (length(dependent)==1)
			{
				type = "LOGIT";
			}

#			elimination.pValue <- pvalue; 	# To test if the variable is part of the model 
			unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="zIDI",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome);
			univariate <- unirank$orderframe;
			featureSize <- nrow(univariate);
			
			unitPvalues <- (1.0-pnorm(univariate$ZUni));
			names(unitPvalues) <- univariate$Name;
			adjPvalues <- p.adjust(unitPvalues,"BH");
			varMax <- sum(univariate$ZUni >= filter.z.value);
			if (categorizationType == "Raw")
			{
				gadjPvalues <- adjPvalues[adjPvalues < 2*filter.p.value]			
				noncornames <- correlated_Remove(data,names(gadjPvalues),thr=0.99);
				if (length(noncornames) > 1) featureSize <- featureSize*length(noncornames)/length(gadjPvalues);
#				cat(length(noncornames),":",length(gadjPvalues),":",length(noncornames)/length(gadjPvalues),"\n");
			}
			pvarMax <- sum(adjPvalues < 2*filter.p.value);
			sizeM <- min(c(pvarMax,varMax));
			if (sizeM < 5) sizeM = min(c(5,nrow(univariate)));
			if (varMax >  nrow(univariate)) varMax = nrow(univariate);
			if (varMax < 5) varMax = min(c(5,nrow(univariate)));

			redlist <- adjPvalues < cutpvalue;
			totlist <- min(sum(1*redlist),100);
			cat("Unadjusted size:",sum(univariate$ZUni >= filter.z.value)," Adjusted Size:",pvarMax," Cut size:",sum(1*redlist),"\n")
			if (totlist<10)
			{
				redlist <- c(1:min(10,nrow(univariate)))
				totlist <- length(totlist);
			}
			
			
			cat("\n Z: ",filter.z.value,", Features to test: ",sizeM,",Adjust Size:",featureSize,"\n");
			shortUniv <- univariate[redlist,]

			
			if (CVfolds>1)
			{
				if (categorizationType!="RawRaw")  
				{
					rownames(variables) <- variables[,1];
#					unirank$variableList <- variables[unique(as.character(univariate[redlist,2])),]
				}
				cvObject <- crossValidationFeatureSelection_Bin(sizeM,fraction,c(pvalue,filter.p.value),loops,acovariates,Outcome,timeOutcome,NULL,data,maxTrainModelSize,type,selectionType,startOffset,elimination.bootstrap.steps,trainFraction,trainRepetition,bootstrap.steps,nk,unirank,print=print,plots=plots,lambda=lambda,equivalent=equivalent,bswimsCycles=bswimsCycles,usrFitFun,featureSize=featureSize);
				firstModel <- cvObject$forwardSelection;
				UpdatedModel <- cvObject$updateforwardSelection;
				reducedModel <- cvObject$BSWiMS;
				bootstrappedModel <- cvObject$FullBSWiMS.bootstrapped;
				BSWiMS.models <- cvObject$BSWiMS.models;
			}
			else
			{
				BSWiMS.models <- BSWiMS.model(formula=formula,data=data,type=type,testType=selectionType,pvalue=pvalue,variableList=shortUniv,size=sizeM,loops=loops,elimination.bootstrap.steps=bootstrap.steps,fraction=1.0,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize,NumberofRepeats=repeats);
				firstModel <- BSWiMS.models$forward.model;
				UpdatedModel <- BSWiMS.models$update.model;
				reducedModel <- BSWiMS.models$BSWiMS.model;
				bootstrappedModel <- reducedModel$bootCV;
			
			}

		}
		if (OptType == "Residual")
		{	
#			elimination.pValue <- pvalue; 	# To test if the variable is part of the model
			if (testType=="zIDI")
			{
				if ((testType=="zIDI")&&(length(theScores)>10))
				{
					warning("Switching to Regresion, More than 10 scores");
					testType = "Ftest";
				}
				else
				{
					cat("Doing a Ordinal Fit with zIDI Selection\n");
					cat("Ordinal Fit will be stored in BSWiMS.models$oridinalModels\n");
					cat("Use predict(BSWiMS.models$oridinalModels,testSet) to get the ordinal prediction on a new dataset \n");
				}
			}
			if (length(dependent)==1)
			{
				if ((length(theScores)>2)||(min(data[,Outcome])<0))
				{
					type = "LM";
					unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Regression",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
					if ((length(theScores)<=10)&&(testType=="zIDI"))
					{
						type = "LOGIT";
					}
				}
				else 
				{
					if (type == "LM") type = "LOGIT";
					unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
				}
			}
			else
			{
				    unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
			}
			univariate <- unirank$orderframe;
			featureSize <- nrow(univariate);
			unitPvalues <- (1.0-pnorm(univariate$ZUni));
			names(unitPvalues) <- univariate$Name;
			adjPvalues <- p.adjust(unitPvalues,"BH");
			varMax <- sum(univariate$ZUni >= filter.z.value);	
			if (categorizationType == "Raw")
			{
				gadjPvalues <- adjPvalues[adjPvalues < 2*filter.p.value]			
				noncornames <- correlated_Remove(data,names(gadjPvalues),thr=0.99);
				if (length(noncornames) > 1) featureSize <- featureSize*length(noncornames)/length(gadjPvalues);
#				cat(length(noncornames),":",length(gadjPvalues),":",length(noncornames)/length(gadjPvalues),"\n");
			}
			pvarMax <- sum(adjPvalues < 2*filter.p.value);
			sizeM <- min(c(pvarMax,varMax));
			if (sizeM < 5) sizeM = min(c(5,nrow(univariate)));
			if (varMax >  nrow(univariate)) varMax = nrow(univariate);
			if (varMax < 5) varMax = min(c(5,nrow(univariate)));

			bootstrappedModel = NULL;
			redlist <- adjPvalues < cutpvalue;
			totlist <- min(sum(1*redlist),100);
			

			cat("Features to test:",sizeM," Adjusted Size:",featureSize,"\n");
			if (totlist<10)
			{
				redlist <- c(1:min(10,nrow(univariate)))
				totlist <- length(totlist);
			}

			cat("\n Z: ",filter.z.value," Var Max: ",featureSize,"FitType: ",type," Test Type: ",testType,"\n");
			shortUniv <- univariate[redlist,]

			if (CVfolds>1)
			{
				if (categorizationType != "RawRaw")  
				{
					rownames(variables) <- variables[,1];
#					unirank$variableList <- variables[unique(as.character(univariate[redlist,2])),]
				}
				cvObject <- crossValidationFeatureSelection_Res(size=sizeM,fraction=fraction,pvalue=c(pvalue,filter.p.value),loops=loops,covariates=acovariates,Outcome=Outcome,timeOutcome=timeOutcome,variableList=unirank$variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,startOffset=startOffset,elimination.bootstrap.steps=elimination.bootstrap.steps,trainFraction=trainFraction,trainRepetition=trainRepetition,setIntersect=setIntersect,unirank=unirank,print=print,plots=plots,lambda=lambda,equivalent=equivalent,bswimsCycles=bswimsCycles,usrFitFun=usrFitFun,featureSize=featureSize);
				firstModel <- cvObject$forwardSelection;
				UpdatedModel <- cvObject$updatedforwardModel;
				reducedModel <- cvObject$BSWiMS;
				bootstrappedModel <- cvObject$BSWiMS$bootCV;
				BSWiMS.models <- cvObject$BSWiMS.models;
			}
			else
			{
				BSWiMS.models <- BSWiMS.model(formula=formula,data=data,type=type,testType=testType,pvalue=pvalue,variableList=shortUniv,size=sizeM,loops=loops,elimination.bootstrap.steps=bootstrap.steps,fraction=1.0,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots,featureSize=featureSize,NumberofRepeats=repeats);
				firstModel <- BSWiMS.models$forward.model;
				UpdatedModel <- BSWiMS.models$update.model;
				reducedModel <- BSWiMS.models$BSWiMS.model;
				bootstrappedModel <- reducedModel$bootCV;

			}
		}
	}
	else
	{
		cat("Expecting a formula object\n");
	}
	if (is.null(reducedModel))
	{
		result <- list(BSWiMS.model = NULL,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		forwardModel=firstModel,
		updatedforwardModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax,
#		independenSize=adjsize,
		call=cl);
	}
	else
	{
	
		eq <- NULL;		
		bagg <- NULL;
		if ((length(reducedModel$back.model$coefficients) > 1 ) && equivalent) 
		{
			collectFormulas <- BSWiMS.models$forward.selection.list;
			bagg <- baggedModel(collectFormulas,data,type,Outcome,timeOutcome,univariate=univariate,useFreq=loops);
			shortcan <- bagg$frequencyTable[(bagg$frequencyTable >= (loops*0.05))];
			modeltems <- attr(terms(reducedModel$back.model),"term.labels");
			eshortlist <- unique(c(names(shortcan),str_replace_all(modeltems,":","\\*")));
			eshortlist <- eshortlist[!is.na(eshortlist)];
			if (length(eshortlist)>0)
			{
				nameslist <- c(all.vars(BSWiMS.models$bagging$bagged.model$formula),as.character(univariate[eshortlist,2]));
				nameslist <- unique(nameslist[!is.na(nameslist)]);
				if (categorizationType != "RawRaw") 
				{
					eqdata <- data[,nameslist];
				}
				else
				{
					eqdata <- data;
				}
				eq <- reportEquivalentVariables(reducedModel$back.model,pvalue = 0.25*pvalue,
							  data=eqdata,
							  variableList=cbind(eshortlist,eshortlist),
							  Outcome = Outcome,
							  timeOutcome=timeOutcome,								  
							  type = type,osize=featureSize,
							  method="BH");
			}
		}

		result <- list(BSWiMS.model = BSWiMS.models$bagging$bagged.model,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		forwardModel=firstModel,
		updatedforwardModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax,
		bagging=bagg,
		eBSWiMS.model=eq,
		BSWiMS.models=BSWiMS.models,
		call=cl
		);
	}
	return (result);
}


