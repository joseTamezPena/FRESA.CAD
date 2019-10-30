BSWiMS.model <-
function(formula=formula,data=NULL,type=c("Auto","LM","LOGIT","COX"),testType=c("Auto","zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"),pvalue=0.05,variableList=NULL,size=0,loops=32,elimination.bootstrap.steps=100,fraction=1.0,maxTrainModelSize=20,maxCycles=20,print=FALSE,plots=FALSE,featureSize=0,NumberofRepeats=1)
{
	type <- match.arg(type);
	testType <- match.arg(testType);

	a = as.numeric(Sys.time());
	set.seed(a);
	
#	cat(featureSize," <- F Size\n");

	#	print(colnames(data))
#	print(variableList[,1])
    fforward.model <- NULL
	fupdate.model <- NULL;
	forward.selection.list <- NULL;
	if (featureSize==0) featureSize = ncol(data)-1;
	nfeat <- ncol(data)-1;
#	cat(size,":",pvalue,":",update.pvalue[1],":",elimination.pValue,": fs=",featureSize,"\n")
	if (class(formula)=="character")
	{
		formula <- str_replace_all(formula,"[.]","1");
		baseformula <- formula;
		formula <- formula(formula);
	}
	else
	{
		baseformula <- as.character(formula);
		baseformula[3] <- str_replace_all(baseformula[3],"[.]","1");
		baseformula <- paste(baseformula[2],"~",baseformula[3]);
		formula <- formula(baseformula);
	}
	varsmod <- all.vars(formula);
		
	varlist <- attr(terms(formula),"variables")
	termslist <- attr(terms(formula),"term.labels");
	setIntersect <- attr(terms(formula),"intercept");
	if (setIntersect == 0)
	{
		covariates = "0";
	}
	else
	{
		covariates = "1";
	}
	startOffset = length(termslist);
	acovariates <- covariates[1];
	if (length(termslist)>0)
	{
		for (i in 1:length(termslist))
		{
			covariates <- paste(covariates,"+",termslist[i]);
			acovariates <- append(acovariates,termslist[i]);
		}
	}
		
	dependent <- as.character(varlist[[2]])
	unitype ="Regression";
	rankingTest="Ztest"
	timeOutcome = NULL;
	Outcome = NULL;
	if (length(dependent)==3)
	{
		timeOutcome = dependent[2];
		Outcome = dependent[3];
		type = "COX";
		if (testType[1] == "Auto") testType="zIDI";
		unitype ="Binary";
		rankingTest="zIDI";
		univType = type;
	}
	else
	{
		Outcome = dependent[1];
		outcomeTable <- table(data[,Outcome]);
		theScores <- as.numeric(names(outcomeTable))
		totScores <- length(theScores);
		univType = type;
		if (type[1] == "Auto")
		{
			if (totScores>2)
			{
				type = "LM";
				univType = "LM";
				if (testType[1]=="Auto") 
				{
					testType="Ftest";
					if ((totScores <= 10) && (min(outcomeTable) >= 10))
					{
						type = "LOGIT";
						testType ="zIDI"
						univType = "LM";
						warning("Ordinal Model Fit\n")
					}
				}
			}
			else
			{
				univType = "LOGIT";
				if (min(data[,Outcome]) != 0)
				{
					type = "LM";
					if (testType[1]=="Auto") testType="Ftest";
				}
				else
				{
					type = "LOGIT";
					if (testType[1]=="Auto") testType="zIDI";
					unitype ="Binary";
					rankingTest="zIDI";
				}
			}
#			cat(testType[1],"<-");
		}
		if (testType[1]=="Auto") testType="Ftest";
	}
	theScores <- as.numeric(names(table(data[,Outcome])))
	totScores <- length(theScores);
#	print(theScores);

	unirank <- NULL;
	if (is.null(variableList))
	{
		vnames <- names(data);
		names(vnames) <- names(data);
		namesinformula <- vnames %in% all.vars(formula);
		vnames <- vnames[!namesinformula];
		variableList <- cbind(vnames,vnames)
		colnames(variableList) <- c("Name","Description");
		unirank <- uniRankVar(variableList,formula=baseformula,Outcome=Outcome,data=data,categorizationType="Raw",type=univType,rankingTest=rankingTest,uniType=unitype,FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
		variableList <- unirank$orderframe;
		rownames(variableList) <- variableList$Name;
#		print(variableList[1:10,]);
		unirank <- unirank$orderframe;
#		if (is.null(unitPvalues)) 
		{
			unitPvalues <- (1.0-pnorm(variableList$ZUni));
			names(unitPvalues) <- variableList$Name;
		}
		if (size==0)
		{
			featureSize <- max(featureSize,nrow(variableList));
			unitPvalues <- p.adjust(unitPvalues,"BH");
			unitPvalues <- unitPvalues[unitPvalues < 4*pvalue]; # 4 times the pvalue
#			print(unitPvalues);
			
#			filtered <- names(unitPvalues);
			filtered <- correlated_Remove(data,names(unitPvalues),thr=0.99);
			if (length(filtered) > 1) featureSize <- featureSize*length(filtered)/length(unitPvalues);

			size <- length(filtered);
#			cat ("Removed:", size,"\n");

			if (size<5)
			{
				size <- min(5,nrow(variableList));
			}
			else
			{
				variableList <- variableList[filtered ,]
#				print(variableList);
				if (is.null(timeOutcome))
				{
					data <- data[,c(Outcome,acovariates[-1],unique(as.character(variableList[,2])))];
				}
				else
				{
					data <- data[,c(Outcome,timeOutcome,acovariates[-1],unique(as.character(variableList[,2])))];
				}
			}
			if (print) cat(nrow(variableList),": Number of variables to test:",size,"\n");
		}
	}
	invariableList <- variableList;
	if (size<5)
	{
		size <- min(5,nrow(variableList));
	}
	if (print) cat(mean(data[,Outcome])," Repeats: ",NumberofRepeats,". Number of Features: ",nrow(variableList),": Number of features to test:",size,"\n");
	BSWiMS.model <- NULL;
	forward.model <- NULL;
	update.model <- NULL;
	formula.list <- character();
	forward.selection.list <- character();
	startVariableList <- variableList;
	startSize <- size;
	vartoTest <- variableList[1:size,]
	selectedVariableList <- NULL;
	theOutcome <- data[,Outcome];
	oridinalModels <- NULL
	if ((testType == "zIDI") && (totScores > 2) && (length(dependent) < 3))
	{
		oridinalModels <- list(theScores=theScores,data=data,formulas=NULL)
		class(oridinalModels) <- c("fitFRESA","ordinalFit");
	}
	halfSocres <- as.integer(totScores/2+0.5);
	IIRMetricPDF <- NULL;
	sdOutcome <- sd(theOutcome);
	infraction <- 0;
	cat("[");
	equivalent = FALSE;
	if (NumberofRepeats <= 0) 
	{
		equivalent = TRUE;
		NumberofRepeats = abs(NumberofRepeats) + 1*(NumberofRepeats == 0);
	}
	equiMaxFreq <- 0;
	addedEquFreq <- 0;
	for (nrep in 1:NumberofRepeats)
	{
		forward.selection.list <- character();
		firstModel <- NULL;
		metric <- 0;
		cycles <- 0;
		variableList <- startVariableList;
		size <- startSize;
		isInferior <- FALSE;
		while ( ( !isInferior || (infraction < 0.975) ) && (cycles<maxCycles) && (size>1))
		{
			isInferior <-  TRUE;
			infraction <- 1.0;
			ordinalFormulas <- NULL;
			if ((testType=="zIDI") || (testType=="zNRI") )
			{
				if (type == "LM")
				{
					type = "LOGIT"
				}
				if (!is.null(oridinalModels))
				{
					stdata <- data;
					sa <- theScores[length(theScores)];
					max.currentMeanAUC <- 0;
					for (s in theScores[-1])
					{
						stdata[,Outcome] <- 1*(data[,Outcome] >= s);
						oforward.model <- ForwardSelection.Model.Bin(size=size,fraction=fraction,pvalue,loops,acovariates,Outcome,variableList,stdata,maxTrainModelSize,type,timeOutcome,selectionType=testType,featureSize=featureSize);
						oupdate.model <- oforward.model$update.model;
						oBSWiMS.model <- bootstrapVarElimination_Bin(object=oupdate.model$final.model,pvalue=oforward.model$theZthr,Outcome=Outcome,data=stdata,startOffset=startOffset,type=type,selectionType=testType,loops=elimination.bootstrap.steps,print=print,plots=plots);
						if (length(all.vars(formula(oBSWiMS.model$back.formula)))>1)
						{
							ordinalFormulas <- append(ordinalFormulas,oBSWiMS.model$back.formula);
							currentMeanAUC <- (oBSWiMS.model$bootCV$blind.sensitivity + oBSWiMS.model$bootCV$blind.specificity)/2;
							if (currentMeanAUC > max.currentMeanAUC)
							{
								max.currentMeanAUC <- currentMeanAUC;
								BSWiMS.model <- oBSWiMS.model;
								forward.model <- oforward.model;
								update.model <- oforward.model$update.model;
								if (print)
								{
									cat("Score:",s," AUC :",currentMeanAUC,"\n")
								}
							}
						}
						else
						{
							ordinalFormulas <- append(ordinalFormulas,oupdate.model$formula);
							if (max.currentMeanAUC ==0) BSWiMS.model <- oBSWiMS.model;
						}
					}
				}
				else
				{
					forward.model <- ForwardSelection.Model.Bin(size=size,fraction=fraction,pvalue,loops,acovariates,Outcome,variableList,data,maxTrainModelSize,type,timeOutcome,selectionType=testType,featureSize=featureSize);
					update.model <- forward.model$update.model;
					BSWiMS.model <- bootstrapVarElimination_Bin(object=update.model$final.model,pvalue=forward.model$theZthr,Outcome=Outcome,data,startOffset=startOffset,type=type,selectionType=testType,loops=elimination.bootstrap.steps,print=print,plots=plots);
				}
				if (length(forward.model$var.names)>0)
				{
					if (elimination.bootstrap.steps>1)
					{
						if (!is.null(firstModel))
						{
							currentMeanAUC <- (BSWiMS.model$bootCV$blind.sensitivity + BSWiMS.model$bootCV$blind.specificity)/2;
							metric <- currentMeanAUC;
							if (length(attr(terms(formula(BSWiMS.model$back.formula)),"term.labels"))>0)
							{
								curAUC <-  (BSWiMS.model$bootCV$sensitivity + BSWiMS.model$bootCV$specificity)/2;
								firstMedAUC <- median(IIRMetricPDF);
								currentMedAUC <- median(curAUC);
								firstCount <- sum(currentMedAUC >= IIRMetricPDF)
								curCount <- sum(firstMedAUC <= curAUC)
								supchance <- sum(0.5 <= curAUC)/length(curAUC);
								infraction <- 1.0 - 0.5*firstCount/length(IIRMetricPDF)-0.5*curCount/length(curAUC);
								simTest <- ks.test(IIRMetricPDF,curAUC + rnorm(length(curAUC),0,1e-10))$p.value
								if (print) 
								{
#									hist(IIRMetricPDF)
#									hist(curAUC)
									cat(BSWiMS.model$back.formula,": Base AUC: ",firstMedAUC,"Current Blind AUC: ",currentMedAUC," Inferior Count:",firstCount," Tests:",length(IIRMetricPDF)," Fraction:",infraction," KStest:",simTest,"\n");
								}
								if (supchance < 0.75) infraction <- 1.0;
								isInferior <- (infraction > 0.975);
								if ( !isInferior && (cycles < 3) && (simTest > 0.05) )
								{
									IIRMetricPDF <- c(IIRMetricPDF,curAUC[sample(length(curAUC),(1.0-infraction)*length(curAUC))]);
								}
							}
						}
					}
					else
					{
						BSWiMS.model <- backVarElimination_Bin(object=update.model$final.model,pvalue=forward.model$theZthr,Outcome=Outcome,data=stdata,startOffset=startOffset,type=type,selectionType=testType);
					}
				}
			}
			else
			{
				forward.model <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=acovariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,featureSize=featureSize);
				update.model <- forward.model$update.model;
				BSWiMS.model <- bootstrapVarElimination_Res(object=update.model$final.model,pvalue=forward.model$p.thresholds,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots);
				if (print)
				{
					cat("Final Formula:",BSWiMS.model$back.formula,"\n");
				}

				if (length(forward.model$var.names)>0)
				{
					if (elimination.bootstrap.steps>1)
					{
						metric <- BSWiMS.model$bootCV$testRMSE;
						if (!is.null(firstModel))
						{
							if (length(attr(terms(formula(BSWiMS.model$back.formula)),"term.labels"))>0)
							{
								sdOutcome <- BSWiMS.model$bootCV$outcomeSD;
								curRMS <- BSWiMS.model$bootCV$testSampledRMSE;
								firstMedRMS <- median(IIRMetricPDF);
								curMedRMS <- median(curRMS);
								firstCount <- sum(curMedRMS <= IIRMetricPDF);
								curCount <- sum(firstMedRMS >= curRMS);
								supchance <- sum(sdOutcome >= curRMS)/length(curRMS);
								infraction <- 1.0 - 0.5*firstCount/length(IIRMetricPDF) - 0.5*curCount/length(curRMS);
								simTest <- ks.test(IIRMetricPDF,curRMS + rnorm(length(curRMS),0,1e-10))$p.value
								if (print) 
								{
									cat("Sd:", sdOutcome,"(",supchance,")",BSWiMS.model$back.formula,": Base: ",firstMedRMS,"(",max(IIRMetricPDF),") Current: ",BSWiMS.model$bootCV$testRMSE,"(",min(BSWiMS.model$bootCV$testSampledRMSE),") Inferior Count:",firstCount," Tests:",length(firstModel$bootCV$testSampledRMSE)," Fraction:",infraction,"\n");
								}
								if (supchance < 0.75) infraction <- 1.0;
								isInferior <- (infraction > 0.975);
								if ( !isInferior && (cycles < 3) && (simTest > 0.05) )
								{
									IIRMetricPDF <- c(IIRMetricPDF,curRMS[sample(length(curRMS),(1.0-infraction)*length(curRMS))]);
								}
							}
						}
					}
					else
					{
						BSWiMS.model <- backVarElimination_Res(object=update.model$final.model,pvalue=forward.model$p.thresholds,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect);
					}
				}
			}
			if (is.null(firstModel))
			{
				if (print)
				{
					cat("First Formula:",BSWiMS.model$back.formula,"\n");
				}
				firstModel <- BSWiMS.model;
				fforward.model <- forward.model;
				fupdate.model <- update.model;
				sdOutcome <- BSWiMS.model$bootCV$outcomeSD;
				isInferior <- (length(attr(terms(formula(firstModel$back.formula)),"term.labels"))==0);
				selectedVariableList <- rownames(variableList[as.numeric(rownames(forward.model$ranked.var)),]);
				if ((testType=="zIDI") || (testType=="zNRI"))
				{
					IIRMetricPDF <- (firstModel$bootCV$sensitivity + firstModel$bootCV$specificity)/2;
					IIRMetricPDF <- IIRMetricPDF + rnorm(length(IIRMetricPDF),0,1e-10);
					if (!isInferior) 
					{
						isInferior <- ( (sum(0.5 <= IIRMetricPDF)/length(IIRMetricPDF)) < 0.5 );
					}
				}
				else
				{
					IIRMetricPDF <- firstModel$bootCV$testSampledRMSE;
					IIRMetricPDF <- IIRMetricPDF + rnorm(length(IIRMetricPDF),0,1e-10);
					if (print)
					{
						cat(length(IIRMetricPDF),": Sup std Outcome",sdOutcome,":",(sum(sdOutcome >= IIRMetricPDF)/length(IIRMetricPDF)),"\n");
					}
					if (!isInferior) 
					{
						isInferior <- ( (sum(sdOutcome >= IIRMetricPDF)/length(IIRMetricPDF)) < 0.5);
					}
				}
				if (isInferior)
				{
					infraction <- 1.0;
				}
			}
			if (!isInferior) #removing the models variables
			{
				cat("+");
				if (print) cat(cycles,":",size,":",nrow(variableList),":",metric,":",infraction,":",BSWiMS.model$back.formula,"\n");
				
#				formula.list <- append(formula.list,BSWiMS.model$back.formula);
				termslist <- attr(terms(formula(BSWiMS.model$back.formula)),"term.labels");
				if (equivalent)
				{
#					print(variableList[as.integer(names(forward.model$ranked.var)),]);
					equmod <- reportEquivalentVariables(BSWiMS.model$back.model,pvalue = 0.25*pvalue,
							  data=data,
							  variableList = variableList[as.integer(names(forward.model$ranked.var)),],
							  Outcome = Outcome,
							  timeOutcome = timeOutcome,								  
							  type = type, osize=nfeat,
							  method="BH");
					formula.list <- append(formula.list,equmod$formula.list);
					if (print) print(equmod$formula.list);
					equiMaxFreq <- equiMaxFreq+max(equmod$bagged$frequencyTable);
					addedEquFreq <- addedEquFreq + 1;
					termslist <- attr(terms(formula(equmod$bagged$formula)),"term.labels");
				}
				else
				{
					formula.list <- append(formula.list,BSWiMS.model$back.formula);
				}
				forward.selection.list <- append(forward.selection.list,forward.model$formula.list);
				if (!is.null(oridinalModels))
				{
#					print(ordinalFormulas);
					oridinalModels$formulas <- append(oridinalModels$formulas,ordinalFormulas);
					for (fn in 1:length(ordinalFormulas))
					{
						termslist <- append(termslist,attr(terms(formula(ordinalFormulas[fn])),"term.labels"));
					}
					termslist <- unique(termslist);
				}
				
				selectedVariableList <- unique(c(selectedVariableList,rownames(variableList[as.numeric(rownames(forward.model$ranked.var)),])));

				if (cycles == max(as.integer(maxCycles/3),5))
				{
					included <- rownames(variableList) %in% selectedVariableList;
					variableList <- variableList[included,]
					size <- nrow(variableList);
				}
				else
				{
					size <- size - length(termslist);
				}

				included <- rownames(variableList) %in% termslist;
#				print(included);
				variableList <- variableList[!included,]
				size <- min(c(size,nrow(variableList)));
			}
			else
			{
				cat("-");
				if (print) cat(cycles,":",size,"\n");
				termslist <- attr(terms(formula(BSWiMS.model$back.formula)),"term.labels");
				if (length(termslist) > 0)
				{
					included <- rownames(variableList) %in% termslist;
					variableList <- variableList[!included,]
					size <- min(c(size,nrow(variableList)));
				}
				if (cycles == 0) 
				{
					firstModel <- NULL;					
					infraction <- 0.5;
				}
				cycles <- cycles + 1;
			}
			cycles <- cycles + 1;
#			if (equivalent) cycles <- maxCycles;
		}
		if (length(formula.list) == 0)
		{
			formula.list <- append(formula.list,baseformula);
			forward.selection.list <- append(forward.selection.list,forward.model$formula.list);
		}
		if (NumberofRepeats>1) formula.list <- append(formula.list,"=-=End=-=");
	}
	cat("]");
	if(is.null(unirank))
	{
		unirank <- invariableList;
		unirank$ZUni <- (nrow(invariableList):1)
	}
	data[,Outcome] <- theOutcome;
	bagg <- NULL;
	if (!is.null(oridinalModels)) 
	{
		bagg <- baggedModel(oridinalModels$formulas,data,"LM",univariate=unirank,useFreq=FALSE); 
		totmodels <- length(theScores)-1;
		repet <- length(oridinalModels$formulas)/totmodels - 1;
		theBaggs <- list();
		theBaggs2 <- list();
		theClassBaggs <- list();
		mc <- 1; 
		for (s in theScores)
		{
			sdata <- data;
			sdata[,Outcome] <- 1*(sdata[,Outcome] == s);
			theClassBaggs[[mc]] <- baggedModel(oridinalModels$formulas,sdata,type="LOGIT",univariate=unirank,useFreq=TRUE);
			
			if (mc <= totmodels)
			{
				olist <- (0:repet)*totmodels+mc;
				sdata <- data
				sdata[,Outcome] <- 1*(sdata[,Outcome] > s);
				theBaggs[[mc]] <- baggedModel(oridinalModels$formulas[olist],sdata,type="LOGIT",univariate=unirank,useFreq=(NumberofRepeats>1));
				sdata <- rbind(subset(data,get(Outcome) == s),subset(data,get(Outcome) == theScores[mc + 1]));
				sdata[,Outcome] <- 1*(sdata[,Outcome] > s);
				theBaggs2[[mc]] <- baggedModel(oridinalModels$formulas[olist],sdata,type="LOGIT",univariate=unirank,useFreq=(NumberofRepeats>1));
			}
			mc <- mc + 1;
		}
		oridinalModels$theBaggedModels <- theBaggs; 
		oridinalModels$redBaggedModels <- theBaggs2; 
		oridinalModels$theClassBaggs <- theClassBaggs;
#		print(oridinalModels$formulas);

		if (!requireNamespace("MASS", quietly = TRUE)) {
			install.packages("MASS", dependencies = TRUE)
		} 
		data[,Outcome] <- as.factor(theOutcome);
		frma <- as.character(bagg$formula);
#		print(frma);
		oridinalModels$formula <- frma;
		oridinalModels$polr <- try(MASS::polr(frma,data));
		if (inherits(oridinalModels$polr, "try-error"))
		{
			cat(frma,"\n");
			warning (paste(frma,": No ordinal model\n"));
		}
	}
	else
	{	
		if (addedEquFreq>0) 
		{
			equiMaxFreq <- equiMaxFreq/addedEquFreq;
		}
		bagg <- baggedModel(formula.list,data,type,univariate=unirank,useFreq=FALSE,equifreqCorrection=equiMaxFreq,n_bootstrap=1);
	}

	

	result <- list(BSWiMS.model=firstModel,
		forward.model=fforward.model,
		update.model=fupdate.model,
		univariate=unirank,
		bagging=bagg,
		formula.list=formula.list,
		forward.selection.list=forward.selection.list,
		oridinalModels=oridinalModels,
		equivalent=equivalent
	);
	class(result) <- c("fitFRESA","BSWiMS");
	
	return (result);
}