ForwardSelection.Model.Bin <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),timeOutcome="Time",selectionType=c("zIDI","zNRI"),cores=6,randsize = 0,featureSize=0)
{
#	    R_CStackLimit = -1;
	type <- match.arg(type)
	seltype <- match.arg(selectionType)
	if (featureSize==0) featureSize = max(c(ncol(data)-1,featureSize,nrow(variableList)));
	
#	print(tracemem(data))


#	cat(featureSize," <- F Size\n");


	Outcome<-as.character(Outcome);

	if (type == "COX")
		timeOutcome<-as.character(timeOutcome)
	else
		timeOutcome=".";
	if (is.na(size))
	{
	  stop("Number of variables to be used is not defined\n")
	}


#		cat(timeOutcome," <- Time Outcome\n");
	
	vnames <- as.vector(variableList[,1]);
	mcnt=0;
	i=1;
	if (length(vnames)<size) size = length(vnames)
	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}

	while ((mcnt==0)&&(i<=size))
	{
		mcnt = mcnt+str_count(vnames[i],"\\*");
		i = i + 1;
	}
	
	if (mcnt>0)
	{
		if (nrow(variableList)<size) size = nrow(variableList)
		if (timeOutcome == ".") 
		{
			frm <- paste(Outcome,"~",acovariates);
		}
		else
		{
			frm <- paste(Outcome,"~",acovariates,"+",timeOutcome);
		}
		for (i in 1:size)
		{
			frm <- paste(frm,"+",vnames[i])
		}
#			cat(frm,"\n")
		modelFrame <- model.frame(formula(frm),data);
	}
	else
	{
		varz <- unique(c(Outcome,timeOutcome,covariates,vnames[1:size]));
		varz <- varz[varz %in% colnames(data)];
		modelFrame <- data[,varz];	
	}

	
	colNames=colnames(modelFrame);
#	print(colNames);
#	print(variableList[,1]);
#	print(covariates);
	
	if (randsize >= 0)
	{
		output<-.Call("ReclassificationFRESAModelCpp",size, fraction, pvalue, loops, covariates, Outcome,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,data.matrix(modelFrame),colNames,featureSize,cores);
	}
	else
	{
		if (timeOutcome != ".") modelFrame[,timeOutcome] <- runif(nrow(modelFrame));
		output <-.Call("ReclassificationFRESAModelCpp",size, fraction, pvalue, loops, "1", paste("RANDOM",Outcome,sep="") ,as.vector(variableList[,1]), maxTrainModelSize, type, timeOutcome, seltype,data.matrix(modelFrame),colNames,featureSize,cores);
	}

	random.fraction <- 1.0
	if (randsize<0)
	{
		covcount <- 1;
		mcount <- 0;
		pfind <- 0;
		randsize <- 0;
#			print(output$formula.list);
		for (i in 1:loops)
		{
			plusc = str_count(output$formula.list[i],"\\+");
			if (plusc>=covcount)
			{
				randsize <- randsize + plusc - covcount;
				mcount <- mcount+1;
				pfind <- pfind + 1*(plusc>covcount);
			}
		}
		randsize = (nrow(variableList)/size)*(randsize/loops);
		random.fraction <- pfind/loops;
		cat ("\n Vars:",nrow(variableList),"Size:",size,sprintf(", Fraction= %6.4f,  Average random size = %6.2f, Size:%6.2f",random.fraction,randsize,randsize/pvalue),"\n");
	}
	else
	{
		if (randsize==0) randsize = pvalue*nrow(variableList);
	}
	
	base.Zvalues <- output$Base.values;
	names(base.Zvalues) <- vnames[1:nrow(base.Zvalues)];
	
	mynames <- output$mynames + 1 
#		print(mynames);
	topvar <- table(mynames);
#		print(topvar);
	if (length(topvar)>1)
	{
		topvar <- topvar[order(-topvar)];
		if (loops > 1)
		{
#			print(topvar);
			oF <- orderFeatures(output$formula.list,univariate=variableList);
#			print(names(oF$VarFrequencyTable));
			oF <- oF$VarFrequencyTable[names(oF$VarFrequencyTable) %in% rownames(variableList)];
#			print(names(oF));
			linspace <- as.character(1:nrow(variableList));
			names(linspace) <- rownames(variableList);
#			print(names(linspace));
			linspace <- linspace[names(oF)];
#			print(names(linspace));
			topvar <- topvar[linspace];
#			print(topvar);
		}
	}

#		print(output$formula.list);
#		print(topvar);

	
	update.model <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,pvalue),VarFrequencyTable=topvar,variableList=variableList,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,zthrs=output$Zthr)
		
	

	ftmp <- formula(update.model$formula);
	bestmodel <- update.model$final.model;

	result <- list(final.model=bestmodel,
	var.names=update.model$var.names,
	formula=ftmp,
	ranked.var=topvar,
	formula.list=output$formula.list,
	random.formula.size=randsize,
	random.fraction = random.fraction,
	variableList=variableList,
	base.Zvalues=base.Zvalues,
	theZthr=output$Zthr,
	update.model = update.model
	);
#		cat ("Final :",frm1,"\n")
	return (result);
}
	