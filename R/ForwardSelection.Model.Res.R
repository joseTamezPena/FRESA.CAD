ForwardSelection.Model.Res <-
function(size=100,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,variableList,data,maxTrainModelSize=20,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),timeOutcome="Time",cores = 6,randsize = 0,featureSize = 0)
{
#	R_CStackLimit = -1;

	if (is.na(size))
	{
		stop("Size: Number of variables to be explored is not defined\n")
	}

#	cat(featureSize," <- Feat1 Size\n");
	if (featureSize==0) featureSize = max(c(ncol(data)-1,featureSize,nrow(variableList)));
#	cat(featureSize," <- Feat2 Size\n");

	type <- match.arg(type)
	testType <- match.arg(testType)
	Outcome<-as.character(Outcome);

	if (type == "COX")
		timeOutcome<-as.character(timeOutcome)
	else
		timeOutcome=".";

	vnames <- as.vector(variableList[,1]);
	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}
	mcnt=0;
	i=1;
	if (length(vnames)<size) size = length(vnames)
	while ((mcnt==0)&&(i<=size))
	{
		mcnt = mcnt+str_count(vnames[i],"\\*");
		i = i + 1;
	}

	if ( mcnt>0 )
	{
#		cat(Outcome," :",nrow(data),":",nrow(variableList)," Model Frames\n")
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
			frm <- paste(frm,"+",vnames[i]);
		}
#		cat(frm,"\n")
		modelFrame <- model.frame(formula(frm),data);
	}
	else
	{
#		modelFrame <- data;
		varz <- unique(c(covariates,Outcome,timeOutcome,vnames[1:size]));
		varz <- varz[varz %in% colnames(data)];
		modelFrame <- data[,varz];	
	}

	colNames=colnames(modelFrame);
	if (randsize >= 0)
	{
		output<-.Call("ForwardResidualModelCpp",size, fraction, pvalue, loops, covariates, Outcome,vnames, maxTrainModelSize, type, timeOutcome, testType,data.matrix(modelFrame),colNames,featureSize,cores);
	}
	else
	{
		if (timeOutcome!=".") modelFrame[,timeOutcome] <- runif(nrow(modelFrame));
		output<-.Call("ForwardResidualModelCpp",size, fraction, pvalue, loops, "1", paste("RANDOM",Outcome,sep=""),vnames, maxTrainModelSize, type, timeOutcome, testType,data.matrix(modelFrame),colNames,featureSize,cores);
	}

	random.fraction <- 1.0
	if (randsize<0)
	{
		covcount <- 1;
		mcount=0;
		pfind <- 0;
		randsize <- 0;
#		print(output$formula.list);
		for (i in 1:loops)
		{
		    plusc = str_count(output$formula.list[i],"\\+");
			if (plusc>=covcount)
			{
				randsize <- randsize + plusc - covcount;
				pfind <- pfind + 1*(plusc>covcount);
				mcount <- mcount+1;
			}
		}
		random.fraction <- pfind/loops;
		randsize = (nrow(variableList)/size)*(randsize/loops);
		cat ("\n Vars:",nrow(variableList),"Size:",size,sprintf(", Fraction= %6.3f,  Average random size = %6.2f, Size:%6.2f",random.fraction,randsize,randsize/pvalue),"\n");
	}
	else
	{
		if (randsize==0) randsize = pvalue*nrow(variableList);
	}

	mynames <- output$mynames + 1;
	formula.list <- output$formula.list
	
#	print(formula.list);

	topvar <- table(mynames);
		
	if (length(topvar)>1)
	{
		topvar <- topvar[order(-topvar)];
		if (loops > 1)
		{
			oF <- orderFeatures(output$formula.list,univariate=variableList);
			oF <- oF$VarFrequencyTable[names(oF$VarFrequencyTable) %in% rownames(variableList)];
#			print(oF$VarFrequencyTable);
			linspace <- as.character(1:nrow(variableList))
			names(linspace) <- rownames(variableList);
			linspace <- linspace[names(oF)];
			topvar <- topvar[linspace];
		}
	}
	

	update.model <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,pvalue),VarFrequencyTable=topvar,variableList=variableList,data=data,type=type,testType=testType,timeOutcome=timeOutcome,p.thresholds=output$p.thresholds)


	ftmp <- update.model$formula;
	bestmodel <- update.model$final.model;
#	print(output$Base.values);

	base.Zvalues <- -1.0*qnorm(as.vector(output$Base.values));
	
#	print(base.Zvalues);
	
	names(base.Zvalues) <- vnames[1:size];

	result <- list(final.model=bestmodel,
	var.names=update.model$var.names,
	formula=ftmp,
	ranked.var=topvar,
	formula.list=formula.list,
	random.formula.size=randsize,
	random.fraction = random.fraction,
	variableList=variableList,
	base.Zvalues=base.Zvalues,
	p.thresholds=output$p.thresholds,
	update.model=update.model
	);
	
	return (result);
}
