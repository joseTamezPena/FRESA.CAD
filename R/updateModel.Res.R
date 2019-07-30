updateModel.Res <-
function(Outcome,covariates="1",pvalue=c(0.025,0.05),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent"), lastTopVariable= 0,timeOutcome="Time",maxTrainModelSize=-1,p.thresholds=NULL)
{
	type <- match.arg(type)

	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}
	covariates <- acovariates;

	varsize = ncol(data)-1;
	loopst = 2;
	
	nvars <- length(VarFrequencyTable);
	
	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(rownames(VarFrequencyTable));
	vnames_model <- vector();
	nsize <- nrow(data)
	if (maxTrainModelSize <= 0)
	{
		maxTrainModelSize = as.integer(nsize/2);
	}

	baseForm = Outcome;
#For Cox  models 
	if (type == "COX")
	{
	  baseForm = paste("Surv(",timeOutcome,",",Outcome,")",sep="");
	}
	varlist <- vector();
	model_ziri <- vector();

	topfreq <- as.integer(0.05*VarFrequencyTable[1]+0.5); #check only features with a 5% bootstrap frequency relative to the top
	if (lastTopVariable < 1) 
	{
		lastTopVariable = sum(1*(VarFrequencyTable>topfreq));
	}
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
#	cat("Top Freq: ",VarFrequencyTable[1],"All Selected Features: ",nvars,"To be tested: ",lastTopVariable,"\n");
		
#	frm1 = paste(baseForm,"~",covariates,"+",vnames[topvarID[1]]);	
	frm1 = paste(baseForm,"~",covariates);
	ftmp <- formula(frm1);
#	varlist <- append(varlist,topvarID[1])
	bestmodel <- modelFitting(ftmp,data,type,TRUE)
#	topvarID[1]=0;
	bestResiduals <- residualForFRESA(bestmodel,data,Outcome);
	termsinserted = 0;
	stdoutput <- sd(data[,Outcome]);
	error <- 1.0;
	tol <- 1.0e-8;
	indexlastinserted <- 1
	for (pval in 1:length(pvalue))
	{
		cthr_a <- pvalue[pval];
		ftmp <- formula(frm1);
		i <- 1;
#		cat("Update: ",frm1,"\n");
		while ((i<=lastTopVariable)&&(error>tol))
		{
			if (nvars>termsinserted) cthr_a <- pvalue[pval]/(nvars-termsinserted);
			if (is.null(p.thresholds))
			{
				cthr <- cthr_a;
			}
			else
			{
				lobs <- termsinserted+1;				
				if (lobs>length(p.thresholds)) lobs <- length(p.thresholds);
				cthr <- p.thresholds[lobs];
				if (cthr>cthr_a) cthr <- cthr_a;
			}
			if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (termsinserted < maxTrainModelSize))
			{
				frma <- paste(frm1,"+",vnames[topvarID[i]]);
				ftmp <- formula(frma);
				newmodel <- modelFitting(ftmp,data,type,TRUE);
				if ( !inherits(newmodel, "try-error"))
				{
					cur_residuals <- residualForFRESA(newmodel,data,Outcome);
					error <- mean(abs(cur_residuals))/stdoutput;
					iprob <- .Call("improvedResidualsCpp",bestResiduals,cur_residuals,testType,0);
					piri <- iprob$p.value;
					if (!is.nan(piri) && !is.na(piri) && (piri<cthr))
					{
						bestmodel <- newmodel;
						bestResiduals <- cur_residuals;
						frm1 <- frma;
						vnames_model <- append(vnames_model,vnames[topvarID[i]]);
						model_ziri <- append(model_ziri,-1*(qnorm(piri)));
						termsinserted = termsinserted + 1;
						if (indexlastinserted <= i) 
						{
							indexlastinserted <- i+1;
						}
						VarFrequencyTable[i]=0;
						topvarID[i]=0;
					}	
				}
			}
			i = i+1;
		}
#		cat(cthr," Update: ",frm1,"\n");
	}
#	cat("Update: ",frm1,"\n");
#    print(p.thresholds[1:termsinserted]);

	environment(bestmodel$formula) <- globalenv()
	environment(bestmodel$terms) <- globalenv()
	
 	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=frm1,
	z.NeRI=model_ziri
	);
  
	return (result);
}
