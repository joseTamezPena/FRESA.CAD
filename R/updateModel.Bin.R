updateModel.Bin <-
function(Outcome,covariates="1",pvalue=c(0.025,0.05),VarFrequencyTable,variableList,data,type=c("LM","LOGIT","COX"), lastTopVariable= 0,timeOutcome="Time",selectionType=c("zIDI","zNRI"),maxTrainModelSize=0,zthrs=NULL)
{
	type <- match.arg(type)
  	seltype <- match.arg(selectionType)

	acovariates <- covariates[1];
	if (length(covariates)>1)
	{
		for (i in 2:length(covariates))
		{	
			acovariates <- paste(acovariates,"+",covariates[i])
		}
	}
	covariates <- acovariates;

	vnames <- as.vector(variableList[,1]);
	topvarID <- as.numeric(names(VarFrequencyTable));
#	print(topvarID);
	
	if (maxTrainModelSize == 0)
	{
		maxTrainModelSize = as.integer(ncol(data)/2);
	}

	nvars <- length(VarFrequencyTable);
	baseForm = Outcome;
	theoutcome <- data[,Outcome];
					
	#For Cox  models 
	if (type == "COX")
	{
		baseForm = paste("Surv(",timeOutcome,",",Outcome,")",sep="");
	}
		
	vnames_model <- vector();
	model_zmin <- vector();

	
	topfreq <- as.integer(0.05*VarFrequencyTable[1]+0.5); #check only features with a 5% bootstrap frequency relative to the top
	if (lastTopVariable < 1) 
	{
		lastTopVariable = sum(1*(VarFrequencyTable > topfreq));
	}
	if (lastTopVariable > length(VarFrequencyTable)) lastTopVariable = length(VarFrequencyTable);
#	cat("Top Freq: ",VarFrequencyTable[1],"All Selected Features: ",nvars,"To be tested: ",lastTopVariable,"\n");
#	print(variableList[topvarID,1]);


	bestmodel <- NULL;
	ftmp <- NULL;
	frm1 <- NULL;
	error <- 1.0;
	tol = 1.0e-8; # if error less than tol exit
	if (lastTopVariable>0)
	{
#		frm1 = paste(baseForm,"~",covariates,"+",vnames[topvarID[1]]);
		frm1 = paste(baseForm,"~",covariates);

#		vnames_model <- append(vnames_model,vnames[topvarID[1]]);
		model_zmin <- append(model_zmin,NA);
#		topvarID[1] = 0;
		termsinserted = 0;
		indexlastinserted = 1;
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE)
		if ( !inherits(bestmodel, "try-error"))
		{
			bestpredict <- predict.fitFRESA(bestmodel,data,'prob');
			for (pidx in 1:length(pvalue))
			{
				cthr_a <- abs(qnorm(pvalue[pidx]));	
				i = 1;
				while ((i<=lastTopVariable)&&(error>tol))
				{
					if (termsinserted<nvars) cthr_a <- abs(qnorm(pvalue[pidx]/(nvars-termsinserted)));
					if (is.null(zthrs))
					{
						cthr <- cthr_a;
					}
					else
					{
						lobs <- termsinserted+1;
						if (lobs>length(zthrs)) lobs <- length(zthrs);
						cthr <- zthrs[lobs];
						if (cthr<cthr_a) cthr <- cthr_a;
					}
					if ((VarFrequencyTable[i]>0) && (topvarID[i]>0) && (termsinserted < maxTrainModelSize))
					{
						frma <- paste(frm1,"+",vnames[topvarID[i]]);
						ftmp <- formula(frma);
						newmodel <- modelFitting(ftmp,data,type,TRUE);
						if ( !inherits(newmodel, "try-error"))
						{
							curpredict <- predict.fitFRESA(newmodel,data,'prob');
							iprob_t <- .Call("improveProbCpp",bestpredict,curpredict,theoutcome);
							if (seltype=="zIDI") 
							{
								zmin <- iprob_t$z.idi;
							}
							else
							{
								zmin <- iprob_t$z.nri;
							}
							if (!is.nan(zmin) && !is.na(zmin))
							{
								error <- mean(abs(curpredict-theoutcome));
								if (zmin>cthr)
								{
									bestpredict <-curpredict;
									bestmodel <- newmodel;
									frm1 <- frma;
									vnames_model <- append(vnames_model,vnames[topvarID[i]]);
									model_zmin <- append(model_zmin,zmin);
									termsinserted = termsinserted + 1;
									if (indexlastinserted <= i) 
									{
										indexlastinserted = i+1;
									}
									topvarID[i] = 0;
								}
							}
						}
					}
					i = i+1;
				}
#				cat (pidx," :",frm1,"\n")
			}
		}
		ftmp <- formula(frm1);
	}
#	print(zthrs[1:termsinserted])

	if (length(vnames_model)==0)
	{
		frm1 = paste(baseForm,"~",covariates,"+",vnames[topvarID[1]]);
		ftmp <- formula(frm1);
		bestmodel <- modelFitting(ftmp,data,type,TRUE)
	}
	else
	{
		bestmodel <- modelFitting(ftmp,data,type,TRUE)
	}
#	print(summary(bestmodel));

	
	environment(bestmodel$formula) <- globalenv()
	environment(bestmodel$terms) <- globalenv()
	
  	result <- list(final.model=bestmodel,
	var.names=vnames_model,
	formula=frm1,
	z.selectionType=model_zmin
	);
  
	return (result);
}
