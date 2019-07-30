bootstrapVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops=64,setIntersect=1,print=TRUE,plots=TRUE) 
{
  	testType <- match.arg(testType)
	pvalue <- as.vector(pvalue);
	boot.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops,setIntersect=1,best.formula=NULL) 
	{
		testType <- match.arg(testType)
		type <- match.arg(type);
		FullModel <- object;
		varsList <- unlist(as.list(attr(terms(object),"variables")))
		termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
		
		removeID = 0;
		outCome = paste(varsList[2],"~",setIntersect);
		frm1 = outCome;
		sizemodel <- length(termList);
		for ( i in 1:sizemodel)
		{
			frm1 <- paste(frm1,"+",termList[i]);
		}
		ftmp <- formula(frm1);
		backfrm <- frm1;
		NeRICVp <- bootstrapValidation_Res(1.0,loops,ftmp,Outcome,data,type,plots=plots,bestmodel.formula=best.formula)
		startSearch = 1 + startOffset;
		who = -1;
		if (sizemodel >= startSearch)
		{		
			ploc <- 1+sizemodel-startSearch;
			if (ploc>length(pvalue)) ploc <- length(pvalue);
			maxPvalue <- pvalue[ploc];
#			cat(maxPvalue,":");
			modelReclas <-getVar.Res(object,data,Outcome,type);
			idlist=startOffset+1;
			frm1 = outCome;
			if (startSearch > 1)
			{
				for ( i in 1:(startSearch-1))
				{
					frm1 <- paste(frm1,"+",termList[i]);
				}
			}
			for ( i in startSearch:sizemodel)
			{	
				{
					switch(testType, 
						tStudent = 
						{
							c0 <- modelReclas$tP.value[idlist];
							ci <- median(NeRICVp$tStudent.pvalues[,idlist], na.rm = TRUE);
							ci2 <- median(NeRICVp$test.tStudent.pvalues[,idlist], na.rm = TRUE);
						},
						Wilcox = 
						{ 
							c0 <- modelReclas$WilcoxP.value[idlist];
							ci <- median(NeRICVp$wilcox.pvalues[,idlist], na.rm = TRUE);
							ci2 <- median(NeRICVp$test.wilcox.pvalues[,idlist], na.rm = TRUE);
						},
						Binomial =
						{ 
							c0 <- modelReclas$BinP.value[idlist];
							ci <- median(NeRICVp$bin.pvalues[,idlist], na.rm = TRUE);
							ci2 <- median(NeRICVp$test.bin.pvalues[,idlist], na.rm = TRUE);
						},
						Ftest =
						{ 
							c0 <- modelReclas$FP.value[idlist];
							ci <- median(NeRICVp$F.pvalues[,idlist], na.rm = TRUE);
							ci2 <- median(NeRICVp$test.F.pvalues[,idlist], na.rm = TRUE);
						},
					)
					if (is.nan(ci) || is.na(ci) ) ci <- c0;
					if (is.nan(ci2) || is.na(ci2) ) ci2 <- ci;
					maxp <- max(c(ci2,ci,c0)); 

					if  (maxp >= maxPvalue) # any test p value superior to threshold remove
					{
						maxPvalue = maxp;
						who = i;
					}
				}
				idlist=idlist+1;
			}
			for ( i in startSearch:sizemodel)
			{
				if (who != i)
				{
					frm1 <- paste(frm1,"+",termList[i]);
				}
				else
				{
					removeID=i;
				}
			}
			ftmp <- formula(frm1);
			FullModel <- modelFitting(ftmp,data,type,TRUE)
			backfrm <- frm1
			if (inherits(FullModel, "try-error"))
			{
				cat("Error: Reduced Formula: ",frm1,". Attempting to remove the last term\n");
				who= sizemodel;
				frm1 = outCome;
				if (startSearch > 1)
				{
					for ( i in 1:(startSearch-1))
					{
						frm1 <- paste(frm1,"+",termList[i]);
					}
				}
				for ( i in startSearch:sizemodel)
				{
					if (who != i)
					{
						frm1 <- paste(frm1,"+",termList[i]);
					}
				}
				ftmp <- formula(frm1);
				backfrm <- frm1
				FullModel <- modelFitting(ftmp,data,type,TRUE);
				if (inherits(FullModel, "try-error"))
				{
					cat("Error: Reduced Formula: ",frm1,"Unsuccessful attempt\n");
				}
			}
		}
		if ((sizemodel==1)&&(who>0))
		{
			who = 0;
			removeID = -removeID;
		}
		if (who>0) NeRICVp <- bootstrapValidation_Res(1.0,loops,ftmp,Outcome,data,type,plots=plots)
		testRMSE <- NeRICVp$testRMSE;
		if (length(testRMSE)==0) testRMSE=1.0e10;
		if (is.na(testRMSE)) testRMSE=1.0e10;
		if (is.null(testRMSE)) testRMSE=1.0e10;
		result <- list(Removed=removeID,backfrm=backfrm,testRMSE=testRMSE,bootVal=NeRICVp);

		return (result)
	}

	model <- NULL;
	modelReclas <- NULL;
	best.formula <- NULL;
	changes=1;
	NeRICV <- bootstrapValidation_Res(1.0,loops,formula(object),Outcome,data,type,plots=plots)

	bestbootRMSE <- as.vector(quantile(NeRICV$testSampledRMSE, probs = c(0.05, 0.5, 0.95), na.rm = TRUE,names = FALSE, type = 7));		
	startRMSE <- bestbootRMSE[2];
	minbootRMSE <- bestbootRMSE[2];

	varsList <- unlist(as.list(attr(terms(object),"variables")))
	termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")

	outCome = paste(varsList[2]," ~ ",setIntersect);
	frm1 = outCome;
	if (length(termList)>0)
	{
		for ( i in 1:length(termList))
		{
			frm1 <- paste(frm1,paste("+",termList[i]));
		}
	}
	model.formula <- frm1;
	if (is.null(best.formula)) 
	{
		best.formula <- frm1;
	}
	
	min.formula <- best.formula;
	loopsAux=0;
    model = object;
	beforeFSCmodel <- object;
	beforeFSC.model.formula <- frm1;
	changes = 1*(1<=length(termList));
	changes2 <- 0
	bk <- NULL;
#	print(pvalue[1:10]);
	while ((changes>0) && (loopsAux<100)) 
	{
		bk = boot.var.NeRISelection(object=model,pvalue=pvalue,Outcome=Outcome,
		startOffset=startOffset,type=type,testType=testType,loops=loops,setIntersect=setIntersect,best.formula=best.formula);


		model.formula <- bk$backfrm;
		nmodel <- modelFitting(model.formula,data,type,TRUE);
		changes = bk$Removed;
		if (changes>0)
		{
			testRMSE <- as.vector(quantile(bk$bootVal$testSampledRMSE, probs = c(0.05, 0.5, 0.95), na.rm = TRUE,names = FALSE, type = 7));

			if (testRMSE[2] <= bestbootRMSE[3])
			{
				best.formula <- bk$backfrm;
				if (testRMSE[2] <= minbootRMSE)			
				{
					min.formula <- bk$backfrm;
					minbootRMSE <- testRMSE[2];
				}
				bestbootRMSE <- testRMSE;
			}
			changes2 <- attr(terms(model),"term.labels")[which(!(attr(terms(model),"term.labels") %in% attr(terms(nmodel),"term.labels")))]
			if (length(changes2)>1)
			{
				changes2<-changes2[2]
			}
		}
		model <- nmodel;
		loopsAux = loopsAux + 1;
	}
	if (length(all.vars(formula(model.formula))) > 1)
	{
		modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type);
		NeRICV <- bootstrapValidation_Res(1.0,2*loops,model.formula,Outcome,data,type,plots=plots);
		# if (is.null(bk))
		# {
			# NeRICV <- bootstrapValidation_Res(1.0,loops,model.formula,Outcome,data,type,plots=plots);
		# }
		# else
		# {
			# NeRICV <- bk$bootVal;
			# if (is.null(NeRICV))
			# {
				# NeRICV <- bootstrapValidation_Res(1.0,loops,model.formula,Outcome,data,type,plots=plots);
			# }
		# }
	}
	else
	{
		model.formula <- outCome;
		NeRICV <- bootstrapValidation_Res(1.0,loops,model.formula,Outcome,data,type,plots=plots);
	}
	if (print==TRUE)
	{
		cat("Before BSC   Mod:",beforeFSC.model.formula,"\n");
		cat("Min RMSE Formula:",min.formula,"\n");
		if (!is.null(bk)) 
		{
			cat("Final Formula:",bk$backfrm,"\n")			
		}
		cat("Start RMSE:",startRMSE,"Min RMSE:",bestbootRMSE[2],"final RMSE:",NeRICV$testRMSE,"\n")
	}

	back.model=modelFitting(formula(model.formula),data,type,TRUE);
	environment(back.model$formula) <- globalenv()
	environment(back.model$terms) <- globalenv()
	at.opt.model = modelFitting(formula(min.formula),data,type,TRUE);
	environment(at.opt.model$formula) <- NULL
	environment(at.opt.model$terms) <- NULL
	
	result <- list(back.model=back.model,
	loops=loopsAux,
	reclas.info=modelReclas,
	bootCV=NeRICV,
	back.formula=model.formula,
	lastRemoved=changes2,
	at.opt.model =at.opt.model,
	beforeFSC.formula = beforeFSC.model.formula,
	at.RMSE.formula = best.formula);
	
	return (result);
}