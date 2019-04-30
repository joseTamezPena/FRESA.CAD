bootstrapVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops=64,print=TRUE,plots=TRUE) 
{
  	seltype <- match.arg(selectionType)
	pvalue <- as.vector(pvalue);

	boot.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI"),loops,best.formula=NULL) 
	{
		seltype <- match.arg(selectionType)
		type <- match.arg(type);
		varsList <- unlist(as.list(attr(terms(object),"variables")))
		termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
		
		if (pvalue[1]<0.5) 
		{
			cthr <- abs(qnorm(pvalue)); 
		}
		else 
		{
			cthr <- pvalue;
		}
		removeID <- 0;

		outCome <- paste(varsList[2]," ~ 1");
		frm1 <- outCome;
		testAUC <- 0.5;
		removedTerm <- NULL;
		who <- 0;
		idiCV <- NULL;
		modsize <- length(termList);
		if (modsize>0)
		{
			for ( i in 1:modsize)
			{
				frm1 <- paste(frm1,"+",termList[i]);
			}
#			print(frm1)
			ftmp <- formula(frm1);

			idiCV <- bootstrapValidation_Bin(1.0,loops,ftmp,Outcome,data,type,plots =plots,best.model.formula=best.formula)
			testAUC <- (idiCV$sensitivity + idiCV$specificity)/2;
			testAUC <- median(testAUC,na.rm = TRUE);
			resuBin <- getVar.Bin(object,data,Outcome,type);
			startSearch <- 1 + startOffset;
			frm1  <- outCome;
			if (startSearch > 1)
			{
				for ( i in 1:(startSearch-1))
				{
					frm1 <- paste(frm1,"+",termList[i]);
				}
			}
			if (startSearch <= modsize)
			{
				ploc <- 1+modsize-startSearch;
				if (ploc>length(cthr)) ploc <- length(cthr);
				minlcl <- cthr[ploc];
				idlist <- startOffset+1;
				for ( i in startSearch:modsize )
				{
					{
						if (seltype=="zIDI")
						{
							c0 <- resuBin$z.IDIs[idlist];
							ci <- median(idiCV$z.IDIs[,idlist], na.rm = TRUE);
							ci2 <- median(idiCV$test.z.IDIs[,idlist], na.rm = TRUE);
						}
						else
						{
							c0 <- resuBin$z.NRIs[idlist];
							ci <- median(idiCV$z.NRIs[,idlist], na.rm = TRUE);
							ci2 <- median(idiCV$test.z.NRIs[,idlist], na.rm = TRUE);
						}
						if (is.nan(ci) || is.na(ci) ) ci <- c0;
						if (is.nan(ci2) || is.na(ci2) ) ci2 <- ci;
						minz <- min(c(c0,ci,ci2));
#						cat(c0,":",ci,":",ci2,":",minlcl,":",minz,":",termList[i],"\n");
						if  (minz < minlcl)
						{
							minlcl = minz;
							who = i;
						}
					}
					idlist=idlist+1;
				}
			}
			for ( i in startSearch:modsize)
			{
				if (who != i)
				{
					if (who != -1)
					{
						frm1 <- paste(frm1,"+",termList[i]);
					}
				}
				else
				{
					removeID=i;
					removedTerm=termList[i];
				}
			}
			if ((modsize == startSearch) && (who == startSearch)) 
			{
				removeID = -removeID;
			}
		}
		ftmp <- formula(frm1);

		if ((who>0) && (modsize>1)) idiCV <- bootstrapValidation_Bin(1.0,loops,ftmp,Outcome,data,type,plots=plots)
		afterTestAUC <- (idiCV$sensitivity + idiCV$specificity)/2;
		afterTestAUC <- median(afterTestAUC,na.rm = TRUE);
		if (is.null(afterTestAUC)) afterTestAUC=0.0;
		if (is.null(testAUC)) testAUC=0.5;
		if (is.na(afterTestAUC)) afterTestAUC=0.0;
		if (is.na(testAUC)) testAUC=0.5;
		
		result <- list(Removed=removeID,BootModelAUC=idiCV$blind.ROCAUC$auc,backfrm=frm1,bootval=idiCV,afterTestAUC=afterTestAUC,beforeTestAUC=testAUC,removedTerm=removedTerm);

		return (result)
	}

	bkobj <- NULL;
	
	bestAccuracy <- c(0.5,0.5,0.5);
	best.formula=NULL;

	startAccuracy = bestAccuracy;
	maxAccuracy <- startAccuracy[2];
	
	changes=1;
	loopsAux=0;
    model <- object;
	modelReclas <- NULL;
	myOutcome <- Outcome;

	varsList <- unlist(as.list(attr(terms(object),"variables")))
	termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
	outCome = paste(varsList[2]," ~ 1");
	frm1 = outCome;
	if (length(termList) > 0)
	{
		for ( i in 1:length(termList))
		{
			frm1 <- paste(frm1,paste("+",termList[i]));
		}
	}
	beforeFSCmodel.formula <- frm1;
	model.formula <- frm1;
	if (is.null(best.formula))
	{
		best.formula <- frm1;
	}
	min.formula <- best.formula;
	
	
	beforeFSCmodel <- object;
	beforeFormula <- frm1;
	bk <- NULL;
	changes2 <- 0;
#	print(pvalue[1:10]);
	while ((changes>0) && (loopsAux<100))
	{
		bk <- boot.var.IDISelection(model,pvalue,Outcome=myOutcome,startOffset,type,seltype,loops,best.formula);
		beforeFormula <- bk$backfrm;
		nmodel = modelFitting(formula(bk$backfrm),data,type,TRUE);
		if (!is.null(bk$bootval))
		{
			testAccuracy <-as.vector(quantile(bk$bootval$accuracy, probs = c(0.05, 0.5, 0.95), na.rm = TRUE,names = FALSE, type = 7));
			if (loopsAux == 0) startAccuracy <- bk$beforeTestAUC;
		}

		if ((bk$Removed>0) && (!inherits(nmodel, "try-error")))
		{
			if (!is.null(bk$bootval))
			{	
				if (!is.na(testAccuracy) && !is.null(testAccuracy))
				{
					if (testAccuracy[2] >= bestAccuracy[1])
					{
						best.formula <- bk$backfrm;
						if (testAccuracy[2] >= maxAccuracy)
						{
							min.formula <- bk$backfrm;
							maxAccuracy <- testAccuracy[2];
						}
						bestAccuracy <- testAccuracy;
					}
				}
			}
			if (changes>0)
			{
				changes2<- attr(terms(model),"term.labels")[which(!(attr(terms(model),"term.labels") %in% attr(terms(nmodel),"term.labels")))]

				if (length(changes2)>1)
				{
					changes2<-changes2[2]
				}
			}
		}
		changes = as.integer(bk$Removed);
		model <- nmodel;
		model.formula <- bk$backfrm;
		
		loopsAux = loopsAux + 1
	}
	idiCV <- NULL;
	if (length(all.vars(formula(model.formula))) > 1)
	{
		modelReclas <- getVar.Bin(model,data=data,Outcome=myOutcome,type);
		idiCV <- bootstrapValidation_Bin(1.0000,2*loops,formula(model.formula),myOutcome,data,type,plots=plots);
		# if (is.null(bk))
		# {
			# idiCV <- bootstrapValidation_Bin(1.0000,loops,formula(model.formula),myOutcome,data,type,plots=plots);
		# }
		# else
		# {
			# idiCV <- bk$bootval;
			# if (is.null(idiCV))
			# {
				# idiCV <- bootstrapValidation_Bin(1.0000,loops,formula(model.formula),myOutcome,data,type,plots=plots);
			# }
		# }
	}
	else
	{
		model.formula <- outCome;
		idiCV <- bootstrapValidation_Bin(1.0000,loops,formula(model.formula),myOutcome,data,type,plots=plots);
	}
	testAccuracy <-as.vector(quantile(idiCV$accuracy, probs = c(0.05, 0.5, 0.95), na.rm = TRUE,names = FALSE, type = 7));
	if (print == TRUE)
	{
		cat("Before FSC Mod:",beforeFSCmodel.formula,"\n")
		cat("At Acc  Model :",min.formula,"\n")
		cat("Reduced Model :",model.formula,"\n")
		cat("Start AUC:",startAccuracy,"last AUC:",idiCV$blind.ROCAUC$auc,"Accuracy:",testAccuracy[2],"\n")
	}

	back.model<-modelFitting(formula(model.formula),data,type,TRUE);
	environment(back.model$formula) <- globalenv()
	environment(back.model$terms) <- globalenv()


	at.opt.model<-modelFitting(formula(min.formula),data,type,TRUE);
	environment(at.opt.model$formula) <- NULL
	environment(at.opt.model$terms) <- NULL

	result <- list(back.model=back.model,
	loops=loopsAux,
	reclas.info=modelReclas,
	bootCV=idiCV,
	back.formula=model.formula,
	lastRemoved=changes2,
	at.opt.model=at.opt.model,
	beforeFSC.formula=beforeFSCmodel.formula,
	at.Accuracy.formula=best.formula);
	return (result);
}