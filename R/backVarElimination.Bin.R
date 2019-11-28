backVarElimination_Bin <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
{
  	seltype <- match.arg(selectionType)

	back.var.IDISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),selectionType=c("zIDI","zNRI")) 
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
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		frm1 = outCome;
		if (length(termList)>0)
		{
			for ( i in 1:length(termList))
			{
				frm1 <- paste(frm1,paste("+",termList[i]));
			}
		}
		else
		{
			frm1 <- paste(frm1," 1 ")
		}
		
		ftmp <- formula(frm1);
		bckform <- frm1;
		FullModel <- modelFitting(ftmp,data,type,TRUE)
		startSearch = 1 + startOffset;
		if ( !inherits(FullModel, "try-error"))
		{
			FullPredict <- predict.fitFRESA(FullModel,data,'prob');
			if (length(termList)>startSearch)
			{
				ploc <- 1+length(termList)-startSearch;
				if (ploc>length(cthr)) ploc <- length(cthr);
				minlcl = cthr[ploc];
				for ( i in startSearch:length(termList))
				{
				
					frm1 = outCome;
					for ( j in 1:length(termList))
					{
						if (i!=j)
						{
							frm1 <- paste(frm1,paste("+",termList[j]));
						}
					}
					ftmp <- formula(frm1);
					redModel <- modelFitting(ftmp,data,type,TRUE)
					if ( !inherits(redModel, "try-error"))
					{
						redPredict <- predict.fitFRESA(redModel,data,'prob');
						iprob <- .Call("improveProbCpp",redPredict,FullPredict,data[,Outcome]);
						if (seltype=="zIDI") 
						{
							ztst = iprob$z.idi;
						}
						else
						{
							ztst = iprob$z.nri;
						}
						if (is.na(ztst)) ztst=0;
						if (ztst<minlcl)
						{
							minlcl = ztst;
							removeID = i;
						}
					}
				}
			}
		}

		if ((length(termList) == startSearch) && (removeID == startSearch)) 
		{
			removeID = -1;
		}
		
		if (removeID > 0)
		{
			frm1 = outCome;
			for ( i in 1:length(termList))
			{
				if (i != removeID)
				{
					frm1 = paste(frm1,paste("+",termList[i]));
				}
			}
			ftmp <- formula(frm1);
			bckform <- frm1;
			FullModel <- modelFitting(ftmp,data,type,TRUE)
		}
		result <- list(Model=FullModel,Removed=removeID,backfrm=bckform);

		return (result)
	}

	bkobj <- NULL;
	beforeFSC.formula <- NULL;

	changes=1;
	loops=0;
    model <- object;
	beforeFSCmodel <- object;
	mydataFrame <- data;
	myOutcome <- Outcome;
	changes2 <- 0;
	while ((changes>0) && (loops<100))
	{
		bk <- back.var.IDISelection(model,pvalue,Outcome=myOutcome,data=mydataFrame,startOffset,type,seltype);

		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  changes2<- attr(terms(model),"term.labels")[which(!(attr(terms(model),"term.labels") %in% attr(terms(bk$Model),"term.labels")))]
		  model = bk$Model;
		  if (length(changes2)>1)
			{
				changes2<-changes2[2]
			}
		}
		if (changes < 0)
		{
			changes2<- changes
		}
		model = bk$Model;
		  loops = loops + 1;
	}
	modelReclas <- getVar.Bin(model,data=mydataFrame,Outcome=myOutcome,type);
	result <- list(back.model=model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=bk$backfrm,
	lastRemoved=changes2,
	at.opt.model=beforeFSCmodel,
#	string.formula=bk$backfrm,
	beforeFSC.formula=formula(beforeFSC.formula));
	return (result);
}