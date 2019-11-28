backVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
{

	back.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),setIntersect=1) 
	{
		type <- match.arg(type);
	  
		varsList <- unlist(as.list(attr(terms(object),"variables")))
		termList <- str_replace_all(attr(terms(object),"term.labels"),":","\\*")
		
		modsize <- length(termList);
		
		removeID = 0;

		outCome = paste(varsList[2]," ~ ");
		if (setIntersect==0) 
		{
			outCome = paste(outCome," 0  ");
		}
		else
		{
			outCome = paste(outCome," 1  ");
		}
		frm1 = outCome;
		if (length(termList)>0)
		{
			for ( i in 1:length(termList))
			{
				frm1 <- paste(frm1,paste("+",termList[i]));
			}
		}
#		cat ("Len: ",length(termList)," : ",frm1,"\n")
		ftmp <- formula(frm1);
		bckform <- frm1;
		startSearch = 1 + startOffset;
		if (length(termList)>1)
		{
			FullModel <- modelFitting(ftmp,data,type,TRUE)
			FullResiduals <- residualForFRESA(FullModel,data,Outcome);

			if (length(termList)>startSearch)
			{
				ploc <- 1+length(termList)-startSearch;
				if (ploc>length(pvalue)) ploc <- length(pvalue);
				cpv <- pvalue[ploc];
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
					if (inherits(redModel, "try-error"))
					{
						redModel <- FullModel
					}
					

					redResiduals <- residualForFRESA(redModel,data,Outcome);
					iprob <- improvedResiduals(redResiduals,FullResiduals,testType);
					if (iprob$p.value>cpv)
					{
						cpv = iprob$p.value;
						removeID = i;
					}
				}
			}
			if ((length(termList) == startSearch) && (removeID == startSearch)) 
			{
				removeID = -1;
			}

			if (removeID>0)
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
		}
		else
		{
			FullModel <- object;
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
	changes2<-0
	while ((changes>0) && (loops<100))
	{

		bk <- back.var.NeRISelection(model,pvalue,Outcome=Outcome,data=data,startOffset,type,testType,setIntersect);
		changes = as.integer(bk$Removed);
		if (changes>0)
		{
		  loops = loops + 1;
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
	}

	modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type=type);
	
	result <- list(back.model= model,
	loops=loops,
	reclas.info=modelReclas,
	back.formula=bk$backfrm,
	bootCV=NULL,
	lastRemoved=changes2,
	at.opt.model=beforeFSCmodel,
#	string.formula=bk$backfrm,
	beforeFSC.formula=formula(beforeFSC.formula));
	return (result);
}