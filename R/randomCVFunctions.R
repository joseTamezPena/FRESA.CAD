#Filter Methods


FRESAcacheEnv <- new.env();

correlated_Remove <- function(data= NULL,fnames= NULL,thr=0.999,isDataCorMatrix=FALSE)
{
#	cat(thr,":",length(fnames)," ->");
	if (length(fnames)>1)
	{
#		print(fnames);
		if (!isDataCorMatrix)
		{
			colsd <- apply(data[,fnames],2,sd,na.rm = TRUE);
			fnames <- fnames[colsd > 0];
		}
		
		if (length(fnames)>1)
		{
			corm <- NULL;
			if (isDataCorMatrix)
			{
				corm <- abs(data[,fnames]);
				corm <- corm[fnames,];
			}
			else
			{
				corm <- abs(cor(data[,fnames],method="spearman"));
				diag(corm) <- 0;
			}
			keep <- numeric(length(fnames));
			for (i in length(fnames):1)
			{
#				plot(corm[,i],main=fnames[i])
				keep[i] <- (sum(corm[,i] > thr) == 0);
				if (!keep[i])  
				{
					corm[i,] <- 0;
				}
			}
			attributes(fnames) <- list(removed=fnames[keep == 0]);
			fnames <- fnames[keep == 1];
			attr(fnames,"CorrMatrix") <- corm;
		}
	}
#	cat(length(fnames),"\n");

	return (fnames);
}

correlated_RemoveToLimit <- function(data,unitPvalues,limit=0,thr=0.975,maxLoops=50,minCorr=0.50)
{
	if ((limit >= 0) && (thr < 1.0))
	{
		ntest <- 0;
		cthr <- thr;
		if (limit > 0)
		{
			slimit <- limit;
			if (limit <= 1)
			{
				slimit <- as.integer(limit*nrow(data));
			}
			if (slimit < 2) 
			{
				slimit <- 2;
			}
			if (length(unitPvalues) > slimit)
			{
				pvalmin <- min(unitPvalues)
				pvalatlimin <- unitPvalues[order(unitPvalues)][slimit]
				maxpvalue <- max(c(1000*pvalmin,100*pvalatlimin,1.0e-9));
				unitPvalues <- unitPvalues[unitPvalues <= maxpvalue];
				cormat <- correlated_Remove(data,names(unitPvalues),cthr)
				unitPvalues <- unitPvalues[cormat];
				cormat <- attr(cormat,"CorrMatrix");
				while ( (length(unitPvalues) > slimit) && (ntest < maxLoops) && (cthr > minCorr) )
				{
					unitPvalues <- unitPvalues[correlated_Remove(cormat,names(unitPvalues),cthr,isDataCorMatrix=TRUE)];
					ntest = ntest + 1;
					cthr = cthr*thr;
				}
				if (length(unitPvalues) > slimit)
				{
					unitPvalues <- unitPvalues[1:slimit];
				}
			}
		}
		else
		{
			unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues),cthr)];
		}
	}
	return (unitPvalues);

}


univariate_Logit <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", uniTest=c("zIDI","zNRI"),limit=0,...,n = 0)
{
	varlist <- colnames(data);
	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		  baseformula <- as.character(Outcome);
		  baseformula[3] <- str_replace_all(baseformula[3],"[.]","1");
		  baseformula <- paste(baseformula[2],"~",baseformula[3]);
		  Outcome <- formula(baseformula);
		  olist <- attr(terms(Outcome),"variables")
		  dependent <- as.character(olist[[2]])
		  
		  if (length(dependent)==3)
		  {

			  timeOutcome = dependent[2];
			  status = dependent[3];
	#		  print(c(status,timeOutcome));

			  varlist <- varlist[!varlist %in% c(timeOutcome,status)];
			  varlist <- cbind(varlist,varlist);
			  univ <- ForwardSelection.Model.Bin(nrow(varlist),1.0,0.0,1,"1",status,varlist,data,1,type="COX",selectionType=uniTest,timeOutcome=timeOutcome);
			  unitPvalues <- 2.0*(1.0 - pnorm(univ$base.Zvalues));
			  unitPvalues[unitPvalues > 1.0] <- 1.0;
		  }
		  else
		  {
			status = dependent[1];
 		    varlist <- varlist[status != varlist]
		    varlist <- cbind(varlist,varlist);
			univ <- ForwardSelection.Model.Bin(nrow(varlist),1.0,0.0,1,"1",status,varlist,data,1,type="LOGIT",selectionType=uniTest);
			unitPvalues <- 2.0*(1.0 - pnorm(univ$base.Zvalues));
			unitPvalues[unitPvalues > 1.0] <- 1.0;
		  }
	}
	else
	{ ## It is a standard binary outcome
	#    if (inherits(data[,Outcome],"factor")) data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
		varlist <- varlist[Outcome != varlist];
		nvar <- length(varlist);
		varlist <- cbind(varlist,varlist);
		uniTest <- match.arg(uniTest);
		outcomes <- data[,Outcome];
		outlist <- unique(outcomes)
		pvalue <- pvalue/sqrt(length(outlist) - 1.0)
		unitPvalues <- rep(1.0,nvar);
		for (oft in outlist[-1])
		{
			datat <- rbind(subset(data,outcomes == oft),subset(data,outcomes != oft))
			datat[,Outcome] <- 1*(datat[,Outcome] == oft);
			univ <- ForwardSelection.Model.Bin(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,datat,1,type="LOGIT",selectionType=uniTest);
			pval <- 2.0*(1.0 - pnorm(univ$base.Zvalues));
			pval[pval > 1.0] <- 1.0;
			unitPvalues <- pmin(unitPvalues,pval);
		}
	}

	names(unitPvalues) <-  varlist[,1];
	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
#	print(unitPvalues)
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
   return(unitPvalues);
}

univariate_residual <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",uniTest=c("Ftest","Binomial","Wilcox","tStudent"),type=c("LM","LOGIT"),limit=0,...,n = 0)
{

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}


	varlist <- colnames(data);
#    if (inherits(data[,Outcome],"factor")) data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	varlist <- varlist[Outcome != varlist];
	nvar <- length(varlist);
	varlist <- cbind(varlist,varlist)
	uniTest <- match.arg(uniTest);
	type <- match.arg(type);

	outcomes <- data[,Outcome];
	outlist <- unique(outcomes)
	pvalue <- pvalue/sqrt(length(outlist) - 1.0)
	unitPvalues <- rep(1.0,nvar);
	for (oft in outlist[-1])
	{
		datat <- rbind(subset(data,outcomes == oft),subset(data,outcomes != oft))
		datat[,Outcome] <- 1*(datat[,Outcome] == oft);
		univ <- ForwardSelection.Model.Res(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,datat,1,type=type,testType=uniTest);
		pval <- 2.0*(1.0 - pnorm(univ$base.Zvalues));
		pval[pval > 1.0] <- 1.0;
		unitPvalues <- pmin(unitPvalues,pval);
	}



#	print(unitPvalues);
	names(unitPvalues) <- varlist[,1];
	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
	return(unitPvalues);
}

univariate_DTS <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{
if (!requireNamespace("twosamples", quietly = TRUE)) {
	 install.packages("twosamples", dependencies = TRUE)
} 

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}

	varlist <-colnames(data);
    if (inherits(data[,Outcome], "factor")) data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;
	
	nboots = max(2000,10.0*length(varlist));

	for (j in varlist) 
	{
		tb <- table(data[,j],data[,Outcome])
		if (nrow(tb) > 5)
		{
			 pval <- twosamples::dts_test(control[,j],case[,j],nboots = nboots)["P-Value"]
		}
		else
		{
			pval <- prop.test(tb)$p.value
		}
		 if (inherits(pval, "try-error")) {pval <- 1.0;}
		 if (is.null(pval)) { pval <- 1.0; }
		 if (is.na(pval)) {pval <- 1.0;}
		 unitPvalues[j] <- pval;
	}
  	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
   return(unitPvalues);
}


univariate_filter <- function(data=NULL, Outcome=NULL, pvalue=0.2,pvalueMethod=wilcox.test, adjustMethod="BH",limit=0,...,n = 0)
{

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}


	varlist <-colnames(data);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- rep(1.0,length(varlist));
	names(unitPvalues) <- varlist;

	outcomes <- as.character(data[,Outcome]);
	outlist <- unique(outcomes)
	pvalue <- pvalue/sqrt(length(outlist) - 1.0)
	for (oft in outlist[-1])
	{
		case <- subset(data,outcomes == oft);
		control <- subset(data,outcomes != oft);
		for (j in varlist) 
		{
			tb <- table(data[,j],outcomes == oft)
			if (nrow(tb) > 5)
			{
				 pval <- pvalueMethod(control[,j],case[,j],na.action = na.exclude)$p.value; 
			}
			else
			{
				pval <- prop.test(tb)$p.value
			}
			if (inherits(pval, "try-error")) {pval <- 1.0;}
			if (is.null(pval)) { pval <- 1.0; }
			if (is.na(pval)) {pval <- 1.0;}
			unitPvalues[j] <- min(unitPvalues[j],pval);
		}
	}

  	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
   return(unitPvalues);
}


univariate_KS <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{
	unitPvalues <- univariate_filter(data=data,
										Outcome=Outcome,
										pvalue=pvalue,
										pvalueMethod=ks.test,
										adjustMethod=adjustMethod,
										limit=limit,
										...,
										n = n)


   return(unitPvalues);
}


univariate_Wilcoxon <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	unitPvalues <- univariate_filter(data=data,
										Outcome=Outcome,
										pvalue=pvalue,
										pvalueMethod=wilcox.test,
										adjustMethod=adjustMethod,
										limit=limit,
										...,
										n = n)

   return(unitPvalues);
}

univariate_tstudent <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	unitPvalues <- univariate_filter(data=data,
										Outcome=Outcome,
										pvalue=pvalue,
										pvalueMethod=t.test,
										adjustMethod=adjustMethod,
										limit=limit,
										...,
										n = n)

    return(unitPvalues);
}

univariate_correlation <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", method = "kendall",limit=0,...,n = 0)
{

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}
	varlist <- colnames(data);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- rep(1.0,length(varlist));
	names(unitPvalues) <- varlist;
	
    if (inherits(data[,Outcome], "factor")) 
	{
		outcomes <- as.character(data[,Outcome]);
		outlist <- unique(outcomes);
		pvalue <- pvalue/sqrt(length(outlist) - 1.0);
		for (ots in outlist[-1])
		{
			for (j in varlist) 
			{
				pval <- try(cor.test(1*(outcomes==ots),data[,j],na.action = na.exclude,method = method)$p.value);
				if (inherits(pval, "try-error")) {pval <- 1.0;}
				if (is.null(pval)) { pval <- 1.0; }
				if (is.na(pval)) {pval <- 1.0;}
				unitPvalues[j] <- min(unitPvalues[j],pval);
			}
		}
	}
	else
	{
		for (j in varlist) 
		{
			pval <- try(cor.test(data[,Outcome],data[,j],na.action = na.exclude,method = method)$p.value);
			if (inherits(pval, "try-error")) {pval <- 1.0;}
			if (is.null(pval)) { pval <- 1.0; }
			if (is.na(pval)) {pval <- 1.0;}
			unitPvalues[j] <-  pval;
		}
	}
  	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
   return(unitPvalues);
}

mRMR.classic_FRESA <- function(data=NULL, Outcome=NULL,feature_count=0,...)
{
	if (!requireNamespace("mRMRe", quietly = TRUE)) {
	   install.packages("mRMRe", dependencies = TRUE)
	} 

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}


    if (inherits(data[,Outcome], "factor")) data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	if (feature_count == 0)
	{
		feature_count <- min(c(nrow(data)-1,ncol(data)-1));
	}
	fs <- try(mRMRe::mRMR.classic(data = mRMRe::mRMR.data(data), target_indices = grep(Outcome, colnames(data)), feature_count = feature_count,...));
	if ( !inherits(fs, "try-error"))
	{
		result <- as.vector(mRMRe::scores(fs)[[1]]);
		names(result) <- colnames(data)[mRMRe::solutions(fs)[[1]]];
		result <- result[!is.na(result)];
		gz <- (result >= 0);
		if (sum(1*gz) > 5)	result <- result[gz];
	}
	else
	{
		warning("mRMR.classic Error. I'll make all columns numeric and try again\n")
		data[,1:ncol(data)] <- sapply(data,as.numeric)
		fs <- try(mRMRe::mRMR.classic(data = mRMRe::mRMR.data(data), target_indices = grep(Outcome, colnames(data)), feature_count = feature_count,...));
		if ( !inherits(fs, "try-error"))
		{
			result <- as.vector(mRMRe::scores(fs)[[1]]);
			names(result) <- colnames(data)[mRMRe::solutions(fs)[[1]]];
			result <- result[!is.na(result)];
			gz <- (result >= 0);
			if (sum(1*gz) > 5)	result <- result[gz];
		}
		else
		{
			warning("mRMR.classic Error. Returning the first elements of the data-frame\n")
			result <- 1:feature_count;
			varlist <- colnames(data);
			varlist <- varlist[Outcome != varlist];
			names(result) <- varlist[result];
			result <- result/feature_count;
		}
	}
	return(result);
}

univariate_BinEnsemble <- function(data,Outcome,pvalue=0.2,limit=0,adjustMethod="BH",...)
{
  allf <- numeric();
  varcount <- numeric(ncol(data));
  rankVar <- numeric(ncol(data));
  names(varcount) <- colnames(data);
  names(rankVar) <- colnames(data);

#  if (inherits(data[,Outcome], "factor")) data[,Outcome] <- as.numeric(as.character(data[,Outcome]));

	if (inherits(Outcome,"formula")) ## We assume that it is Survival Object
	{
		dependent <- all.vars(Outcome)
		Outcome = dependent[1];
		if (length(dependent) == 3)
		{
			Outcome = dependent[2];
		}
	}

  
  data <- data[,c(Outcome,correlated_Remove(data,colnames(data)[!(colnames(data) %in% Outcome)]))]

  pvallist <- list();
  pvaltest <- univariate_Logit(data,Outcome,pvalue=pvalue,limit=-1,uniTest="zNRI",adjustMethod=adjustMethod);
#  cat("zNRI")
  geomMeanpVal <- attr(pvaltest,"Unadjusted");
  maxPval <- geomMeanpVal
  minPval <- geomMeanpVal

  pvallist$LogitNRI <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)] + 1;
  rankVar[names(pvaltest)] <- c(1:length(pvaltest))

  wilcxf <- univariate_Wilcoxon(data,Outcome,pvalue=pvalue,limit=-1,adjustMethod=adjustMethod);
#  cat("->Wilcoxon")
  geomMeanpVal <- geomMeanpVal*(attr(wilcxf,"Unadjusted")[names(geomMeanpVal)]);
  maxPval <- pmax(maxPval,attr(wilcxf,"Unadjusted")[names(geomMeanpVal)]);
  minPval <- pmin(minPval,attr(wilcxf,"Unadjusted")[names(geomMeanpVal)]);
  pvallist$Wilcox <- wilcxf;
  varcount[names(wilcxf)] <- varcount[names(wilcxf)]+1;
  rankVar[names(wilcxf)] <- rankVar[names(wilcxf)] + c(1:length(wilcxf));

  features <- intersect(names(pvaltest),names(wilcxf));
  both <- pmin(pvaltest[features],wilcxf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],wilcxf[!(names(wilcxf) %in% features)],both);


	pvaltest <- univariate_residual(data,Outcome,pvalue=pvalue,limit=-1,uniTest="tStudent",type="LOGIT",adjustMethod=adjustMethod);
	#cat("->tstudent")
	 geomMeanpVal <- geomMeanpVal*(attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
	 maxPval <- pmax(maxPval,attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
	 minPval <- pmin(minPval,attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
	 pvallist$tstudent <- pvaltest;
	 varcount[names(pvaltest)] <- varcount[names(pvaltest)] + 1;
	 rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
	 features <- intersect(names(pvaltest),names(allf));
	 both <- pmin(pvaltest[features],allf[features]);
	 allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);
  
  pvaltest <- univariate_KS(data,Outcome,pvalue=pvalue,limit=-1,adjustMethod=adjustMethod)
#  cat("->KS")
  geomMeanpVal <- geomMeanpVal*(attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
  maxPval <- pmax(maxPval,attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
  minPval <- pmin(minPval,attr(pvaltest,"Unadjusted")[names(geomMeanpVal)]);
  pvallist$KS <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)] + 1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

# Removing the largest p-value for the estimation of the geometric mean p-value
  expgeom <- 1.0/(length(pvallist)-1);
  maxPval[maxPval < 1.0e-16] <- 1.0e-16;
  geomMeanpVal <- (geomMeanpVal/maxPval)^expgeom;

# Sorting and top selection
  geomMeanpVal <- geomMeanpVal[order(geomMeanpVal-varcount[names(geomMeanpVal)])];
  allf <- p.adjust(geomMeanpVal,adjustMethod);
  adjusptedp <- allf;
  top <- allf[1];
  if ((top < 0.45) && (geomMeanpVal[1] < 0.05))
  {
	top <- allf[allf <= 1.01*allf[1]];
  }

  allf <- allf[allf <= pvalue]; # Removing after adjusting
  
  varcount <- varcount[names(geomMeanpVal)];
  rankVar <- rankVar[names(geomMeanpVal)];

#  print(names(allf))
  allf <- correlated_RemoveToLimit(data,allf,limit=limit,...);
  if (length(allf) < 2) 
  {
	allf <- c(allf,top[!(names(top) %in% names(allf))]);
  }
#  print(names(allf))

  attr(allf,"varcount") <- varcount;
  attr(allf,"geomMeanpVal") <- geomMeanpVal;
  attr(allf,"adjusptedpGeom") <- adjusptedp;
  attr(allf,"Pvalues") <- pvallist;
  attr(allf,"rankVar") <- rankVar;
  return (allf);
}


univariate_Strata <- function(data,Outcome,pvalue=0.2,limit=0,adjustMethod="BH",unifilter=univariate_BinEnsemble,strata="Gender",...)
{
	statatable <- table(data[,strata]);
	pvalues <- unifilter(data=data,Outcome=Outcome,pvalue=pvalue,limit=limit,adjustMethod=adjustMethod,...);
	for (st in as.integer(names(statatable)))
	{
		spvalues <- unifilter(data=data[data[,strata]==st,],Outcome=Outcome,pvalue=pvalue,limit=limit,adjustMethod=adjustMethod,...);
		features <- intersect(names(spvalues),names(pvalues));
		both <- pmin(pvalues[features],spvalues[features]);
		pvalues <- c(both,pvalues[!(names(pvalues) %in% features)],spvalues[!(names(spvalues) %in% features)]);
	}
	pvalues <- pvalues[order(pvalues)]
	top <- pvalues[1];
	if (top < 0.45)
	{
		top <- pvalues[pvalues <= 1.01*pvalues[1]];
	}
	pvalues <- pvalues[pvalues <= pvalue];
	if (length(pvalues) < 2)
	{
		pvalues <- top;
	}
	else
	{
		pvalues <- correlated_RemoveToLimit(data,pvalues,limit=limit,...);
	}
	return (pvalues)

}


sperman95ci <- function(datatest,nss=4000)
{
	sz <- nrow(datatest)
	sesci <- c(0,0,0);
	if (sz>2)
	{
		ses <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace=TRUE),];
		  ses[i] <- cor(bootsample[,1],bootsample[,2],method="spearman");
		}
		sesci <- quantile(ses, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (sesci);
}

MAE95ci <- function(datatest,nss=4000)
{
	sz <- nrow(datatest)
	sesci <- c(0,0,0);
	if (sz>1)
	{
		ses <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace=TRUE),];
		  ses[i] <- mean(abs(bootsample[,1]-bootsample[,2]));
		}
		sesci <- quantile(ses, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (sesci);
}
	
ClassMetric95ci <- function(datatest,nss=4000)
{
	  sz <- nrow(datatest)
	  accci <- c(0,0,0);
	  senci <- c(0,0,0);
	  aucci <- c(0,0,0);
	  berci <- c(0,0,0);
	  preci <- c(0,0,0);
	  F1ci <- c(0,0,0);
	  tbstat <- as.matrix(table(datatest[,1]));
	  nscores <- nrow(tbstat);
	  scores <- as.numeric(rownames(tbstat));
	  if (sz > 1)
	  {
		acc <- numeric(nss);
		sen <- numeric(nss);
		csen <- numeric(nscores);
		cspe <- numeric(nscores);
		auc <- numeric(nss);
		ber <- numeric(nss);
		pre <- numeric(nss);
		F1 <- numeric(nss);
		cpre <- numeric(nscores);
		cF1 <- numeric(nscores);
		for (i in 1:nss)
		{
		  bootsample <- datatest[sample(sz,sz,replace = TRUE),];
		  acc[i] <- mean(bootsample[,1] == bootsample[,2]);
		  for (n in 1:nscores)
		  {
			scoresamples <- sum(bootsample[,1] == scores[n])
			allPosPredictions <- sum(bootsample[,2] == scores[n])
			tp <- sum( (bootsample[,1] == scores[n]) & (bootsample[,2] == scores[n]) )

			if (scoresamples > 0)
			{
				csen[n] <- tp/scoresamples;
			}
			else
			{
				csen[n] <- 0;
			}
			if (scoresamples < sz)
			{
				cspe[n] <- sum( (bootsample[,1] != scores[n]) & (bootsample[,2] != scores[n]) )/(sz - scoresamples);
			}
			else
			{
				cspe[n] <- 0;
			}
			if (allPosPredictions > 0)
			{			
				cpre[n] <- tp/allPosPredictions;
			}
			else
			{
				cpre[n] <- NA;
			}
			if (is.na(csen[n]) || is.na(cpre[n]))
			{
				cF1[n] <- NA;
			}
			else
			{
				cF1[n] <- 2.0*csen[n]*cpre[n]/(csen[n]+cpre[n]);
			}
			sen[i] <- sen[i] + csen[n];
			auc[i] <- auc[i] + (csen[n] + cspe[n])/2;
			ber[i] <- ber[i] + 1.0-(csen[n] + cspe[n])/2;
			pre[i] <- pre[i] + cpre[n];
			F1[i] <- F1[i] + cF1[n];
		  }
		  sen[i] <- sen[i]/nscores;
		  auc[i] <- auc[i]/nscores;
		  ber[i] <- ber[i]/nscores;
		  pre[i] <- pre[i]/nscores;
		  F1[i] <- F1[i]/nscores;
		}
		accci <- quantile(acc, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		senci <- quantile(sen, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		aucci <- quantile(auc, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		berci <- quantile(ber, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		F1ci <- quantile(F1, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
		preci <- quantile(pre, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	  }
	  return(list(accci = accci,senci = senci,aucci = aucci,berci = berci,preci = preci,F1ci = F1ci));
}


