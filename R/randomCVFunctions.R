#Filter Methods


FRESAcacheEnv <- new.env();

correlated_Remove <- function(data= NULL,fnames= NULL,thr=0.999)
{
#	cat(thr,":",length(fnames)," ->");
	if (length(fnames)>1)
	{
#		print(fnames);
		colsd <- apply(data[,fnames],2,sd,na.rm = TRUE);
		fnames <- fnames[colsd > 0];
		
		if (length(fnames)>1)
		{
			corm <- abs(cor(data[,fnames],method="spearman"));
			diag(corm) <- 0;
			keep <- numeric(length(fnames));
			for (i in length(fnames):1)
			{
				keep[i] <- (sum(corm[,i] > thr) == 0);
				if (!keep[i])  
				{
					corm[i,] <- 0;
				}
			}
			attributes(fnames) <- list(removed=fnames[keep == 0]);
			fnames <- fnames[keep == 1];
		}
	}
#	cat(length(fnames),"\n");

	return (fnames);
}

correlated_RemoveToLimit <- function(data,unitPvalues,limit=0,thr=0.975,maxLoops=50,minCorr=0.80)
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
				slimit <- limit*nrow(data);
			}
			if (slimit < 2) 
			{
				slimit <- 2;
			}
			if (slimit >= nrow(data))
			{
				slimit = nrow(data) - 1;
			}
			if (length(unitPvalues) > (10*slimit))
			{
				unitPvalues <- unitPvalues[1:(10*slimit)];
			}
			while ( (length(unitPvalues) > slimit) && (ntest < maxLoops) && (cthr > minCorr) )
			{
				unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues),cthr)];
				ntest = ntest + 1;
				cthr = cthr*thr;
			}
			if (length(unitPvalues) > slimit)
			{
				unitPvalues <- unitPvalues[1:slimit];
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
	varlist <-colnames(data);
  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	varlist <- varlist[Outcome != varlist];
	varlist <- cbind(varlist,varlist);
	uniTest <- match.arg(uniTest);
	univ <- ForwardSelection.Model.Bin(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,data,1,type="LOGIT",selectionType=uniTest);
	unitPvalues <- (1.0 - pnorm(univ$base.Zvalues))*2.0;
	names(unitPvalues) <-  varlist[,1];
	unitPvalues <- unitPvalues[order(unitPvalues)];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
#	print(unitPvalues)
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

univariate_residual <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",uniTest=c("Ftest","Binomial","Wilcox","tStudent"),type=c("LM","LOGIT"),limit=0,...,n = 0)
{
	varlist <- colnames(data);
  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	varlist <- varlist[Outcome != varlist];
	varlist <- cbind(varlist,varlist)
	uniTest <- match.arg(uniTest);
	type <- match.arg(type);

	univ <- ForwardSelection.Model.Res(nrow(varlist),1.0,0.0,1,"1",Outcome,varlist,data,1,type=type,testType=uniTest);
	unitPvalues <- (1.0 - pnorm(univ$base.Zvalues))*2.0;
#	print(unitPvalues);
	names(unitPvalues) <- varlist[,1];
	unitPvalues <- unitPvalues[order(unitPvalues)];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	return(unitPvalues);
}

univariate_KS <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	varlist <-colnames(data);
  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;

	for (j in varlist) 
	{
		 pval <- ks.test(control[,j],case[,j],na.action = na.exclude)$p.value; 
		 if (inherits(pval, "try-error")) {pval <- 1.0;}
		 if (is.null(pval)) { pval <- 1.0; }
		 if (is.na(pval)) {pval <- 1.0;}
		 unitPvalues[j] <- pval;
	}
  	unitPvalues <- unitPvalues[order(unitPvalues)];
  	unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues < pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
   return(unitPvalues);
}


univariate_Wilcoxon <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	varlist <-colnames(data);
  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;

	for (j in varlist) 
	{
		 pval <- wilcox.test(control[,j],case[,j],na.action = na.exclude)$p.value; 
		 if (inherits(pval, "try-error")) {pval <- 1.0;}
		 if (is.null(pval)) { pval <- 1.0; }
		 if (is.na(pval)) {pval <- 1.0;}
		 unitPvalues[j] <- pval;
	}
  
	unitPvalues <- unitPvalues[order(unitPvalues)];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	
   return(unitPvalues);
}

univariate_tstudent <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",limit=0,...,n = 0)
{

	varlist <-colnames(data);
  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	case <- subset(data,get(Outcome) == 1);
	control <- subset(data,get(Outcome) == 0);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;
	for (j in varlist) 
	{
		pval <- try(t.test(control[,j],case[,j],na.action = na.exclude)$p.value);
		if (inherits(pval, "try-error")) {pval <- 1.0;}
		if (is.null(pval)) { pval <- 1.0; }
		if (is.na(pval)) {pval <- 1.0;}
		unitPvalues[j] <-  pval;	
	}
	unitPvalues <- unitPvalues[order(unitPvalues)];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

univariate_correlation <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", method = "kendall",limit=0,...,n = 0)
{

  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
	varlist <-colnames(data);
	varlist <- varlist[Outcome != varlist];
	unitPvalues <- numeric(length(varlist));
	names(unitPvalues) <- varlist;

	for (j in varlist) 
	{
		pval <- try(cor.test(data[,Outcome],data[,j],na.action = na.exclude,method = method)$p.value);
		if (inherits(pval, "try-error")) {pval <- 1.0;}
		if (is.null(pval)) { pval <- 1.0; }
		if (is.na(pval)) {pval <- 1.0;}
		unitPvalues[j] <-  pval;
	}
  	unitPvalues <- unitPvalues[order(unitPvalues)];
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
#		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
   return(unitPvalues);
}

mRMR.classic_FRESA <- function(data=NULL, Outcome=NULL,feature_count=0,...)
{
	if (!requireNamespace("mRMRe", quietly = TRUE)) {
	   install.packages("mRMRe", dependencies = TRUE)
	} 

  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));
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

  if (class(data[,Outcome]) == "factor") data[,Outcome] <- as.numeric(as.character(data[,Outcome]));

#  tpvalue <- max(c(sqrt(pvalue)/2.0,pvalue));
  tpvalue <- min(c(2.0*pvalue,0.499));
  
  pvallist <- list();
  pvaltest <- univariate_Logit(data,Outcome,pvalue=tpvalue,limit=-1,uniTest="zNRI",adjustMethod=adjustMethod);
  pvallist$LogitNRI <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- c(1:length(pvaltest))

  wilcxf <- univariate_Wilcoxon(data,Outcome,pvalue=tpvalue,limit=-1,adjustMethod=adjustMethod);
  pvallist$Wilcox <- wilcxf;
  varcount[names(wilcxf)] <- varcount[names(wilcxf)]+1;
  rankVar[names(wilcxf)] <- rankVar[names(wilcxf)] + c(1:length(wilcxf));

  features <- intersect(names(pvaltest),names(wilcxf));
  both <- pmin(pvaltest[features],wilcxf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],wilcxf[!(names(wilcxf) %in% features)],both);

  pvaltest <- univariate_Logit(data,Outcome,pvalue=tpvalue,limit=-1,uniTest="zIDI",adjustMethod=adjustMethod);
  pvallist$LogitIDI <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  pvaltest <- univariate_residual(data,Outcome,pvalue=tpvalue,limit=-1,uniTest="Wilcox",type="LOGIT",adjustMethod=adjustMethod);
  pvallist$LogitWilcox <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  pvaltest <- univariate_residual(data,Outcome,pvalue=tpvalue,limit=-1,uniTest="Ftest",type="LOGIT",adjustMethod=adjustMethod);
  pvallist$LogitFtest <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  pvaltest <- univariate_residual(data,Outcome,pvalue=tpvalue,limit=-1,uniTest="Binomial",type="LOGIT",adjustMethod=adjustMethod);
  pvallist$LogitBin <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  pvaltest <- univariate_correlation(data,Outcome,pvalue=tpvalue,limit=-1, method = "kendall",adjustMethod=adjustMethod)
  pvallist$CorKendall <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  pvaltest <- univariate_correlation(data,Outcome,pvalue=tpvalue,limit=-1,..., method = "spearman",adjustMethod=adjustMethod)
  pvallist$CorSpearman <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);
  
  pvaltest <- univariate_KS(data,Outcome,pvalue=tpvalue,limit=-1,adjustMethod=adjustMethod)
  pvallist$KS <- pvaltest;
  varcount[names(pvaltest)] <- varcount[names(pvaltest)]+1;
  rankVar[names(pvaltest)] <- rankVar[names(pvaltest)] + c(1:length(pvaltest));
  unadjustedKS <- attr(pvaltest,"Unadjusted");
  features <- intersect(names(pvaltest),names(allf));
  both <- pmin(pvaltest[features],allf[features]);
  allf <- c(pvaltest[!(names(pvaltest) %in% features)],allf[!(names(allf) %in% features)],both);

  limitmrmr <- limit
  if ((limitmrmr < 1) && (limitmrmr > 0))
  {
    limitmrmr = as.integer(nrow(data)*limitmrmr + 0.5);
  }
  if (limit == 0) limit = length(allf);
  if ( limitmrmr < 2 ) limitmrmr = 2;

  mRMRf <- mRMR.classic_FRESA(data,Outcome,feature_count = limitmrmr)

  mRMRf <- mRMRf[mRMRf > 0];
  varcount[names(mRMRf)] <- varcount[names(mRMRf)] + 1.0;
  rankVar[names(mRMRf)] <- rankVar[names(mRMRf)] + c(1:length(mRMRf))/length(mRMRf);
  mRMRf <- unadjustedKS[names(mRMRf)];
  mRMRf <- mRMRf[mRMRf <= 0.05];
#  mRMRf <- mRMRf/(1 + varcount[names(mRMRf)]);
#  mRMRf <- p.adjust(mRMRf,adjustMethod,n=length(unadjustedKS))
  pvallist$mRMR <- mRMRf;
  if (length(mRMRf)>0)
  {
    missing <- !(names(mRMRf) %in% names(allf))
    allf <- c(allf,mRMRf[missing])
  }

  varcount <- varcount[names(allf)];
  rankVar <- rankVar[names(allf)]/varcount;
  allf <- allf[order(rankVar + allf - ncol(data)*varcount)];
  varcount <- varcount[names(allf)];
  rankVar <- rankVar[names(allf)];

  top <- allf[1];

  padjs <- sqrt(varcount);
  allf <- allf/padjs;
  allf <- allf[allf <= pvalue];
  if (length(allf) > 1) 
  {
  	  parameters <- list(...);
	  thr = 0.975;
	  if (!is.null(parameters$thr))
	  {
	  	thr <- parameters$thr;
	  }
	  allf <- allf[correlated_Remove(data,names(allf),thr)];
	  allf <- correlated_RemoveToLimit(data,allf,limit=limit,...);
  }
  else 
  {
	allf <- top;
  }
#  allf <- unadjustedKS[names(allf)];
  rankVar <- rankVar[names(allf)];

  attr(allf,"varcount") <- varcount;
  attr(allf,"UnadjustedKS") <- unadjustedKS;
  attr(allf,"Pvalues") <- pvallist;
  attr(allf,"rankVar") <- rankVar;
  return (allf);
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
			if (scoresamples > 0)
			{
				csen[n] <- sum( (bootsample[,1] == scores[n]) & (bootsample[,2] == scores[n]) )/scoresamples;
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
			tp <- sum( (bootsample[,2] == scores[n]) );
			if (tp > 0)
			{			
				cpre[n] <- sum( (bootsample[,1] == scores[n]) & (bootsample[,2] == scores[n]) )/tp ;
			}
			else
			{
				cpre[n] <- NA;
			}
			if (is.na(csen) || is.na(cpre))
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


