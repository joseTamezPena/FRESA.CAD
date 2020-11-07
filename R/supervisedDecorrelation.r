featureDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,loops=c(20,10),unipvalue=0.05,method=NULL,...)
{

  dataAdjusted <- data;
  if (is.null(refdata))
  {
    refdata <- data;
  }
  parameters <- list(...);
  thr <- 0.75;
  if (!is.null(parameters$thr))
  {
    thr <- parameters$thr;
  }
  totuncorrelated <- character()
  topFeatures <- character()
  baseFeatures <- character()
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  if (is.null(method)) method ="LM";
  countf <- numeric(ncol(refdata));
  names(countf) <- colnames(refdata);
  tsum <- 10;
  if (length(loops) > 1) tsum <- loops[2];
  varincluded <- !(colnames(refdata) %in% Outcome);
  addfeaturematrix <- as.data.frame(matrix(0,nrow=length(varincluded),ncol=length(varincluded)));
  colnames(addfeaturematrix) <- colnames(refdata)[varincluded];
  rownames(addfeaturematrix) <- colnames(addfeaturematrix);
  unipvalue = 2.0*unipvalue/(nrow(data)-1);
  while ((addedlist > 0) && (lp < loops[1]))
  {
    lp = lp + 1;
    addedlist <- 0;
	cormat <- abs(cor(refdata[,varincluded],method="spearman"))
    diag(cormat) <- 0;
    maxcor <- apply(cormat,2,max)
    mmaxcor <- max(maxcor);
    topfeat <- colnames(cormat);
    thr2 <- thr;    
    if (is.null(Outcome))
    {
      ordcor <- 0.99*maxcor + 0.01*apply(cormat,2,mean)
      topfeat <- topfeat[order(-ordcor)];
      names(topfeat) <- topfeat;
      topfeat <- c(topfeat[topFeatures],topfeat[!(topfeat %in% topFeatures)]);

      if (thr2 < mmaxcor)
      {
        thr2 <- 0.25*mmaxcor + 0.75*thr;
      }
      topfeat <- topfeat[maxcor[topfeat] >= thr];
      if (length(topfeat) > 0)
      {
          topfeat <- correlated_Remove(refdata,topfeat,thr = thr2);
          topfeat <- topfeat[order(-ordcor[topfeat])];
      }
    }
    else
    {
      topfeat <- names(univariate_KS(refdata,Outcome,...))
    }
    lastuncorrelatedFetures <- uncorrelatedFetures;
    uncorrelatedFetures <- character();
    if (length(topfeat)>0)
    {
      intopfeat <- character();
      for (feat in topfeat)
      {
        corlist <- cormat[,feat];
        corlist <- corlist[corlist >= thr]
#        cat(feat,":");
#        print(corlist)
        varlist <- names(corlist)
        varlist <- varlist[!(varlist %in% baseFeatures)]
        varlist <- varlist[!(varlist %in% topfeat)]
        varlist <- varlist[!(varlist %in% uncorrelatedFetures)]
		varlist <- varlist[countf[varlist] < tsum]
        if (length(varlist) > 0)
        {
           dvarlist <- cbind(varlist,varlist)
#          print(corlist[varlist])
           dataAdjusted <- featureAdjustment(dvarlist,
                                          baseModel=feat,
                                          data=dataAdjusted,
                                          referenceframe=refdata,
                                          type=method,
                                          pvalue=unipvalue);
           refdata <- featureAdjustment(dvarlist,
                                             baseModel=feat,
                                             data=refdata,
                                             referenceframe=refdata,
                                             type=method,
                                             pvalue=unipvalue);
		  adjusted <- numeric(length(varlist)) == 1;
		  names(adjusted) <- varlist;
		  models <- attr(refdata,"models")
		  if (length(models) > 0)
		  {
			  for (vl in 1:length(models))
			  {
				adjusted[models[[vl]]$feature] <- (models[[vl]]$pval < unipvalue);
			  }
              intopfeat <- c(intopfeat,feat);      
		  }
#          cat("The adjusted :",feat,":");
#          print(varlist);
		  varlist <- varlist[adjusted];
          addfeaturematrix[feat,varlist] <- addfeaturematrix[feat,varlist] + 1;
#         print(varlist);
		  countf[varlist] <- countf[varlist] + 1;
          uncorrelatedFetures <- unique(c(uncorrelatedFetures,varlist));
          addedlist <- length(uncorrelatedFetures) + 1.0*(thr2 > 1.05*thr);
        }
      }
#      print(uncorrelatedFetures)
      if (length(lastuncorrelatedFetures) == length(uncorrelatedFetures))
      {
        addedlist <- sum(lastuncorrelatedFetures != uncorrelatedFetures);
      }
      cat (addedlist,":")
      if (length(baseFeatures) == 0)
      {
        baseFeatures <- intopfeat;
      }
      topFeatures <- unique(c(topFeatures,intopfeat));
    }
    totuncorrelated <- unique(c(totuncorrelated,uncorrelatedFetures));
  }
  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"featureMatrix") <- addfeaturematrix;
  return(dataAdjusted)
}