steerDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,loops=20,unipvalue=0.05,method=NULL,...)
{

  dataAdjusted <- data;
  if (is.null(refdata))
  {
    refdata <- data;
  }
  parameters <- list(...);
  thr <- 0.375;
  if (!is.null(parameters$thr)) thr <- 0.5*parameters$thr;
  totuncorrelated <- character()
  topFeatures <- character()
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  if (is.null(method)) method ="LM";
  while ((addedlist > 0) && (lp < loops))
  {
    lp = lp + 1;
    addedlist <- 0;
    topfeat <- univariate_KS(refdata,Outcome,...)
    topFeatures <- unique(c(topFeatures,names(topfeat)));
#    cat(lp,":",topFeatures,":")
    lastuncorrelatedFetures <- uncorrelatedFetures;
    uncorrelatedFetures <- character();
    if (length(topFeatures)>0)
    {
      datacor <- refdata[,!(colnames(refdata) %in% Outcome)]
      cormat <- abs(cor(datacor,method="spearman"))
      for (feat in topFeatures)
      {
        corlist <- cormat[feat,]
        corlist <- corlist[corlist >= thr]
        varlist <- names(corlist)
        varlist <- varlist[!(varlist %in% topFeatures)]
        varlist <- varlist[!(varlist %in% uncorrelatedFetures)]
        if (length(varlist) > 0)
        {
           dvarlist <- cbind(varlist,varlist)
#           print(c(feat,varlist))
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
				adjusted[models[[vl]]$feature] <- models[[vl]]$pval < unipvalue;
			  }
		  }
		  varlist <- varlist[adjusted];
          uncorrelatedFetures <- c(uncorrelatedFetures,varlist)
          addedlist <- addedlist + length(varlist);
        }
      }
#      print(uncorrelatedFetures)
      cat (addedlist,":")
      if (length(lastuncorrelatedFetures) == length(uncorrelatedFetures))
      {
        addedlist <- sum(lastuncorrelatedFetures != uncorrelatedFetures);
      }
    }
    totuncorrelated <- unique(c(totuncorrelated,uncorrelatedFetures));
  }
  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- topFeatures;
  attr(dataAdjusted,"Uncorrelated") <- totuncorrelated;
  return(dataAdjusted)
}