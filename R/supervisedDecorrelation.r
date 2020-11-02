supervisedDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,loops=10,unipvalue=0.1,...)
{

  dataAdjusted <- data;
  if (is.null(refdata))
  {
    refdata <- data;
  }
  parameters <- list(...);
  thr <- 0.5;
  if (!is.null(parameters$thr)) thr <- parameters$thr;
  totuncorrelated <- character()
  topFeatures <- character()
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  while ((addedlist > 0) && (lp < loops))
  {
    lp = lp + 1;
    addedlist <- 0;
    topfeat <- univariate_KS(refdata,Outcome,...)
    topFeatures <- unique(c(topFeatures,names(topfeat)));
#    cat(lp,":")
#    print(topFeatures)
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
          addedlist <- addedlist + length(varlist);
          uncorrelatedFetures <- c(uncorrelatedFetures,varlist)
           varlist <- cbind(varlist,varlist)
#           print(c(feat,varlist))
           dataAdjusted <- featureAdjustment(varlist,
                                          baseModel=feat,
                                          data=dataAdjusted,
                                          referenceframe=refdata,
                                          type="LM",
                                          pvalue=unipvalue);
           refdata <- featureAdjustment(varlist,
                                             baseModel=feat,
                                             data=refdata,
                                             referenceframe=refdata,
                                             type="LM",
                                             pvalue=unipvalue);
           
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