featureDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,loops=c(20,10),thr=0.75,unipvalue=0.05,...)
{

  dataids <- rownames(data)

  if (is.null(refdata))
  {
    refdata <- data;
  }
  refdataids <- rownames(refdata);
  dataAdjusted <- rbind(data,refdata[!(refdataids %in% dataids),]);
  totuncorrelated <- character()
  topFeatures <- character()
  baseFeatures <- character()
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  countf <- numeric(ncol(refdata));
  names(countf) <- colnames(refdata);
  tsum <- 10;
  if (length(loops) > 1) tsum <- loops[2];
  varincluded <- colnames(refdata)[!(colnames(refdata) %in% Outcome)];
  addfeaturematrix <- as.data.frame(matrix(0,nrow=length(varincluded),ncol=length(varincluded)));
  colnames(addfeaturematrix) <- varincluded;
  rownames(addfeaturematrix) <- varincluded;
  wmax = 0.5;
  

  while ((addedlist > 0) && (lp < loops[1]))
  {
    lp = lp + 1;
    addedlist <- 0;
    refdata <- dataAdjusted[refdataids,];
    lastuncorrelatedFetures <- uncorrelatedFetures;
    uncorrelatedFetures <- character();

	cormat <- abs(cor(refdata[,varincluded],method="spearman"))
    diag(cormat) <- 0;
    maxcor <- apply(cormat,2,max)
    mmaxcor <- max(maxcor);
    topfeat <- colnames(cormat);
    thr2 <- thr;    
    if (thr2 < mmaxcor)
    {
      thr2 <- wmax*mmaxcor + (1.0-wmax)*thr;
    }
    ordcor <- 0.9999*maxcor + 0.0001*apply(cormat,2,mean)
    topfeat <- topfeat[order(-ordcor)];
    names(topfeat) <- topfeat;
    topfeat <- c(topfeat[topFeatures],topfeat[!(topfeat %in% topFeatures)]);
    if (!is.null(Outcome))
    {
      topfeat <- names(univariate_correlation(refdata,Outcome,method="spearman",limit=-1,pvalue=0.45));
      topfeat <- c(topfeat[topfeat %in% topFeatures],topfeat[!(topfeat %in% topFeatures)]);
    }
    topfeat <- topfeat[maxcor[topfeat] >= thr];
    if (length(topfeat)>0)
    {
      topfeat <- correlated_Remove(refdata,topfeat,thr = thr2);
      topfeat <- topfeat[order(-ordcor[topfeat])];
      intopfeat <- character();
      for (feat in topfeat)
      {
        corlist <- cormat[,feat];
        corlist <- corlist[corlist >= thr]
#        cat(lp,":",feat,":");
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
                                          pvalue = unipvalue,
                                          ...
                                          );
#          refdata <- dataAdjusted[refdataids,];
		  models <- attr(dataAdjusted,"models")
		  if (length(models) > 0)
		  {
              adjusted <- numeric(length(varlist)) == 1;
              names(adjusted) <- varlist;
			  for (vl in 1:length(models))
			  {
				adjusted[models[[vl]]$feature] <- (models[[vl]]$pval < unipvalue);
#                cat("(",models[[vl]]$feature,":",models[[vl]]$pval,")");
			  }
  #          cat("The adjusted :",feat,":");
  #
  #          print(varlist);
            varlist <- varlist[adjusted];
            if (length(varlist) > 0)
            {
                intopfeat <- c(intopfeat,feat);      
                addfeaturematrix[feat,varlist] <- addfeaturematrix[feat,varlist] + 1;
#                print(varlist);
                countf[varlist] <- countf[varlist] + 1;
                uncorrelatedFetures <- unique(c(uncorrelatedFetures,varlist));
            }
		  }
        }
        addedlist <- length(uncorrelatedFetures) + 1.0*(thr2 > 1.001*thr);
      }
      if (thr2 > 1.001*thr)
      {
        wmax <- 0.5*wmax;
      }
#      print(uncorrelatedFetures)
#      if (length(lastuncorrelatedFetures) == length(uncorrelatedFetures))
#      {
#        addedlist <- sum(lastuncorrelatedFetures != uncorrelatedFetures);
#      }
#      cat (addedlist,":")
      if (length(baseFeatures) == 0)
      {
        baseFeatures <- intopfeat;
      }
      topFeatures <- unique(c(topFeatures,intopfeat));
    }
    totuncorrelated <- unique(c(totuncorrelated,uncorrelatedFetures));
  }
  dataAdjusted <- dataAdjusted[dataids,];
#  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"featureMatrix") <- addfeaturematrix;
  return(dataAdjusted)
}