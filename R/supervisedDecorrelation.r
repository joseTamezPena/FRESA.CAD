featureDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,loops=c(10,10),thr=0.75,unipvalue=0.05,...)
{

  dataids <- rownames(data)

  if (is.null(refdata))
  {
    refdata <- data;
  }
  refdataids <- rownames(refdata);
  dataAdjusted <- rbind(data,refdata[!(refdataids %in% dataids),]);
#  print(c(nrow(dataAdjusted),ncol(dataAdjusted),nrow(refdata),ncol(refdata)));
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
  colsd <- apply(refdata[,varincluded],2,sd,na.rm = TRUE);
  varincluded <- varincluded[colsd > 0];
#  print(length(varincluded));

  
  wmax = 0.5;
  models <- NULL
  
  cormat <- abs(cor(refdata[,varincluded],method="spearman"))
  diag(cormat) <- 0;
  maxcor <- apply(cormat,2,max)
  varincluded <- names(maxcor)[maxcor >= thr];
  addfeaturematrix <- as.data.frame(matrix(0,nrow=length(varincluded),ncol=length(varincluded)));
  colnames(addfeaturematrix) <- varincluded;
  rownames(addfeaturematrix) <- varincluded;
  
  cormat <- cormat[,varincluded];
  cormat <- cormat[varincluded,];

  while ((addedlist > 0) && (lp < loops[1]))
  {
    lp = lp + 1;
    addedlist <- 0;
#    lastuncorrelatedFetures <- uncorrelatedFetures;
    uncorrelatedFetures <- character();

    maxcor <- apply(cormat,2,max)
    mmaxcor <- max(maxcor);
    topfeat <- colnames(cormat);
    names(topfeat) <- topfeat;
    thr2 <- thr;
    if (thr2 < mmaxcor)
    {
      thr2 <- wmax*mmaxcor + (1.0-wmax)*thr;
    }
    ordcor <- maxcor
    ordcor <- 0.9999*maxcor + 0.0001*apply(cormat,2,mean)
    if (!is.null(Outcome))
    {
      outcomep <- univariate_correlation(refdata[,c(Outcome,varincluded)],Outcome,method="spearman",limit=-1,pvalue=0.20) # the top associated features to the outcome
      selectfeat <- names(outcomep);
      ordcor[selectfeat] <- ordcor[selectfeat] + 0.000001*(1.0 - outcomep); # To sort by associated features to the outcome
    }
    topfeat <- topfeat[maxcor[topfeat] >= thr];
    if (length(topfeat)>0)
    {
      topfeat <- topfeat[order(-ordcor[topfeat])];
#      topfeat <- c(topfeat[topfeat %in% topFeatures],topfeat[!(topfeat %in% topFeatures)]);
      topfeat <- correlated_Remove(refdata[,topfeat],topfeat,thr = thr2);
#      topfeat <- topfeat[order(-ordcor[topfeat])];
      intopfeat <- character();
      for (feat in topfeat)
      {
        corlist <- cormat[,feat];
        corlist <- corlist[corlist >= thr];
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
           if (length(varincluded) > 1000) 
           {
             if ((length(uncorrelatedFetures) %% 100) == 0) cat("|");
           }
           adataframe <- featureAdjustment(dvarlist,
                                          baseModel=feat,
                                          data=dataAdjusted[,c(feat,varlist)],
                                          referenceframe=refdata[,c(feat,varlist)],
                                          pvalue = unipvalue,
                                          ...
                                          );
		  models <- attr(adataframe,"models")
          attr(adataframe,"models") <- NULL
          dataAdjusted[,c(feat,varlist)] <- adataframe;
		  if (length(models) > 0)
		  {
              refdata <- dataAdjusted[refdataids,];

              adjusted <- numeric(length(varlist)) == 1;
              names(adjusted) <- varlist;
			  for (vl in 1:length(models))
			  {
				adjusted[models[[vl]]$feature] <- (models[[vl]]$pval < unipvalue);
#                cat("(",models[[vl]]$feature,":",models[[vl]]$pval,")");
			  }
            adjusted[is.na(adjusted)] <- FALSE;
#            cat("The adjusted :",feat,":");
  #
            varlist <- varlist[adjusted];
#            print(varlist);
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
      if (length(varincluded) > 1000) cat (addedlist,":")
#       cat (addedlist,":")
#      cat ("|")
      if (length(baseFeatures) == 0)
      {
        baseFeatures <- intopfeat;
      }
      cormat <- abs(cor(refdata[,varincluded],method="spearman"))
      diag(cormat) <- 0;
      topFeatures <- unique(c(topFeatures,intopfeat));
    }

    totuncorrelated <- unique(c(totuncorrelated,uncorrelatedFetures));
  }
  dataAdjusted <- dataAdjusted[dataids,];
#  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"featureMatrix") <- addfeaturematrix;
  attr(dataAdjusted,"varincluded") <- varincluded;
  return(dataAdjusted)
}