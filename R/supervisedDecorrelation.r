featureDecorrelation <- function(data=NULL, Outcome=NULL,refdata=NULL,baseFeatures=NULL,loops=20,thr=0.80,unipvalue=0.05,useWhite=TRUE,...)
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
  if (is.null(baseFeatures))
  {
    baseFeatures <- character()
  }
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  countf <- numeric(ncol(refdata));
  names(countf) <- colnames(refdata);
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


  if (length(varincluded) > 1000) cat (length(varincluded),":")

  if ( !is.null(Outcome) && length(baseFeatures)==0 )
  {
      outcomep <- univariate_correlation(refdata[,c(Outcome,varincluded)],Outcome,method="spearman",limit=0,pvalue=0.20,thr=thr) # the top associated features to the outcome
      baseFeatures <- names(outcomep);
#      print(baseFeatures);
  }
  
  whiteningmatrix <- diag(length(varincluded));
  colnames(whiteningmatrix) <- varincluded;
  rownames(whiteningmatrix) <- varincluded;

  cormat <- cormat[,varincluded];
  cormat <- cormat[varincluded,];

  if (length(varincluded) > 1000) cat (length(baseFeatures),":")


  while ((addedlist > 0) && (lp < loops))
  {
    lp = lp + 1;
    addedlist <- 0;

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
    topfeat <- topfeat[maxcor[topfeat] >= thr];
    if (length(topfeat)>0)
    {
      uncorrelatedFetures <- character();
      betamatrix <- diag(length(varincluded));
      colnames(betamatrix) <- varincluded;
      rownames(betamatrix) <- varincluded;
      topfeat <- topfeat[order(-ordcor[topfeat])];
      topfeat <- correlated_Remove(cormat,topfeat,thr = thr2,isDataCorMatrix=TRUE);
      intopfeat <- character();
      if (length(varincluded) > 1000) 
      {
          cat("(",length(topfeat),")");
      }
      for (feat in topfeat)
      {
#        if (length(varincluded) > 1000) cat(".");
        corlist <- cormat[,feat];
        corlist <- corlist[corlist >= thr];
#        cat(lp,":",feat,":");
#        print(corlist)
        varlist <- names(corlist)
        varlist <- varlist[!(varlist %in% topfeat)]
        varlist <- varlist[!(varlist %in% uncorrelatedFetures)]
        varlist <- varlist[!(varlist %in% baseFeatures)]
        if (length(varlist) > 0)
        {
           dvarlist <- cbind(varlist,varlist)
#          print(corlist[varlist])

#          if (length(varincluded) > 1000) 
#          {
#              if ((length(uncorrelatedFetures) %% 10) == 0) cat("<",length(varlist),"\\");
#          }
           adataframe <- featureAdjustment(dvarlist,
                                          baseModel=feat,
                                          data=dataAdjusted[,c(feat,varlist)],
                                          referenceframe=refdata[,c(feat,varlist)],
                                          pvalue = unipvalue,
                                          ...
                                          );
		  models <- attr(adataframe,"models")
          attr(adataframe,"models") <- NULL
#          if (length(varincluded) > 1000) 
#          {
#             if ((length(uncorrelatedFetures) %% 10) == 0) cat("|");
#          }
		  if (length(models) > 0)
		  {
              adjusted <- numeric(length(varlist)) == 1;
              names(adjusted) <- varlist;
			  for (vl in 1:length(models))
			  {
                if (models[[vl]]$pval < unipvalue)
                {
                  adjusted[models[[vl]]$feature] <- TRUE;
                  if (is.null(models[[vl]]$model$coef))
                  {
                    betamatrix[feat,models[[vl]]$feature] <- 1.0;
                    useWhite <- FALSE;
                  }
                  else
                  {
                    if (!is.na(models[[vl]]$model$coef[2]))
                    {
                      betamatrix[feat,models[[vl]]$feature] <- -1.0*models[[vl]]$model$coef[2];
                    }
                  }
                }
#                cat("(",models[[vl]]$feature,":",models[[vl]]$pval,")");
			  }
            adjusted[is.na(adjusted)] <- FALSE;
#            cat("The adjusted :",feat,":");
  #
            varlist <- varlist[adjusted];
#            print(varlist);
            if (length(varlist) > 0)
            {
                dataAdjusted[,c(feat,varlist)] <- adataframe[,c(feat,varlist)];
                refdata[,c(feat,varlist)] <- adataframe[refdataids,c(feat,varlist)];
                intopfeat <- c(intopfeat,feat);      
                countf[varlist] <- countf[varlist] + 1;
                uncorrelatedFetures <- unique(c(uncorrelatedFetures,varlist));
            }
#             if (length(varincluded) > 1000) 
#             {
#               if ((length(uncorrelatedFetures) %% 10) == 0) cat("/");
#             }
		  }
#          if (length(varincluded) > 1000) 
#          {
#              if ((length(uncorrelatedFetures) %% 10) == 0) cat(">");
#          }
        }
      }
      if (length(varincluded) > 1000) cat("{");
      addedlist <- length(uncorrelatedFetures) + 1.0*(thr2 > 1.005*thr);
      if (length(uncorrelatedFetures) > 0)
      {
         whiteningmatrix[,uncorrelatedFetures] <-  whiteningmatrix %*% betamatrix[,uncorrelatedFetures];
      }
      betamatrix <- NULL;
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
      if (length(varincluded) > 1000) cat("}");
      if (length(varincluded) > 1000) cat("[");
      cormat <- abs(cor(refdata[,varincluded],method="spearman"))
      diag(cormat) <- 0;
      topFeatures <- unique(c(topFeatures,intopfeat));
      if (length(varincluded) > 1000) cat("]");
    }

    totuncorrelated <- unique(c(totuncorrelated,uncorrelatedFetures));
  }
  if (useWhite)
  {
    dataAdjusted <- data
    dataAdjusted[,varincluded] <- as.matrix(data[,varincluded]) %*% whiteningmatrix;
  }
  else
  {
    dataAdjusted <- dataAdjusted[dataids,];
  }
  
  
#  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"whiteningmatrix") <- whiteningmatrix;
  attr(dataAdjusted,"varincluded") <- varincluded;
  attr(dataAdjusted,"baseFeatures") <- baseFeatures;
  attr(dataAdjusted,"useWhite") <- useWhite;
  return(dataAdjusted)
}