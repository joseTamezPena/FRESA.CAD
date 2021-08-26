featureDecorrelation <- function(data=NULL,thr=0.80,refdata=NULL,Outcome=NULL,baseFeatures=NULL,unipvalue=0.05,useDeCorr=TRUE,maxLoops=20,wmax=0.0,verbose=FALSE,...)
{

  dataids <- rownames(data)

  if (is.null(refdata))
  {
    refdata <- data;
  }
  refdataids <- rownames(refdata);
  dataAdjusted <- rbind(data,refdata[!(refdataids %in% dataids),]);
  topFeatures <- character()
  correlatedToBase <- character();
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

  
  models <- NULL
  
  cormat <- abs(cor(refdata[,varincluded],method="spearman"))
  diag(cormat) <- 0;
  maxcor <- apply(cormat,2,max)
  varincluded <- names(maxcor)[maxcor >= thr];
  DeCorrmatrix <- NULL;
  if ( !is.null(Outcome) && length(baseFeatures)==0 )
  {
      outcomep <- univariate_correlation(data,Outcome,method="spearman",limit=0,pvalue=0.20,thr = 0.5) # the top associated features to the outcome
      baseFeatures <- names(outcomep);
  }
  
  if (length(varincluded) > 1)
  {
    if (verbose) cat (length(varincluded),":")

    
    DeCorrmatrix <- diag(length(varincluded));
    colnames(DeCorrmatrix) <- varincluded;
    rownames(DeCorrmatrix) <- varincluded;

    cormat <- cormat[,varincluded];
    cormat <- cormat[varincluded,];

    if (verbose) cat (length(baseFeatures),":")


    while ((addedlist > 0) && (lp < maxLoops))
    {
      lp = lp + 1;
      addedlist <- 0;

      maxcor <- apply(cormat,2,max)
      mmaxcor <- max(maxcor);
      topfeat <- colnames(cormat);
      names(topfeat) <- topfeat;
      thr2 <- thr;
      if (verbose) 
      {
          cat("maxCor:",mmaxcor," thr:",thr2,":");
      }
      if (thr2 < mmaxcor)
      {
        thr2 <- wmax*mmaxcor + (1.0-wmax)*thr;
      }
      ordcor <- maxcor
      ordcor <- 0.95*maxcor + 0.05*apply(cormat,2,mean)
      if (length(baseFeatures) > 0)
      {
        ordcor[topfeat %in% baseFeatures] <- ordcor[topfeat %in% baseFeatures] + 1;
      }
      topfeat <- topfeat[maxcor[topfeat] >= thr];
      if (length(topfeat)>0)
      {
        decorrelatedFetureList <- character();
        betamatrix <- diag(length(varincluded));
        colnames(betamatrix) <- varincluded;
        rownames(betamatrix) <- varincluded;
        topfeat <- topfeat[order(-ordcor[topfeat])];
        topfeat <- unique(c(correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE),topfeat[topfeat %in% baseFeatures]));
        intopfeat <- character();
        if (verbose) 
        {
            cat("(",length(topfeat),")");
        }
        for (feat in topfeat)
        {
          corlist <- cormat[,feat];
          corlist <- corlist[corlist >= thr2];
          varlist <- names(corlist)
          varlist <- varlist[!(varlist %in% topfeat)]
          varlist <- varlist[!(varlist %in% decorrelatedFetureList)]
          varlist <- varlist[!(varlist %in% baseFeatures)]
          if (length(varlist) > 0)
          {
             dvarlist <- cbind(varlist,varlist)
             adataframe <- featureAdjustment(dvarlist,
                                            baseModel=feat,
                                            data=dataAdjusted[,c(feat,varlist)],
                                            referenceframe=refdata[,c(feat,varlist)],
                                            pvalue = unipvalue,
                                            ...
                                            );
            models <- attr(adataframe,"models")
            attr(adataframe,"models") <- NULL
            if (length(models) > 0)
            {
                adjusted <- numeric(length(varlist)) == 1;
                names(adjusted) <- varlist;
                for (vl in 1:length(models))
                {
                  if (models[[vl]]$pval <= unipvalue)
                  {
                    adjusted[models[[vl]]$feature] <- TRUE;
                    if (is.null(models[[vl]]$model$coef))
                    {
                      betamatrix[feat,models[[vl]]$feature] <- 1.0;
                      useDeCorr <- FALSE;
                    }
                    else
                    {
                      if (!is.na(models[[vl]]$model$coef[2]))
                      {
                        betamatrix[feat,models[[vl]]$feature] <- -1.0*models[[vl]]$model$coef[2];
                      }
                    }
                  }
                }
              adjusted[is.na(adjusted)] <- FALSE;
              varlist <- varlist[adjusted];
              if (length(varlist) > 0)
              {
                  dataAdjusted[,c(feat,varlist)] <- adataframe[,c(feat,varlist)];
                  refdata[,c(feat,varlist)] <- adataframe[refdataids,c(feat,varlist)];
                  intopfeat <- c(intopfeat,feat);      
                  countf[varlist] <- countf[varlist] + 1;
                  decorrelatedFetureList <- unique(c(decorrelatedFetureList,varlist));
              }
            }
          }
        }
        if (verbose) cat("{");
        addedlist <- length(decorrelatedFetureList) + 1.0*(thr2 > 1.05*thr);
        if (length(decorrelatedFetureList) > 1)
        {
           DeCorrmatrix[,decorrelatedFetureList] <-  DeCorrmatrix %*% betamatrix[,decorrelatedFetureList];
        }
        betamatrix <- NULL;
        if (thr2 > 1.05*thr)
        {
          wmax <- 0.5*wmax;
        }
        else
        {
          wmax <- 0.0;
        }
        if (verbose) cat (addedlist,":")
        if ( (lp == 1) && is.null(Outcome) )
        {
          baseFeatures <- unique(c(baseFeatures,intopfeat));
        }
        if (verbose) cat("}");
        if (verbose) cat("[");
        colsd <- apply(refdata[,varincluded],2,sd,na.rm = TRUE);
        if (sum(colsd==0) > 0)
        {
          if (verbose) cat(sum(colsd==0));
          zerovar <- varincluded[colsd==0];
          for (zcheck in zerovar)
          {
            refdata[,zcheck] <- refdata[,zcheck] + rnorm(nrow(refdata),0,1e-10);
          }
        }
        cormat <- abs(cor(refdata[,varincluded],method="spearman"))
        diag(cormat) <- 0;
        topFeatures <- unique(c(topFeatures,intopfeat));
        if (verbose) cat("]");
      }

    }
    if (useDeCorr)
    {
      dataAdjusted <- data
      if (length(varincluded) > 1)
      {
        dataAdjusted[,varincluded] <- as.matrix(data[,varincluded]) %*% DeCorrmatrix;
        correlatedToBase <- varincluded[varincluded %in% baseFeatures];
        if (length(correlatedToBase) > 1)
        {
          correlatedToBase <- apply(DeCorrmatrix[correlatedToBase,] != 0,2,sum)
          correlatedToBase <- names(correlatedToBase[correlatedToBase > 0])
        }
        else
        {
          if (length(correlatedToBase) == 1)
          {
            correlatedToBase <- varincluded[DeCorrmatrix[correlatedToBase,] != 0]
          }
        }
        correlatedToBase <- correlatedToBase[!(correlatedToBase %in% baseFeatures)];
      }
    }
  }
  dataAdjusted <- dataAdjusted[dataids,];
  
  
  
#  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"DeCorrmatrix") <- DeCorrmatrix;
  attr(dataAdjusted,"varincluded") <- varincluded;
  attr(dataAdjusted,"baseFeatures") <- baseFeatures;
  attr(dataAdjusted,"useDeCorr") <- useDeCorr;
  attr(dataAdjusted,"correlatedToBase") <- correlatedToBase;
  return(dataAdjusted)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"DeCorrmatrix")))
  {
    decorMat <- attr(decorrelatedobject,"DeCorrmatrix")
    testData[,colnames(decorMat)] <- as.matrix(testData[,colnames(decorMat)]) %*% decorMat
  }
  else
  {
    warning("The object does not have a decorrelation Matrix\n")
  }
  return (testData)
}

