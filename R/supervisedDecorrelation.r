featureDecorrelation <- function(data=NULL,
                                  thr=0.80,
                                  refdata=NULL,
                                  Outcome=NULL,
                                  baseFeatures=NULL,
                                  unipvalue=0.05,
                                  useDeCorr=TRUE,
                                  maxLoops=20,
                                  verbose=FALSE,
                                  method=c("spearman","pearson","kendall","fast"),
                                  ...)
{

	if (!requireNamespace("Rfast", quietly = TRUE)) {
		install.packages("Rfast", dependencies = TRUE)
	} 

  method <- match.arg(method);
  useFastCor <- method=="fast";
  if (useFastCor) 
  {
    useDeCorr <- TRUE;
    method <- "pearson";
  }

getAllBetaCoefficients <- function(feat,varlist=NULL)
{
  
  featVector <- refdata[,feat]
#  print(featVector[1:10])
  getBetaCoefficient <- function(x)
  {
    betaCoef <- 0;
    if (class(x) != "factor")
    {
      modellm <- try(lm(x~basefeature,data=as.data.frame(cbind(x=x,basefeature=featVector)),model = FALSE,na.action=na.exclude))
      if (!inherits(modellm, "try-error"))
      {	
        f <- summary(modellm)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=FALSE);
        if (is.na(p)) p <- 1.0;
        if (p < unipvalue)
        {
          betaCoef <- modellm$coef[2]
        }
      }
    }
    return (betaCoef)
  }
#  cat(feat,"(",length(varlist),")","\n")
  allBetaCoef <- apply(as.data.frame(refdata[,varlist]),2,getBetaCoefficient)
  allBetaCoef[is.na(allBetaCoef)] <- 0;
  allBetaCoef[is.nan(allBetaCoef)] <- 0;
  allBetaCoef[is.infinite(allBetaCoef)] <- 0;
  names(allBetaCoef) <- varlist
  return (allBetaCoef)
}



  dataids <- rownames(data)

  if (is.null(refdata))
  {
    refdata <- data;
  }
  refdataids <- rownames(refdata);
  dataAdjusted <- as.data.frame(rbind(data,refdata[!(refdataids %in% dataids),]));
  topFeatures <- character()
  correlatedToBase <- character();
  if (is.null(baseFeatures))
  {
    baseFeatures <- character()
  }
  AbaseFeatures <- character()
  addedlist <- 1;
  lp = 0;
  uncorrelatedFetures <- character();
  countf <- numeric(ncol(refdata));
  names(countf) <- colnames(refdata);
  if (!is.null(Outcome))
  {
    varincluded <- colnames(refdata)[!(colnames(refdata) %in% Outcome)];
  }
  else
  {
    varincluded <- colnames(refdata);
  }
  colsd <- apply(refdata[,varincluded],2,sd,na.rm = TRUE);
  varincluded <- varincluded[colsd > 0];
  totFeatures <- length(varincluded);

  
  models <- NULL
  
  if (useFastCor)
  {
    cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded])))
  }else
  {
    cormat <- abs(cor(refdata[,varincluded],method=method))
  }
  diag(cormat) <- 0;
  maxcor <- apply(cormat,2,max)
  totcorr <- sum(cormat >= thr); 
  varincluded <- names(maxcor)[maxcor >= 0.5*thr];
  DeCorrmatrix <- NULL;
  outcomep <- numeric();
  if ( !is.null(Outcome) && length(baseFeatures)==0 )
  {
      outcomep <- univariate_correlation(data,Outcome,method=method,limit=0,pvalue=0.20,thr = 0.99*thr) # the top associated features to the outcome
      baseFeatures <- names(outcomep);
  }
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  if (length(varincluded) > 1)
  {
    unipvalue <- min(0.05,4.0*unipvalue*sqrt(totcorr)/totFeatures); ## Adjusting for false association
    bvarincluded <- character();
    baseIncluded <- character();
    if (length(baseFeatures) > 0)
    {    
      baseIncluded <- varincluded[varincluded %in% baseFeatures];
      names(baseIncluded) <- baseIncluded;
      if (length(baseIncluded) > 0)
      {
        bcormat <-  cormat[baseIncluded,];
        bmaxcor <- bcormat;
        if (length(baseIncluded) > 1)
        {
          bmaxcor <- apply(bcormat,2,max)
        }
        bvarincluded <- names(bmaxcor)[bmaxcor >= thr];
        bcormat <- NULL;
      } 
    }
#    if ((length(varincluded) > 1) || is.null(Outcome) )
    {
      if (verbose) cat ("\n Included:",length(varincluded),", Uni p:",unipvalue,", Base:",length(baseFeatures),", In Included:",length(baseIncluded),", Base Cor:",length(bvarincluded))
      
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;

      cormat <- cormat[,varincluded];
      cormat <- cormat[varincluded,];

      while ((addedlist > 0) && (lp < maxLoops))
      {
        lp = lp + 1;
        addedlist <- 0;

        maxcor <- apply(cormat,2,max)
        mmaxcor <- max(maxcor);
        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        ordcor <- maxcor
        ordcor <- 0.99*maxcor + 0.009*apply(1*(cormat >= thr),2,mean) + 0.0001*apply(cormat,2,mean)
        if (length(baseFeatures) > 0)
        {
          ordcor[baseIncluded] <- ordcor[baseIncluded] + 2;
        }
        if (length(AbaseFeatures) > 0)
        {
          ordcor[AbaseFeatures] <- ordcor[AbaseFeatures] + 1;
        }
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (length(topfeat)>0)
        {
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          topfeat <- topfeat[order(-ordcor[topfeat])];
#          topfeat <- unique(c(correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE),topfeat[topfeat %in% baseIncluded]));
          topfeat <- unique(c(baseIncluded[baseIncluded %in% topfeat],
                              correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE),
                              AbaseFeatures[AbaseFeatures %in% topfeat]));
          intopfeat <- character();
          if (verbose) 
          {
            cat(", Top:",length(topfeat));
#            if (length(topfeat)<10) print(topfeat)
          }
          for (feat in topfeat)
          {
            corlist <- cormat[,feat];
            corlist <- corlist[corlist >= thr];
            varlist <- names(corlist)
#            if (verbose) 
#            {
#              if (length(topfeat)<10) print(varlist)
#            }
            varlist <- varlist[!(varlist %in% topfeat)]
            varlist <- varlist[!(varlist %in% decorrelatedFetureList)]
            varlist <- varlist[!(varlist %in% baseIncluded)]
            varlist <- varlist[!(varlist %in% AbaseFeatures)]
            if ((lp < 4) && !(feat %in% baseIncluded))
            {
              varlist <- varlist[!(varlist %in% bvarincluded)]
            }
            varlist <- names(corlist)
            # if (verbose) 
            # {
              # if (length(topfeat) < 10) print(varlist)
            # }
            if (length(varlist) > 0)
            {
               if (useFastCor)
               {
                  prebetas <- getAllBetaCoefficients(feat,varlist);
                  betamatrix[feat,varlist] <- -1.0*prebetas;
                  varlist <- varlist[prebetas != 0];
                  if (length(varlist) > 0)
                  {
                    intopfeat <- c(intopfeat,feat);
                    countf[varlist] <- countf[varlist] + 1;
                    decorrelatedFetureList <- unique(c(decorrelatedFetureList,varlist));
                  }
               }
               else
               {
                 dvarlist <- cbind(c(varlist),c(varlist))
                 adataframe <- featureAdjustment(dvarlist,
                                                baseFormula=feat,
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
                      adjusted[models[[vl]]$feature] <- (models[[vl]]$pval < unipvalue);
                      # if (verbose) 
                      # {
                        # if ((length(topfeat)<10)&& (vl < 3)) 
                        # {
                          # cat("Base: ",feat," Cor: ",models[[vl]]$feature," Pvalue=",models[[vl]]$pval,":",models[[vl]]$model$coef[2],"\n")
                        # }
                      # }
                      if (adjusted[models[[vl]]$feature])
                      {
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
                          else
                          {
                              adjusted[models[[vl]]$feature] <- FALSE;
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
            if (verbose && (length(intopfeat) %% 100 == 99)) cat(".")
          }
          addedlist <- length(decorrelatedFetureList);
          if (verbose) cat(",(",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
          if (addedlist > 0)
          {
             mbetas <- c(intopfeat,decorrelatedFetureList);
             totused <- unique(c(totused,mbetas));
             totalpha <- unique(c(totalpha,intopfeat));
             allused <- varincluded[varincluded %in% totused];
             colused <- varincluded[varincluded %in% decorrelatedFetureList];
             alphaused <- varincluded[varincluded %in% intopfeat];
             totalphaused <- varincluded[varincluded %in% totalpha];
#             DeCorrmatrix[,decorrelatedFetureList] <-  DeCorrmatrix %*% as.matrix(betamatrix[,decorrelatedFetureList]);
              if ((addedlist == 1) && (length(totalphaused) == 1))
              {
                DeCorrmatrix[totalphaused,colused] <-  sum(as.numeric(DeCorrmatrix[totalphaused,allused])*
                                                               as.numeric(betamatrix[allused,colused]));
              }
              else
              {
                if (length(totalphaused) > 1)
                {
                  DeCorrmatrix[totalphaused,colused] <-  Rfast::mat.mult( as.matrix(DeCorrmatrix[totalphaused,allused]),
                                                                                as.matrix(betamatrix[allused,colused])
                                                                                );
                }
                else
                {
                   cat("..no Fast..[",length(totalphaused),"][",length(allused),"][",length(colused),"]")
                   mtx1 <- t(as.matrix(DeCorrmatrix[totalphaused,allused]))
                   mtx2 <- as.matrix(betamatrix[allused,colused])
#                   cat("(",nrow(mtx1),",",ncol(mtx1),")(",nrow(mtx2),",",ncol(mtx2),")");
                  DeCorrmatrix[totalphaused,colused] <- mtx1  %*% mtx2;
                }
              }
#              { 
#                cat("..no Fast..")
#                  cat("*")
#                  DeCorrmatrix[,decorrelatedFetureList] <-  DeCorrmatrix %*% as.matrix(betamatrix[,decorrelatedFetureList]);
#                }
#                DeCorrmatrix[totalphaused,colused] <-  sum(as.numeric(DeCorrmatrix[totalphaused,allused])*
#                                                               as.numeric(betamatrix[allused,colused]));
#              }
             if (useFastCor)
             {
                if (verbose) cat("|")
#                refdata[,varincluded] <- as.matrix(dataAdjusted[refdataids,varincluded]) %*% DeCorrmatrix
                if (( length(colused) > 1) || (length(alphaused) > 1 ))
                {
                  refdata[,colused] <- refdata[,colused] + Rfast::mat.mult(as.matrix(refdata[,alphaused]),as.matrix(betamatrix[alphaused,colused]))
                }
                else
                {
                  refdata[,colused] <- refdata[,colused] + sum(as.numeric(refdata[,alphaused])*as.numeric(betamatrix[alphaused,colused]))
                }
             }
          }
          if (verbose) cat(">")
#          betamatrix <- NULL;
          if (lp == 1)
          {
            AbaseFeatures <- intopfeat;
            names(AbaseFeatures) <- AbaseFeatures;
          }
          colsd <- apply(refdata[,varincluded],2,sd,na.rm = TRUE);
          if (sum(colsd==0) > 0)
          {
            zerovar <- varincluded[colsd==0];
            for (zcheck in zerovar)
            {
              refdata[,zcheck] <- refdata[,zcheck] + rnorm(nrow(refdata),0,1e-10);
            }
          }
          topFeatures <- unique(c(topFeatures,intopfeat));
          if (useFastCor)
          {
            cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded])))
          }
          else
          {
            cormat <- abs(cor(refdata[,varincluded],method=method))
          }
          diag(cormat) <- 0;
          if (verbose) 
          {
            cat (", Tot Used:",length(totused),", Added:",addedlist,", Zero Std:",sum(colsd==0),", Max Cor:",max(cormat));
#            if (lp > 5)
#           {
#              cat(",",decorrelatedFetureList[1],":",intopfeat[1]);
#            }
          }
          if ( (length(intopfeat) == length(lastintopfeat)) && (length(decorrelatedFetureList) == length(lastdecorrelated)) )
          {
            if ((sum(intopfeat %in% lastintopfeat) == length(intopfeat)) && (sum(decorrelatedFetureList %in% lastdecorrelated) == length(decorrelatedFetureList)))
            {
              addedlist <- 0;
            }
          }          
          lastdecorrelated <- decorrelatedFetureList;
          lastintopfeat <- intopfeat;
        }

      }
      betamatrix <- NULL;
      tmparincluded <- varincluded
      varincluded <- tmparincluded[tmparincluded %in% totused];
#      alphaincluded <- tmparincluded[tmparincluded %in% totalpha];
      alphaincluded <- tmparincluded[tmparincluded %in% totused];
      if (useDeCorr)
      {
        DeCorrmatrix <- DeCorrmatrix[alphaincluded,varincluded]
        dataAdjusted <- data
        if (length(varincluded) > 1)
        {
#          dataAdjusted[,varincluded] <- as.matrix(data[,varincluded]) %*% DeCorrmatrix;
          dataAdjusted[,varincluded] <- Rfast::mat.mult(as.matrix(data[,alphaincluded]),DeCorrmatrix);
          colsd <- apply(dataAdjusted[,varincluded],2,sd,na.rm = TRUE);
          if (sum(colsd==0) > 0)
          {
            zerovar <- varincluded[colsd==0];
            for (zcheck in zerovar)
            {
              dataAdjusted[,zcheck] <- dataAdjusted[,zcheck] + rnorm(nrow(dataAdjusted),0,1e-10);
            }
          }

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
          if (verbose) 
          {
             if (useFastCor)
             {
              cormat <- abs(Rfast::cora(as.matrix(dataAdjusted[,varincluded])))
             }
             else
             {
              cormat <- abs(cor(dataAdjusted[,varincluded],method=method))
             }
              diag(cormat) <- 0;
              cat (",",max(cormat),". Cor to Base:",length(correlatedToBase),", ABase:",length(AbaseFeatures),"\n")
          }
          newnames <- colnames(dataAdjusted)
          newnames[newnames %in% c(AbaseFeatures,baseFeatures)] <- paste("Ba_",newnames[newnames %in% c(AbaseFeatures,baseFeatures)],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          colnames(dataAdjusted) <- newnames
        }
      }
      else
      {
        dataAdjusted <- dataAdjusted[dataids,];
      }
    }
  }
  
  
  
#  cat ("\n")
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"DeCorrmatrix") <- DeCorrmatrix;
  attr(dataAdjusted,"varincluded") <- varincluded;
  attr(dataAdjusted,"baseFeatures") <- baseFeatures;
  attr(dataAdjusted,"useDeCorr") <- useDeCorr;
  attr(dataAdjusted,"correlatedToBase") <- correlatedToBase;
  attr(dataAdjusted,"AbaseFeatures") <- AbaseFeatures;
  return(dataAdjusted)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"DeCorrmatrix")))
  {
    decorMat <- attr(decorrelatedobject,"DeCorrmatrix")
#    testData[,colnames(decorMat)] <- as.matrix(testData[,rownames(decorMat)]) %*% decorMat
    testData[,colnames(decorMat)] <- Rfast::mat.mult(as.matrix(testData[,rownames(decorMat)]),decorMat);
    AbaseFeatures <- attr(decorrelatedobject,"AbaseFeatures")
    baseFeatures <- attr(decorrelatedobject,"baseFeatures")
    varincluded <- attr(decorrelatedobject,"varincluded")
    newnames <- colnames(testData)
    newnames[newnames %in% c(AbaseFeatures,baseFeatures)] <- paste("Ba_",newnames[newnames %in% c(AbaseFeatures,baseFeatures)],sep="") 
    newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
    colnames(testData) <- newnames
  }
  else
  {
    warning("The object does not have a decorrelation Matrix\n")
  }
  return (testData)
}

