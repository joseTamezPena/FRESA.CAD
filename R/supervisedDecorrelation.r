GDSTMDecorrelation <- function(data=NULL,
                                  thr=0.80,
                                  refdata=NULL,
                                  Outcome=NULL,
                                  baseFeatures=NULL,
                                  unipvalue=0.05,
                                  useDeCorr=TRUE,
                                  maxLoops=100,
                                  verbose=FALSE,
                                  method=c("fast","pearson","spearman","kendall"),
                                  skipRelaxed=FALSE,
                                  corRank=TRUE,
#                                  pcl=NULL,
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
  minUniqueValues <- 5 # more than 5 values and is not factor
  

featVector <<- data[,1];

getAllBetaCoefficients <- function(feat,varlist=NULL)
{
  allBetaCoef <- numeric(length(varlist));
  featVector <<- refdata[,feat]
#  print(featVector[1:10])


  getBetaCoefficient <- function(x)
  {
    betaCoef <- 0;
    modellm <- try(lm(x~featVector,model = FALSE,na.action=na.exclude))
    if (!inherits(modellm, "try-error"))
    {	
      f <- summary(modellm)$coefficients
      if (nrow(f)>1)
      {
        betaCoef <- modellm$coef[2]*(f[2,4] < unipvalue);
      }
    }
    return (betaCoef)
  }
#  cat(feat,"(",length(varlist),")","\n")
  if (length(varlist) > 1)
  {
#    if ((length(varlist) > 31) && !is.null(pcl))
#    {
#      clusterExport(cl=pcl, 'featVector')
#      allBetaCoef <- parApply(pcl,refdata[,varlist],2,getBetaCoefficient)
#    }
#    else
#    {
      allBetaCoef <- apply(refdata[,varlist],2,getBetaCoefficient)
#    }
  }
  else
  {
    allBetaCoef <- getBetaCoefficient(refdata[,varlist])
  }
  allBetaCoef[is.na(allBetaCoef)] <- 0;
  allBetaCoef[is.nan(allBetaCoef)] <- 0;
  allBetaCoef[is.infinite(allBetaCoef)] <- 0;
  names(allBetaCoef) <- varlist
  return (allBetaCoef)
}

  if (thr < 0.01) # The smallest correlation
  {
    thr <- 0.01 
  }
  if (thr >= 1.0) #largest correlation
  {
    thr <- 0.999
  }

  dataids <- rownames(data)

  if (is.null(refdata))
  {
    refdata <- data;
  }
  refdataids <- rownames(refdata);
  dataTransformed <- as.data.frame(rbind(data,refdata[!(refdataids %in% dataids),]));
  topFeatures <- character()
  correlatedToBase <- character();
  if (is.null(baseFeatures))
  {
    baseFeatures <- character()
  }
  AbaseFeatures <- character()
  bfeat <- unique(c(baseFeatures,AbaseFeatures))
  
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


  isFactor <- sapply(refdata[,varincluded],class) == "factor";
  isContinous <- unlist(lapply(lapply(refdata[,varincluded],unique),length)) > minUniqueValues; 

  varincluded <- varincluded[!isFactor & isContinous];
  
  
  fscore <- numeric(length(varincluded));
  names(fscore) <- varincluded;
  
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
  totcorr <- sum(cormat >= max(thr,0.5)); 
  varincluded <- names(maxcor)[maxcor >= 0.95*thr];
  DeCorrmatrix <- NULL;
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  if (length(varincluded) > 1)
  {
      if (!is.null(Outcome))
      { 
          if (length(baseFeatures)==0) 
          {
            if ((class(data[,Outcome]) == "factor") | (length(unique(data[,Outcome])) == 2))
            {
              outcomep <- univariate_correlation(data,Outcome,limit=0,pvalue=unipvalue,thr = 0.95*thr) # the top associated features to the outcome
            }
            else
            {
              outcomep <- univariate_correlation(data,Outcome,method=method,limit=0,pvalue=unipvalue,thr = 0.95*thr) # the top associated features to the outcome
            }
            baseFeatures <- names(outcomep);
            if (length(baseFeatures) > 1)
            {
              baseFeatures <- as.character(correlated_Remove(cormat,baseFeatures,thr = 0.95*thr,isDataCorMatrix=TRUE))
            }
          }
          outcomep <- NULL;
      }
       maxcor <- apply(cormat,2,max)
       AbaseFeatures <- varincluded;
       names(AbaseFeatures) <- AbaseFeatures;
       ordcor <- maxcor;
       cormat[cormat < thr] <- 0;
       if (corRank)
       {       
          ordcor <- maxcor + thr*apply(cormat-thr,2,sum)/(1.0-thr);
       }
      AbaseFeatures <- AbaseFeatures[order(-ordcor[AbaseFeatures])];
      AbaseFeatures <- as.character(correlated_Remove(cormat,AbaseFeatures,thr = 0.95*thr,isDataCorMatrix=TRUE))
      AbaseFeatures <- AbaseFeatures[order(-maxcor[AbaseFeatures])];
      unipvalue <- min(unipvalue,unipvalue*sqrt(totcorr)/totFeatures); ## Adjusting for false association
      baseIncluded <- character();
      if (length(baseFeatures) > 0)
      {    
        baseIncluded <- varincluded[varincluded %in% baseFeatures];
        if (length(baseIncluded)>1)
        {
          names(baseIncluded) <- baseIncluded;
        }
      }
      
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;

      cormat <- cormat[varincluded,varincluded];
      baseFeatures <- baseIncluded;
      bfeat <- unique(c(baseIncluded,AbaseFeatures))
      bfeat <- as.character(correlated_Remove(cormat,bfeat,thr = 0.95*thr,isDataCorMatrix=TRUE));
      if (verbose) cat ("\n Included:",length(varincluded),", Uni p:",unipvalue,", Uncorrelated Base:",length(AbaseFeatures),", Outcome-Driven Size:",length(baseIncluded),", Base Size:",length(bfeat),"\n")

      thr2 <- thr
      thr <- thr2*1.001;
      wthr <- 0.95 - 0.95*skipRelaxed;
      while (((addedlist > 0) || (thr > (thr2*1.0001))) && (lp < maxLoops)) 
      {
        lp = lp + 1;

        addedlist <- 0;

        
        cormat[cormat < thr2] <- 0;
        maxcor <- apply(cormat,2,max)

        if (verbose)  cat("\n\r",lp,"Added:",addedlist,sprintf("<R=%5.3f,w=%5.3f,thr=%5.3f>",max(maxcor),wthr,thr))

        if (max(maxcor) >= thr2)
        {
          thr <- max(c(thr2,thr2 + wthr*(max(maxcor) - thr2)));
        }
        else
        {
          thr <- thr2;
        }

        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        ordcor <- maxcor;
        if (corRank)
        {
          ordcor <- maxcor + thr2*apply(cormat-thr2,2,sum)/(1.0-thr2);
        }

        if (length(bfeat) > 0)
        {
          ordcor[bfeat] <- ordcor[bfeat] + 10.0*ncol(cormat);
        }
        
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (length(topfeat)>0)
        {
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          topfeat <- topfeat[order(-ordcor[topfeat])];
          topfeat <- correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE)
          topfeat <- topfeat[order(-maxcor[topfeat])];
          intopfeat <- character();
          toBeDecorrelated <- length(topfeat)
          if (verbose)  cat(", Top:",toBeDecorrelated,"<",sprintf("%5.3f",thr),">");
          maxcortp <- thr;
          maxnotf <- thr;
          tobereviewed <- topfeat;
          featAdded <- character(1);
          featMarked <- character();
          wcor <- 0.95 - 0.95*skipRelaxed;
          lpct <- 0;
          while (length(featAdded) > 0)
          {
            lpct <- lpct + 1;
            featAdded <- character();
            topfeat <- topfeat[!(topfeat %in% featMarked)]
            for (feat in topfeat)
            {
                falive <- varincluded[!(varincluded %in% decorrelatedFetureList)]
                tobereviewed <- topfeat[!(topfeat %in% c(featMarked,feat))]
                if ((length(tobereviewed) > 1) && (length(falive) > 1))
                {
                  maxnotf <- max(apply(cormat[falive,tobereviewed],2,max))
                  maxcortp <- max(c(wcor*maxnotf,thr))
                }
                else
                {
                  maxcortp <- thr;
                }
                corlist <- cormat[,feat];
                corlist <- corlist[corlist >= maxcortp];

                varlist <- names(corlist)
                varlist <- varlist[!(varlist %in% unique(c(decorrelatedFetureList,bfeat)))]
                ovarlist <- varlist;

                if (length(varlist) > 0)
                {
                   if (verbose && (feat==topfeat[1]))  cat("(",length(varlist),")");
                   if (useFastCor)
                   {
                      prebetas <- getAllBetaCoefficients(feat,varlist);
                      varlist <- varlist[prebetas != 0];
                      if (length(varlist) > 0)
                      {
                        betamatrix[feat,varlist] <- -1.0*prebetas[prebetas != 0];
                        featAdded <- c(featAdded,feat);
                        intopfeat <- unique(c(intopfeat,feat));
                        if (verbose && (length(intopfeat) %% 100 == 99) && (lpct==1)) cat(".")
                        countf[varlist] <- countf[varlist] + 1;
                        decorrelatedFetureList <- c(decorrelatedFetureList,varlist);
                      }
                   }
                   else
                   {
                     adataframe <- featureAdjustment(varlist,
                                                    baseFormula=feat,
                                                    data=dataTransformed[,c(feat,varlist)],
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
                          varlist <- varlist[adjusted];
                          if (length(varlist) > 0)
                          {
                              featAdded <- c(featAdded,feat);
                              dataTransformed[,c(feat,varlist)] <- adataframe[,c(feat,varlist)];
                              refdata[,c(feat,varlist)] <- adataframe[refdataids,c(feat,varlist)];
                              intopfeat <- unique(c(intopfeat,feat));      
                              if (verbose && (length(intopfeat) %% 100 == 99)) cat(".")
                              countf[varlist] <- countf[varlist] + 1;
                              decorrelatedFetureList <- c(decorrelatedFetureList,varlist);
                          }
                      }
                  }
                  fscore[feat] <- fscore[feat] + sum(cormat[varlist,feat]^2)
                  fscore[varlist] <- fscore[varlist] - cormat[varlist,feat]^2
                  cormat[ovarlist,feat] <- 0;

                }
#                if (verbose && (feat==topfeat[1]))  cat(">");
                if ((length(varlist) == 0) && (maxnotf <= thr2))
                {
                    featMarked <- unique(c(featMarked,feat));
                }
            }
#            if (verbose) cat("|",maxnotf,"|<",length(featAdded),">{",length(featMarked),"}");
            if (skipRelaxed || ((toBeDecorrelated == (length(unique(c(featMarked,featAdded))))) && (maxnotf <= thr)))
            {
              featAdded <- character();
            }
            else
            {
              if ((length(featAdded) == 0) && (maxnotf > thr))
              {
                 wcor <- wcor - 0.2;
                 featAdded <- character(1);
              }
            }
          }
          addedlist <- length(decorrelatedFetureList);
          if (verbose) cat("[",lpct,":",length(featMarked),":",sprintf("%5.3f",maxnotf),"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
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
              if ((length(colused) == 1) && (length(allused) == 1))
              {
                DeCorrmatrix[totalphaused,colused] <-  as.numeric(DeCorrmatrix[totalphaused,allused]) *
                                                               as.numeric(betamatrix[allused,colused]);
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
                   mtx1 <- t(as.matrix(DeCorrmatrix[totalphaused,allused]))
                   mtx2 <- as.matrix(betamatrix[allused,colused])
#                   cat("..no Fast..[",ncol(mtx1),"][",nrow(mtx2),"]")
                   DeCorrmatrix[totalphaused,colused] <- mtx1 %*% mtx2;
                }
              }
             if (useFastCor)
             {
                if (verbose) cat("|")
#                refdata[,varincluded] <- as.matrix(dataTransformed[refdataids,varincluded]) %*% DeCorrmatrix
                if (( length(colused) == 1) && (length(alphaused) == 1 ))
                {
                  refdata[,colused] <- refdata[,colused] + as.numeric(refdata[,alphaused])*as.numeric(betamatrix[alphaused,colused])
                }
                else
                {
                  if (length(alphaused) > 1)
                  {
                    refdata[,colused] <- refdata[,colused] + Rfast::mat.mult(as.matrix(refdata[,alphaused]),as.matrix(betamatrix[alphaused,colused]))
                  }
                  else
                  {
                     mtx1 <- as.matrix(refdata[,alphaused])
                     mtx2 <- t(as.matrix(betamatrix[alphaused,colused]))
#                     cat("..no Fast..[",ncol(mtx1),"][",nrow(mtx2),"]")
                     refdata[,colused] <- refdata[,colused] + mtx1 %*% mtx2;
                  }
                }
             }
          }
          if (verbose) cat(">")
          colsd <- apply(refdata[,varincluded],2,sd,na.rm = TRUE);
          if (sum(colsd==0) > 0)
          {
            zerovar <- varincluded[colsd==0];
            for (zcheck in zerovar)
            {
              refdata[,zcheck] <- refdata[,zcheck] + rnorm(nrow(refdata),0,1.0e-10);
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
            cat ("Tot Used:",length(totused),", Added:",addedlist,", Zero Std:",sum(colsd==0),", Max Cor:",sprintf("%5.3f",max(cormat)));
          }
          if ( (length(intopfeat) == length(lastintopfeat)) && (length(decorrelatedFetureList) == length(lastdecorrelated)) )
          {
            if ((sum(intopfeat %in% lastintopfeat) == length(intopfeat)) && (sum(decorrelatedFetureList %in% lastdecorrelated) == length(decorrelatedFetureList)))
            {
#              if (verbose) cat (decorrelatedFetureList);
              addedlist <- 0;
            }
          }          
          lastdecorrelated <- decorrelatedFetureList;
          lastintopfeat <- intopfeat;
        }

        if ((wthr > 0) && ((addedlist <= max(c(1.0,0.25*length(topfeat)))) || (length(intopfeat) <= 1)))
        {
          wthr <- wthr - 0.20;
        }
      }
      betamatrix <- NULL;
      tmparincluded <- varincluded
      varincluded <- tmparincluded[tmparincluded %in% totused];
      if (useDeCorr)
      {
        DeCorrmatrix <- DeCorrmatrix[varincluded,varincluded]
        dataTransformed <- data
        if (length(varincluded) > 1)
        {
#          dataTransformed[,varincluded] <- as.matrix(data[,varincluded]) %*% DeCorrmatrix;
          dataTransformed[,varincluded] <- Rfast::mat.mult(as.matrix(data[,varincluded]),DeCorrmatrix);
          colsd <- apply(dataTransformed[,varincluded],2,sd,na.rm = TRUE);
          if (sum(colsd==0) > 0)
          {
            zerovar <- varincluded[colsd==0];
            for (zcheck in zerovar)
            {
              dataTransformed[,zcheck] <- dataTransformed[,zcheck] + rnorm(nrow(dataTransformed),0,1e-10);
            }
          }

          correlatedToBase <- varincluded[varincluded %in% bfeat];
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
          correlatedToBase <- correlatedToBase[!(correlatedToBase %in% bfeat)];
          if (verbose) 
          {
             if (useFastCor)
             {
              cormat <- abs(Rfast::cora(as.matrix(dataTransformed[,varincluded])))
             }
             else
             {
              cormat <- abs(cor(dataTransformed[,varincluded],method=method))
             }
              diag(cormat) <- 0;
              cat ("\n\r [",lp,"],",max(cormat),". Cor to Base:",length(correlatedToBase),", ABase:",length(AbaseFeatures),"\n\r")
          }
#          banames <- colnames(DeCorrmatrix)[apply(DeCorrmatrix!=0,2,sum)==1]
          newnames <- colnames(dataTransformed)
          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          colnames(dataTransformed) <- newnames
          newnames <- colnames(DeCorrmatrix)
          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          colnames(DeCorrmatrix) <- newnames
          newnames <- names(fscore)
          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          names(fscore) <- newnames
        }
      }
      else
      {
        dataTransformed <- dataTransformed[dataids,];
      }
  }
  
  
  attr(dataTransformed,"topFeatures") <- unique(topFeatures);
  attr(dataTransformed,"TotalAdjustments") <- countf;
  attr(dataTransformed,"GDSTM") <- DeCorrmatrix;
  attr(dataTransformed,"varincluded") <- varincluded;
  attr(dataTransformed,"baseFeatures") <- baseFeatures;
  attr(dataTransformed,"uniqueBase") <- bfeat;
  attr(dataTransformed,"useDeCorr") <- useDeCorr;
  attr(dataTransformed,"correlatedToBase") <- correlatedToBase;
  attr(dataTransformed,"AbaseFeatures") <- AbaseFeatures;
  attr(dataTransformed,"fscore") <- fscore;
  attr(dataTransformed,"unipvalue") <- unipvalue
  return(dataTransformed)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"GDSTM")))
  {
    decorMat <- attr(decorrelatedobject,"GDSTM")
#    testData[,colnames(decorMat)] <- as.matrix(testData[,rownames(decorMat)]) %*% decorMat
    testData[,rownames(decorMat)] <- Rfast::mat.mult(as.matrix(testData[,rownames(decorMat)]),decorMat);
    varincluded <- attr(decorrelatedobject,"varincluded")
    newnames <- colnames(testData)
    bfeat <- attr(decorrelatedobject,"uniqueBase")
    newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
    newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
    colnames(testData) <- newnames
  }
  else
  {
    warning("The object does not have a decorrelation Matrix\n")
  }
  return (testData)
}

getDerivedCoefficients <- function(decorrelatedobject)
{
  CoeffList <- list();
  nonZeronames <- character();
  GDSTM <- attr(decorrelatedobject,"GDSTM")
  namecol <- colnames(GDSTM);

  n=1;
  for (i in 1:ncol(GDSTM))
  {
    associ <- abs(GDSTM[,i])>0;
    if (sum(associ) > 1)
    {
      nonZeronames <- append(nonZeronames,namecol[i]);
      CoeffList[[n]] <- GDSTM[associ,i];
      n=n+1;
    }
  }
  names(CoeffList) <- nonZeronames;
  
  return (CoeffList)
}
