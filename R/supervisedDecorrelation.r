featureDecorrelation <- function(data=NULL,
                                  thr=0.80,
                                  refdata=NULL,
                                  Outcome=NULL,
                                  baseFeatures=NULL,
                                  unipvalue=0.05,
                                  useDeCorr=TRUE,
                                  maxLoops=100,
                                  verbose=FALSE,
                                  method=c("fast","pearson","spearman","kendall"),
                                  skipRelaxed=TRUE,
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

getAllBetaCoefficients <- function(feat,varlist=NULL)
{
  allBetaCoef <- numeric(length(varlist));
  featVector <- refdata[,feat]
#  print(featVector[1:10])

  getBetaCoefficient <- function(x)
  {
    betaCoef <- 0;
    modellm <- try(lm(x~featVector,model = FALSE,na.action=na.exclude))
    if (!inherits(modellm, "try-error"))
    {	
      f <- summary(modellm)$coefficients
      betaCoef <- modellm$coef[2]*(f[2,4] < unipvalue);
    }
    return (betaCoef)
  }
#  cat(feat,"(",length(varlist),")","\n")
  if (length(varlist) > 1)
  {
    allBetaCoef <- apply(as.data.frame(refdata[,varlist]),2,getBetaCoefficient)
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
  fsocre <- numeric(length(varincluded));
  names(fsocre) <- varincluded;


  isFactor <- sapply(refdata[,varincluded],class) == "factor";
  isContinous <- unlist(lapply(lapply(refdata[,varincluded],unique),length)) > minUniqueValues; 

  varincluded <- varincluded[!isFactor & isContinous];
#  print(sum(!isFactor));
#  print(sum(isContinous));
#  print(sum(!isFactor & isContinous));
#  print(length(varincluded));
  
  
  
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
  varincluded <- names(maxcor)[maxcor >= 0.90*thr];
  DeCorrmatrix <- NULL;
  asociatedtoOutcome <- character();
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  if (length(varincluded) > 1)
  {
      if (!is.null(Outcome))
      {
          outcomep <- univariate_correlation(data[,c(Outcome,varincluded)],Outcome,method=method,limit=0,pvalue=2*unipvalue,thr = 0.999*thr) # the top associated features to the outcome
          if (length(baseFeatures)==0) baseFeatures <- names(outcomep);
          asociatedtoOutcome <- attr(outcomep,"Unadjusted")
          asociatedtoOutcome <- p.adjust(asociatedtoOutcome,"BH"); # adjusting for false discovery the unadjusted pvalues
          asociatedtoOutcome <- names(asociatedtoOutcome[asociatedtoOutcome < unipvalue])
          outcomep <- NULL;
    #      print(asociatedtoOutcome)
      }
      else
      {
         maxcor <- apply(cormat,2,max)
         AbaseFeatures <- varincluded;
         names(AbaseFeatures) <- AbaseFeatures;
         ordcor <- maxcor
         cormat[cormat < thr] <- 0;
         ordcor <- 0.95*maxcor + 
                      0.02*(apply(cormat,2,sum)/(0.001 + apply(1*(cormat >= thr),2,sum))) +
                      0.02*apply(cormat,2,mean)
         AbaseFeatures <- AbaseFeatures[order(-ordcor[AbaseFeatures])];
         AbaseFeatures <- correlated_Remove(cormat,AbaseFeatures,thr = thr,isDataCorMatrix=TRUE)
      }
      unipvalue <- min(unipvalue,2.0*unipvalue*sqrt(totcorr)/totFeatures); ## Adjusting for false association
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
      if (verbose) cat ("\n Included:",length(varincluded),", Uni p:",unipvalue,"To Outcome:",length(asociatedtoOutcome),", Base:",length(baseFeatures),", In Included:",length(baseIncluded),", Base Cor:",length(bvarincluded),"\n")
      
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;

      cormat <- cormat[varincluded,varincluded];
      bfeat <- unique(c(baseIncluded,AbaseFeatures))
      thr2 <- thr
      thr <- thr2*1.001;
      wthr <- 0.5 - 0.5*skipRelaxed;
      while (((addedlist > 0) || (thr > (thr2*1.0001))) && (lp < maxLoops)) 
      {
        lp = lp + 1;
        addedlist <- 0;

        maxcor <- apply(cormat,2,max)
        cormat[cormat < thr2] <- 0;
        
        thr <- max(c(thr2,thr2 + wthr*(max(maxcor) - thr2)))

        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        ordcor <- maxcor
        ordcor <- 0.95*maxcor + 
                  0.02*(apply(cormat,2,sum)/(0.001 + apply(1*(cormat >= thr2),2,sum))) +
                  0.02*apply(cormat,2,mean)
        if (length(bfeat) > 0)
        {
          ordcor[bfeat] <- ordcor[bfeat] + 0.01;
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
           intopfeat <- character();
          toBeDecorrelated <- length(topfeat)
          if (verbose)  cat(lp,", Top:",toBeDecorrelated,"<",thr,">");
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
#                if (verbose && (feat==topfeat[1]))  cat("<");
                falive <- varincluded[!(varincluded %in% c(decorrelatedFetureList,bfeat))]
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
                varlist <- varlist[!(varlist %in% decorrelatedFetureList)]
                ovarlist <- varlist;

                if (length(varlist) > 0)
                {
                    if (length(asociatedtoOutcome) > 0)
                    {
                      if (feat %in% asociatedtoOutcome)
                      {
                        varlist <- varlist[(varlist %in% asociatedtoOutcome)]
                      }
                      else
                      {
                        varlist <- varlist[!(varlist %in% asociatedtoOutcome)]
                      }
                      if (length(varlist) == 0)
                      {
#                        if (verbose && (feat==topfeat[1]))  cat("[",feat %in% asociatedtoOutcome,"]");
                        varlist <- names(corlist)
                        varlist <- varlist[!(varlist %in% decorrelatedFetureList)]
                      }
                    }
#                   if (verbose && (feat==topfeat[1]))  cat("(",length(varlist),")");
                   if (useFastCor)
                   {
                      prebetas <- getAllBetaCoefficients(feat,varlist);
                      varlist <- varlist[prebetas != 0];
                      if (length(varlist) > 0)
                      {
                        betamatrix[feat,varlist] <- -1.0*prebetas;
                        featAdded <- c(featAdded,feat);
                        intopfeat <- unique(c(intopfeat,feat));
                        if (verbose && (length(intopfeat) %% 100 == 99) && (lpct==1)) cat(".")
                        countf[varlist] <- countf[varlist] + 1;
                        decorrelatedFetureList <- c(decorrelatedFetureList,varlist);
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
                              dataAdjusted[,c(feat,varlist)] <- adataframe[,c(feat,varlist)];
                              refdata[,c(feat,varlist)] <- adataframe[refdataids,c(feat,varlist)];
                              intopfeat <- unique(c(intopfeat,feat));      
                              if (verbose && (length(intopfeat) %% 100 == 99)) cat(".")
                              countf[varlist] <- countf[varlist] + 1;
                              decorrelatedFetureList <- c(decorrelatedFetureList,varlist);
                          }
                      }
                  }
                  fsocre[feat] <- fsocre[feat] + sum(cormat[varlist,feat]^2)
                  fsocre[varlist] <- fsocre[varlist] - cormat[varlist,feat]^2
                  cormat[ovarlist,feat] <- 0.99*thr;

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
          if (verbose) cat("[",lpct,":",length(featMarked),":",maxnotf,"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
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
                DeCorrmatrix[totalphaused,colused] <-  as.numeric(DeCorrmatrix[totalphaused,allused])*
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
                   cat("..no Fast..[",ncol(mtx1),"][",nrow(mtx2),"]")
                   DeCorrmatrix[totalphaused,colused] <- mtx1 %*% mtx2;
                }
              }
             if (useFastCor)
             {
                if (verbose) cat("|")
#                refdata[,varincluded] <- as.matrix(dataAdjusted[refdataids,varincluded]) %*% DeCorrmatrix
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
                     cat("..no Fast..[",ncol(mtx1),"][",nrow(mtx2),"]")
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
            cat ("Tot Used:",length(totused),", Added:",addedlist,", Zero Std:",sum(colsd==0),", Max Cor:",max(cormat),"\n");
          }
          if ( (length(intopfeat) == length(lastintopfeat)) && (length(decorrelatedFetureList) == length(lastdecorrelated)) )
          {
            if ((sum(intopfeat %in% lastintopfeat) == length(intopfeat)) && (sum(decorrelatedFetureList %in% lastdecorrelated) == length(decorrelatedFetureList)))
            {
              if (verbose) cat (decorrelatedFetureList);
              addedlist <- 0;
            }
          }          
          lastdecorrelated <- decorrelatedFetureList;
          lastintopfeat <- intopfeat;
        }

        if ((wthr > 0) && ((addedlist <= max(c(1.0,0.25*length(topfeat)))) || (length(intopfeat) <= 1)))
        {
          wthr <- wthr - 0.25;
        }
      }
      betamatrix <- NULL;
      tmparincluded <- varincluded
      varincluded <- tmparincluded[tmparincluded %in% totused];
      if (useDeCorr)
      {
        DeCorrmatrix <- DeCorrmatrix[varincluded,varincluded]
        dataAdjusted <- data
        if (length(varincluded) > 1)
        {
#          dataAdjusted[,varincluded] <- as.matrix(data[,varincluded]) %*% DeCorrmatrix;
          dataAdjusted[,varincluded] <- Rfast::mat.mult(as.matrix(data[,varincluded]),DeCorrmatrix);
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
              cat ("[",lp,"],",max(cormat),". Cor to Base:",length(correlatedToBase),", ABase:",length(AbaseFeatures),"\n")
          }
          newnames <- colnames(dataAdjusted)
          newnames[newnames %in% c(AbaseFeatures,baseFeatures)] <- paste("Ba_",newnames[newnames %in% c(AbaseFeatures,baseFeatures)],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          colnames(dataAdjusted) <- newnames
          newnames <- colnames(DeCorrmatrix)
          newnames[newnames %in% c(AbaseFeatures,baseFeatures)] <- paste("Ba_",newnames[newnames %in% c(AbaseFeatures,baseFeatures)],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          colnames(DeCorrmatrix) <- newnames
          fsocre[!(names(fsocre) %in% varincluded)] <- min(fsocre) - 1.0;
          newnames <- names(fsocre)
          newnames[newnames %in% c(AbaseFeatures,baseFeatures)] <- paste("Ba_",newnames[newnames %in% c(AbaseFeatures,baseFeatures)],sep="") 
          newnames[newnames %in% varincluded] <- paste("De_",newnames[newnames %in% varincluded],sep="")
          names(fsocre) <- newnames
        }
      }
      else
      {
        dataAdjusted <- dataAdjusted[dataids,];
      }
  }
  
  
  attr(dataAdjusted,"topFeatures") <- unique(topFeatures);
  attr(dataAdjusted,"TotalAdjustments") <- countf;
  attr(dataAdjusted,"DeCorrmatrix") <- DeCorrmatrix;
  attr(dataAdjusted,"varincluded") <- varincluded;
  attr(dataAdjusted,"baseFeatures") <- baseFeatures;
  attr(dataAdjusted,"useDeCorr") <- useDeCorr;
  attr(dataAdjusted,"correlatedToBase") <- correlatedToBase;
  attr(dataAdjusted,"AbaseFeatures") <- AbaseFeatures;
  attr(dataAdjusted,"fsocre") <- fsocre;
  return(dataAdjusted)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"DeCorrmatrix")))
  {
    decorMat <- attr(decorrelatedobject,"DeCorrmatrix")
#    testData[,colnames(decorMat)] <- as.matrix(testData[,rownames(decorMat)]) %*% decorMat
    testData[,rownames(decorMat)] <- Rfast::mat.mult(as.matrix(testData[,rownames(decorMat)]),decorMat);
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

