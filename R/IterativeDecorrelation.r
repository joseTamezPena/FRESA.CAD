ILAA <- function(data=NULL,
                thr=0.80,
                method=c("pearson","spearman"),
                Outcome=NULL,
                drivingFeatures=NULL,
                maxLoops=100,
                verbose=FALSE,
                ...)
{
  method <- match.arg(method);
  if (method=="pearson") 
  {
    method <- "fast";
  }
  if (verbose)
  {
    cat(method,"|")
  }
  
  result <- IDeA(data=data,
                 thr=thr,
                 method=method,
                 Outcome=Outcome,
                 refdata=NULL,
                 drivingFeatures=drivingFeatures,
                 useDeCorr=TRUE,
                 relaxed=TRUE,
                 corRank=TRUE,
                 maxLoops=maxLoops,
                 unipvalue=0.05,
                 verbose=verbose,
                 ...)
   return(result)

}                                  

IDeA <- function(data=NULL,
                                  thr=0.80,
                                  method=c("fast","pearson","spearman","kendall"),
                                  Outcome=NULL,
                                  refdata=NULL,
                                  drivingFeatures=NULL,
                                  useDeCorr=TRUE,
                                  relaxed=TRUE,
                                  corRank=TRUE,
                                  maxLoops=100,
                                  unipvalue=0.05,
                                  verbose=FALSE,
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
  correlationMeasureEvolution <- numeric()
  sparcityEvolution <- numeric();
  

  featVector <- data[,1];

  getAllBetaCoefficients <- function(feat,varlist=NULL)
  {
    allBetaCoef <- numeric(length(varlist));
    featVector <- refdata[,feat]


    getBetaCoefficient <- function(x)
    {
      betaCoef <- 0;
      modellm <- try(lm(x~featVector,model = FALSE,na.action=na.exclude))
      if (!inherits(modellm, "try-error"))
      {	
        f <- summary(modellm)$coefficients
        if (nrow(f)>1)
        {
          betaCoef <- modellm$coef[2]*(f[2,4] < adjunipvalue);
        }
      }
      return (betaCoef)
    }
    if (length(varlist) > 1)
    {
        allBetaCoef <- apply(refdata[,varlist],2,getBetaCoefficient)
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
  lavariables <- character();
  if (is.null(drivingFeatures))
  {
    drivingFeatures <- character()
  }
  AdrivingFeatures <- character()
  
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
  isContinous <- unlist(lapply(lapply(refdata[,varincluded],unique),length)) >= minUniqueValues; 

  varincluded <- varincluded[!isFactor & isContinous];
  
  
  fscore <- numeric(length(varincluded));
  names(fscore) <- varincluded;
#  fcount <- fscore;
  
  totFeatures <- length(varincluded);
  largeSet <- (totFeatures > 1000) || ((nrow(refdata) > 1000) && (totFeatures > 100))
  

  
  models <- NULL
  
  if (useFastCor)
  {
#    cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded]),large=largeSet))
    cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded])))
  }else
  {
    cormat <- abs(cor(refdata[,varincluded],method=method))
  }
  diag(cormat) <- 0;
  maxcor <- apply(cormat,2,max)
  totcorr <- sum(cormat >= 0.5); ## For False Discovery estimation
  totAtcorr <- sum(cormat >= thr); ## 
  DeCorrmatrix <- NULL;
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  outcomefeatures <- drivingFeatures;
  bfeat <- NULL
#  print(totFeatures)
  adjunipvalue <- min(c(unipvalue,unipvalue*(sqrt(1.0 + totcorr))/totFeatures)); ## Adjusting for false association
  
  ## Min Correlation based on the pearson distribution
  ndf <- nrow(data)-2
#  print(c(ndf,adjunipvalue))
  tvalue <- qt(1.0 - adjunipvalue,ndf)
#  print(tvalue)
  rcrit <- tvalue/sqrt(ndf + tvalue^2)
  if (thr < rcrit) thr <- rcrit
  
  ################################# end #######################
  althr <- 0.10;
  cortoInclude <- max(c(althr*thr,rcrit))
  
  varincluded <- names(maxcor)[maxcor >= cortoInclude];
  if (length(varincluded) > 2000) 
  { 
    althr <- 0.5;
    cortoInclude <- max(c(althr*thr,rcrit))
    varincluded <- names(maxcor)[maxcor >= cortoInclude];
  }    
  allFeatAdded <- character()
  decordim <- 0
  if (length(varincluded) > 1)
  {
      if (!is.null(Outcome))
      { 
          if (length(drivingFeatures)==0) 
          {
            if (inherits(data[,Outcome],"factor") | (length(unique(data[,Outcome])) == 2))
            {
              outcomep <- univariate_correlation(data,Outcome,limit=0,pvalue= 0.1*unipvalue,thr = 0.95*thr) # the top associated features to the outcome
            }
            else
            {
              outcomep <- univariate_correlation(data,Outcome,method=method,limit=0,pvalue= 0.1*unipvalue,thr = 0.95*thr) # the top associated features to the outcome
            }
            outcomefeatures <- attr(outcomep,"Unadjusted");
            outcomefeatures <- names(outcomefeatures[outcomefeatures < unipvalue])
            drivingFeatures <- names(outcomep);
            outcomep <- NULL;
          }
      }
      cormat <- cormat[varincluded,varincluded];
      AdrivingFeatures <- varincluded;
      maxcor <- apply(cormat,2,max)
      if (corRank)
      {       
        axcor <- cormat;
        axcor[axcor < thr] <- 0;
        ordcor <- apply(axcor^2,2,mean);
      }
      else
      {
        ordcor <- apply(cormat,2,mean)
        names(AdrivingFeatures) <- AdrivingFeatures;
        AdrivingFeatures <- AdrivingFeatures[order(-ordcor[AdrivingFeatures])];
        ordcor <- apply(1*(cormat>=thr),2,mean)
        AdrivingFeatures <- AdrivingFeatures[order(-ordcor[AdrivingFeatures])];
        ordcor <- maxcor;
        AdrivingFeatures <- AdrivingFeatures[order(-ordcor[AdrivingFeatures])];
      }
#      print(AdrivingFeatures)
      AdrivingFeatures <- as.character(correlated_Remove(cormat,AdrivingFeatures,thr = 0.75*thr,isDataCorMatrix=TRUE))
      
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;
      decordim <- length(varincluded);
      bfeat <- AdrivingFeatures;
      if (length(drivingFeatures) > 0)
      {    
        drivingFeatures <- drivingFeatures[drivingFeatures %in% varincluded];
        if (length(drivingFeatures)>1)
        {
          drivingFeatures <- drivingFeatures[order(-ordcor[drivingFeatures])];
          names(drivingFeatures) <- drivingFeatures;
          bfeat <- unique(c(drivingFeatures,AdrivingFeatures))
          names(bfeat) <- bfeat;
          bfeat <- as.character(correlated_Remove(cormat,bfeat,thr = 0.75*thr,isDataCorMatrix=TRUE));
          bfeat <- bfeat[order(-ordcor[bfeat])];
          drivingFeatures <- drivingFeatures[drivingFeatures %in% bfeat];
          AdrivingFeatures <- AdrivingFeatures[AdrivingFeatures %in% bfeat];
        }
      }

      if (verbose) cat ("\n Included:",length(varincluded),", Uni p:",adjunipvalue,", Uncorrelated Base:",length(AdrivingFeatures),", Outcome-Driven Size:",length(drivingFeatures),", Base Size:",length(bfeat),"\n")

      thr2 <- thr
      thr <- thr2*1.001;
      thrvalues <- c(0.90,0.75,0.60,0.45,0.20,0.10);
      nextthr <- 1 + (!relaxed)*length(thrvalues);
      nextthra <- 0;
      while (((addedlist > 0) || (thr > (thr2*1.0001))) && (lp < maxLoops)) 
      {
        lp = lp + 1;

        addedlist <- 0;
        
        maxcor <- apply(cormat,2,max)
        mxScor <- max(maxcor);
        correlationMeasureEvolution <- c(correlationMeasureEvolution,mxScor);
        sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

        thr <- thr2;
        nextthr <- sum(thrvalues > mxScor) + 1
        if ((mxScor >= thr2) && (nextthr <= length(thrvalues)))
        {
          thr <- max(thr2,thrvalues[nextthr]);
        }

        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        if (nextthra != nextthr)
        {
          orglentopfeat <- sum(maxcor[topfeat] >= thr);
          nextthra <- nextthr;
        }
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (verbose)  cat("\n\r",lp,sprintf("<R=%5.3f,thr=%5.3f,N=%5d>",mxScor,thr,orglentopfeat))

        if (length(topfeat)>0)
        {
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          if (corRank)
          {
            axcor <- cormat;
            axcor[axcor < thr] <- 0;
            ordcor <- apply(axcor^2,2,mean);
          }
          else
          {
            ordcor <- apply(cormat,2,mean)
            topfeat <- topfeat[order(-ordcor[topfeat])];
            ordcor <- apply(1*(cormat >= thr2),2,mean)
            topfeat <- topfeat[order(-ordcor[topfeat])];
            ordcor <- maxcor;
          }
          ordcor <- ordcor + 1.0*(names(ordcor) %in% allFeatAdded)
          topfeat <- topfeat[order(-ordcor[topfeat])];
          atopbase <- bfeat[bfeat %in% topfeat];
          topfeat <- unique(c(atopbase,topfeat));
          topfeat <- correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE)
          topfeat <- topfeat[order(-maxcor[topfeat])];
          intopfeat <- character();
          toBeDecorrelated <- length(topfeat)
          if (verbose)  cat(", Top:",toBeDecorrelated);
          featAdded <- character(1);
          featMarked <- character();
          lpct <- 0;
          throff <- 0.01;
          shortCor <- NULL;
          themaxthr <- thr;
          while (length(featAdded) > 0)
          {
            lpct <- lpct + 1;
            featAdded <- character();
            topfeat <- topfeat[!(topfeat %in% featMarked)]
            if (length(topfeat) > 0)
            {
              if (length(topfeat) > 1) shortCor <- cormat[topfeat,];
              testedMedian <- FALSE;
              maxthr <- thr;
              for (feat in topfeat)
              {
                  corlist <- cormat[,feat];
                  corlist <- corlist[corlist >= thr];
                  varlist <- names(corlist)
                  varlist <- varlist[!(varlist %in% unique(c(decorrelatedFetureList,bfeat)))]
                  if (length(outcomefeatures) > 1)
                  {
                     ovarlist <- varlist;
                     if (feat %in% outcomefeatures) 
                     {
                        varlist <- varlist[varlist %in% outcomefeatures];
                     }
                     else
                     {
                        varlist <- varlist[!(varlist %in% outcomefeatures)];
                     }
                     if (length(varlist) == 0)
                     {
                        varlist <- ovarlist;
                     }
                  }
                  olength <- length(varlist)
                  notinfeat <- topfeat[topfeat != feat];
                  allvartested <- TRUE;
                  if ((olength > 1) && (length(notinfeat) > 1))
                  {
                      medianthr <- median(apply(shortCor[notinfeat,varlist],2,max));
                      if (medianthr >= (thr + throff))
                      {
                        if (verbose && (feat==topfeat[1]))  cat("[",length(varlist),"]");
                        if (maxthr < medianthr) maxthr <- medianthr;
                        themaxthr <- max(c(themaxthr,maxthr));
                        corlist <- corlist[corlist >= medianthr];
                        varlist <- names(corlist)
                        allvartested <- (length(varlist) == olength);
                        testedMedian <- testedMedian | !allvartested;
                      }
                  }

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
                                                      pvalue = adjunipvalue,
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
                              adjusted[models[[vl]]$feature] <- (models[[vl]]$pval < adjunipvalue);
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
                    fscore[feat] <- fscore[feat] + length(varlist);
                    fscore[varlist] <- fscore[varlist] - 1;
                    cormat[ovarlist,feat] <- 0;

                  }
                  if (allvartested)
                  {
                      featMarked <- unique(c(featMarked,feat));
                  }
              }
              allFeatAdded <- unique(c(allFeatAdded,featAdded))
              if (!testedMedian)
              {
                featAdded <- character();
              }
              if ((maxthr > thr) && testedMedian)
              {
                featAdded <- character(1);
                if (verbose) cat("=")
                throff <- throff + 0.025;
              }
            }
          }
          addedlist <- length(decorrelatedFetureList);
          if (verbose) cat("[",lpct,":",length(featMarked),"Fa=",length(allFeatAdded),":",sprintf("%5.3f",themaxthr),"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
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
              addedlist <- 0;
            }
          }          
          lastdecorrelated <- decorrelatedFetureList;
          lastintopfeat <- intopfeat;
        }

        if ((nextthr <= length(thrvalues)) && ((addedlist <= max(c(2.0,0.5*orglentopfeat))) || (length(intopfeat) <= 1)))
        {
          nextthr <- nextthr + 1; 
        }
      }
      correlationMeasureEvolution <- c(correlationMeasureEvolution,max(cormat));
      sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

      betamatrix <- NULL;
      tmparincluded <- varincluded;
      varincluded <- tmparincluded[tmparincluded %in% totused];
      drivingFeatures <- drivingFeatures[drivingFeatures %in% varincluded];
      AdrivingFeatures <- AdrivingFeatures[AdrivingFeatures %in% varincluded];
      bfeat <- bfeat[bfeat %in% varincluded];
      if (useDeCorr)
      {
        DeCorrmatrix <- DeCorrmatrix[varincluded,varincluded]
        abfeat <- varincluded[apply(DeCorrmatrix!=0,2,sum)==1]
        bfeat <- unique(c(bfeat,abfeat));
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
              cat ("\n\r [",lp,"],",max(cormat),"Decor Dimension:",length(varincluded),"Nused:",length(totused),". Cor to Base:",length(correlatedToBase),", ABase:",length(AdrivingFeatures),", Outcome Base:",length(drivingFeatures),"\n\r")
          }
          changed <- (apply(DeCorrmatrix!=0,2,sum) > 1);
          lavariables <- colnames(DeCorrmatrix)[changed]
          newnames <- colnames(dataTransformed)
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(dataTransformed) <- newnames
          newnames <- colnames(DeCorrmatrix)
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(DeCorrmatrix) <- newnames
          fscore <- fscore[varincluded];

          newnames <- names(fscore)
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          names(fscore) <- newnames
        }
      }
      else
      {
        dataTransformed <- dataTransformed[dataids,];
      }
  }
  IDeAEvolution <-list(Corr=correlationMeasureEvolution,Spar=sparcityEvolution/(decordim^2));
  
  
  attr(dataTransformed,"UPSTM") <- DeCorrmatrix;
  attr(dataTransformed,"fscore") <- fscore;
  attr(dataTransformed,"TotalAdjustments") <- countf;
  attr(dataTransformed,"drivingFeatures") <- drivingFeatures;
  attr(dataTransformed,"unaltered") <- bfeat;
  attr(dataTransformed,"LatentVariables") <- lavariables;
  attr(dataTransformed,"useDeCorr") <- useDeCorr;
  attr(dataTransformed,"unipvalue") <- adjunipvalue
  attr(dataTransformed,"totCorrelated") <- totAtcorr
  attr(dataTransformed,"R.critical") <- rcrit
  attr(dataTransformed,"IDeAEvolution") <- IDeAEvolution
  return(dataTransformed)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"UPSTM")))
  {
    decorMat <- attr(decorrelatedobject,"UPSTM")
#    testData[,rownames(decorMat)] <- as.matrix(testData[,rownames(decorMat)]) %*% decorMat
    testData[,rownames(decorMat)] <- Rfast::mat.mult(as.matrix(testData[,rownames(decorMat)]),decorMat);
#    varincluded <- attr(decorrelatedobject,"varincluded")
    newnames <- colnames(testData)
    lavariables <- attr(decorrelatedobject,"LatentVariables");
#    bfeat <- attr(decorrelatedobject,"unaltered")
#    newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
#    newnames[newnames %in% varincluded] <- paste("La_",newnames[newnames %in% varincluded],sep="")
    newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
    colnames(testData) <- newnames
  }
  else
  {
    warning("The object does not have a decorrelation Matrix\n")
  }
  return (testData)
}


getlatentCharFormula <- function(latentCoeff)
{
  deFromula <- character(length(latentCoeff))
  names(deFromula) <- names(latentCoeff)
  for (dx in names(deFromula))
  {
    coef <- latentCoeff[[dx]]
    cname <- names(latentCoeff[[dx]])
    names(cname) <- cname
    for (cf in names(coef))
    {
      if (cf != dx)
      {
        if (coef[cf] != 0)
        {
          if (abs(log10(abs(coef[cf])))<=2)
          {
            if (coef[cf]>0)
            {
              deFromula[dx] <- paste(deFromula[dx],
                                     sprintf("+ (%5.3f)%s",coef[cf],cname[cf]))
            }
            else
            {
              deFromula[dx] <- paste(deFromula[dx],
                                     sprintf("- (%.3f)%s",abs(coef[cf]),cname[cf]))
            }
          }
          else
          {
            if (coef[cf]>0)
            {
              deFromula[dx] <- paste(deFromula[dx],
                                     sprintf("+ (%.2e)%s",coef[cf],cname[cf]))
            }
            else
            {
              deFromula[dx] <- paste(deFromula[dx],
                                     sprintf("- (%.2e)%s",abs(coef[cf]),cname[cf]))
            }
          }
        }
        deFromula[dx] <- str_replace_all(deFromula[dx],"\\(1.000\\)","")
      }
    }
  }
  return(deFromula)
}

getLatentCoefficients <- function(decorrelatedobject)
{
  CoeffList <- list();
  nonZeronames <- character();
  UPSTM <- attr(decorrelatedobject,"UPSTM")
  namecol <- colnames(UPSTM);
  if (length(namecol) > 0)
  {
    n=1;
    for (i in 1:ncol(UPSTM))
    {
      associ <- abs(UPSTM[,i])>0;
      if (sum(associ) > 1)
      {
        nonZeronames <- append(nonZeronames,namecol[i]);
        CoeffList[[n]] <- UPSTM[associ,i];
        n=n+1;
      }
    }
    names(CoeffList) <- nonZeronames;
  }
  attr(CoeffList,"LatentCharFormulas") <- getlatentCharFormula(CoeffList)

  return (CoeffList)
}

## Get the coeficients associated with each one of the observed features
## The latentModel coefficients input should contain the linear coefficients 

getObservedCoef <- function(decorrelatedobject,latentModel)
{
  latentCoef <- latentModel
  Outcome <- NULL;
  if (!is.null(latentModel$formula))
  {
    Outcome <- as.character(latentModel$formula)[2]    
    latentCoef <- latentModel$coef
  }
  

  UPSTM <- attr(decorrelatedobject,"UPSTM")
  laBetas <- numeric(ncol(UPSTM))
  names(laBetas) <- colnames(UPSTM)
  namesin <- names(latentCoef) %in% colnames(UPSTM)
  lanamesinUPSTM <- names(latentCoef)[namesin]
  lanamesNotinUPSTM <- names(latentCoef)[!namesin]
  
  laBetas[lanamesinUPSTM] <- latentCoef[lanamesinUPSTM]
  obsCoef <-  c(latentCoef[lanamesNotinUPSTM],t(UPSTM %*% laBetas))
  names(obsCoef) <- c(lanamesNotinUPSTM,rownames(UPSTM))
  obsCoef <- obsCoef[obsCoef!=0]

  theformula <- str_flatten(names(obsCoef)," + ")
  theformula <- str_replace_all(theformula,"\\(Intercept\\)","1")
  theformula <- paste(Outcome,"~",theformula)
  
  obsCoef <- list(formula=theformula,
                  coefficients=obsCoef,
                  LatentCharFormulas=attr(getLatentCoefficients(decorrelatedobject),"LatentCharFormulas")
                  )
  class(obsCoef) <- c("ObservedModel")
  return (obsCoef)
}