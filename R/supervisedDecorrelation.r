GDSTMDecorrelation <- function(data=NULL,
                                  thr=0.80,
                                  method=c("fast","pearson","spearman","kendall"),
                                  Outcome=NULL,
                                  refdata=NULL,
                                  drivingFeatures=NULL,
                                  useDeCorr=TRUE,
                                  verbose=FALSE,
                                  relaxed=TRUE,
                                  corRank=FALSE,
                                  maxLoops=100,
                                  unipvalue=0.05,
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
  bfeat <- unique(c(drivingFeatures,AdrivingFeatures))
  
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
#  fcount <- fscore;
  
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
  totcorr <- sum(cormat >= 0.5); ## For False Discovery estimation
  varincluded <- names(maxcor)[maxcor >= 0.35*thr];
  DeCorrmatrix <- NULL;
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  if (length(varincluded) > 1)
  {
      if (!is.null(Outcome))
      { 
          if (length(drivingFeatures)==0) 
          {
            if ((class(data[,Outcome]) == "factor") | (length(unique(data[,Outcome])) == 2))
            {
              outcomep <- univariate_correlation(data,Outcome,limit=0,pvalue=unipvalue,thr = sqrt(thr)) # the top associated features to the outcome
            }
            else
            {
              outcomep <- univariate_correlation(data,Outcome,method=method,limit=0,pvalue=unipvalue,thr = sqrt(thr)) # the top associated features to the outcome
            }
            drivingFeatures <- names(outcomep);
            outcomep <- NULL;
          }
      }
      cormat <- cormat[varincluded,varincluded];
      cormat[cormat < thr] <- 0;
      maxcor <- apply(cormat,2,max)
      AdrivingFeatures <- varincluded;
      names(AdrivingFeatures) <- AdrivingFeatures;
      ordcor <- maxcor;
      if (corRank)
      {       
        ordcor <- apply(cormat^2,2,sum);
      }
      AdrivingFeatures <- AdrivingFeatures[order(-ordcor[AdrivingFeatures])];
      AdrivingFeatures <- as.character(correlated_Remove(cormat,AdrivingFeatures,thr = 0.95*thr,isDataCorMatrix=TRUE))
      unipvalue <- min(unipvalue,unipvalue*sqrt(totcorr)/totFeatures); ## Adjusting for false association
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;
      bfeat <- AdrivingFeatures
      if (length(drivingFeatures) > 0)
      {    
        drivingFeatures <- drivingFeatures[drivingFeatures %in% varincluded];
        if (length(drivingFeatures)>1)
        {
          drivingFeatures <- drivingFeatures[order(-ordcor[drivingFeatures])];
          names(drivingFeatures) <- drivingFeatures;
          bfeat <- unique(c(drivingFeatures,AdrivingFeatures))
          names(bfeat) <- bfeat;
          bfeat <- as.character(correlated_Remove(cormat,bfeat,thr = 0.95*thr,isDataCorMatrix=TRUE));
          bfeat <- bfeat[order(-ordcor[bfeat])];
          drivingFeatures <- drivingFeatures[drivingFeatures %in% bfeat];
          AdrivingFeatures <- AdrivingFeatures[AdrivingFeatures %in% bfeat];
        }
      }
      
      if (verbose) cat ("\n Included:",length(varincluded),", Uni p:",unipvalue,", Uncorrelated Base:",length(AdrivingFeatures),", Outcome-Driven Size:",length(drivingFeatures),", Base Size:",length(bfeat),"\n")

      relaxalpha <- 0.80;
      relaxbeta <- 0.20;
      thr2 <- thr
      thr <- thr2*1.001;
#      wthr <- relaxalpha - relaxalpha*(!relaxed);
      thrvalues <- c(0.95,0.90,0.80,0.70,0.50,0.25,0.10);
      nextthr <- 1 + (!relaxed)*length(thrvalues);
      nextthra <- 0;
      while (((addedlist > 0) || (thr > (thr2*1.0001))) && (lp < maxLoops)) 
      {
        lp = lp + 1;

        addedlist <- 0;

        
        cormat[cormat < thr2] <- 0;
        maxcor <- apply(cormat,2,max)
        mxScor <- max(maxcor);


        if ((mxScor >= thr2) && (nextthr <= length(thrvalues)))
        {
#          thr <- max(c(thr2,thr2 + wthr*(max(maxcor) - thr2)));
           thr <- max(c(thr2,0.5*thrvalues[nextthr] + 0.5*mxScor));
        }
        else
        {
          thr <- thr2;
        }
#        if (verbose)  cat("\n\r",lp,sprintf("<R=%5.3f,w=%5.3f,thr=%5.3f>",max(maxcor),wthr,thr))

        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        if (nextthra != nextthr)
        {
          orglentopfeat <- sum(maxcor[topfeat] >= thr);
          nextthra <- nextthr;
        }
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (verbose)  cat("\n\r",lp,sprintf("<R=%5.3f,w=%3d,N=%5d>",mxScor,nextthr,orglentopfeat))

        if (length(topfeat)>0)
        {
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          ordcor <- maxcor;
          if (corRank)
          {
            ordcor <- apply(cormat^2,2,sum);
          }
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
                    fscore[feat] <- fscore[feat] + length(varlist);
                    fscore[varlist] <- fscore[varlist] - 1;
#                    fscore[feat] <- fscore[feat] + sum(1.0-cormat[varlist,feat]^2)
#                    fscore[varlist] <- fscore[varlist] - cormat[varlist,feat]^2
#                    fcount[feat] <- fcount[feat] + length(varlist);
#                    fcount[varlist] <- fcount[varlist] + 1;
                    cormat[ovarlist,feat] <- 0;

                  }
                  if (allvartested)
                  {
                      featMarked <- unique(c(featMarked,feat));
                  }
              }
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
          if (verbose) cat("[",lpct,":",length(featMarked),":",sprintf("%5.3f",themaxthr),"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
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

        if ((nextthr <= length(thrvalues)) && ((addedlist <= max(c(2.0,0.1*orglentopfeat))) || (length(intopfeat) <= 1)))
        {
#          wthr <- wthr - relaxbeta;
#          if (verbose) 
#          {
#            cat ("orglentopfeat:",orglentopfeat,", addedlist:",addedlist,", intopfeat:",length(intopfeat));
#          }
          nextthr <- nextthr + 1; 
        }
      }
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
        if (verbose) cat("-{",abfeat[!(abfeat %in% bfeat)],"}-")
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
              cat ("\n\r [",lp,"],",max(cormat),"Decor Dimension:",length(varincluded),". Cor to Base:",length(correlatedToBase),", ABase:",length(AdrivingFeatures),", Outcome Base:",length(drivingFeatures),"\n\r")
          }
#          banames <- colnames(DeCorrmatrix)[apply(DeCorrmatrix!=0,2,sum)==1]
          changed <- (apply(DeCorrmatrix!=0,2,sum) > 1);
          lavariables <- colnames(DeCorrmatrix)[changed]
          newnames <- colnames(dataTransformed)
#          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
#          newnames[newnames %in% varincluded] <- paste("La_",newnames[newnames %in% varincluded],sep="")
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(dataTransformed) <- newnames
          newnames <- colnames(DeCorrmatrix)
#          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
#          newnames[newnames %in% varincluded] <- paste("La_",newnames[newnames %in% varincluded],sep="")
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(DeCorrmatrix) <- newnames
          fscore <- fscore[varincluded];
#          fcount <- fcount[varincluded];
#          fscore <- fscore/fcount;

          newnames <- names(fscore)
#          newnames[newnames %in% bfeat] <- paste("Ba_",newnames[newnames %in% bfeat],sep="") 
#          newnames[newnames %in% varincluded] <- paste("La_",newnames[newnames %in% varincluded],sep="")
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          names(fscore) <- newnames
        }
      }
      else
      {
        dataTransformed <- dataTransformed[dataids,];
      }
  }
  
  
  attr(dataTransformed,"GDSTM") <- DeCorrmatrix;
  attr(dataTransformed,"fscore") <- fscore;
  attr(dataTransformed,"TotalAdjustments") <- countf;
  attr(dataTransformed,"drivingFeatures") <- drivingFeatures;
  attr(dataTransformed,"unaltered") <- bfeat;
  attr(dataTransformed,"LatentVariables") <- lavariables;
  attr(dataTransformed,"useDeCorr") <- useDeCorr;
  attr(dataTransformed,"unipvalue") <- unipvalue
  return(dataTransformed)
}

predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"GDSTM")))
  {
    decorMat <- attr(decorrelatedobject,"GDSTM")
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

getLatentCoefficients <- function(decorrelatedobject)
{
  CoeffList <- list();
  nonZeronames <- character();
  GDSTM <- attr(decorrelatedobject,"GDSTM")
  namecol <- colnames(GDSTM);
  if (length(namecol) > 0)
  {
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
  }
  
  return (CoeffList)
}
