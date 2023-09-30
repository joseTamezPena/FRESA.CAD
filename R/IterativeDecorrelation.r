ILAA <- function(data=NULL,
                thr=0.80,
                method=c("pearson","spearman"),
                Outcome=NULL,
                drivingFeatures=NULL,
                maxLoops=100,
                verbose=FALSE,
                bootstrap=0)
{
  method <- match.arg(method);
  type="LM";
  if (method=="pearson") 
  {
    method <- "fast";
  }
  if (method=="spearman") 
  {
    type="RLM";
  }
  if (verbose)
  {
    cat(method,"|",type,"|")
  }
  
  transf <- IDeA(data=data,
                 thr=thr,
                 method=method,
                 Outcome=Outcome,
                 drivingFeatures=drivingFeatures,
                 maxLoops=maxLoops,
                 verbose=verbose,
                 type=type)
  if (bootstrap > 1)
  {
    
    transform <- attr(transf,"UPLTM") 
    fscore <- attr(transf,"fscore");
    countf <- attr(transf,"TotalAdjustments");
    AdrivingFeatures <- attr(transf,"drivingFeatures");
    bfeat<- attr(transf,"unaltered");
    useDeCorr <- attr(transf,"useDeCorr");
    adjunipvalue <- attr(transf,"unipvalue"); 
    totAtcorr <- attr(transf,"totCorrelated"); 
    rcrit <- attr(transf,"R.critical");
    IDeAEvolution <- attr(transf,"IDeAEvolution"); 

    dsize <- nrow(data);
    lavariables <- names(fscore);
    lavariables <- str_remove_all(lavariables,"La_");
    names(fscore) <- lavariables

    taccmatrix <- diag(length(fscore));
    colnames(taccmatrix) <- lavariables
    rownames(taccmatrix) <- lavariables
    colnames(transform) <- str_remove_all(colnames(transform),"La_");
    taccmatrix[rownames(transform),colnames(transform)] <- transform
    if (verbose) 
    {
      cat("bootstrapping \n")
    }
    
    twts <- 1
    for (lp in c(1:bootstrap))
    {
      if (verbose)
      {         
        cat(".")
        if (lp %% 40 == 0) cat("\n")
      }
      if (bootstrap > 100)  ## True bootstrap
      {
        smp <- sample(dsize,dsize,replace = TRUE);
      }
      else ## pseudo bootrstrap 5% will be duplicated
      {
        smp <- c(sample(dsize,0.95*dsize),sample(dsize,0.10*dsize));
      }
      transf <- IDeA(data[smp,],
                   thr=thr,
                   method=method,
                   Outcome=Outcome,
                   drivingFeatures=AdrivingFeatures,
                   maxLoops=maxLoops,
                   type=type)
      transform <- attr(transf,"UPLTM") 
      colnames(transform) <- str_remove_all(colnames(transform),"La_")
      cnames <- rownames(transform)
      wo <- 1.0;
      if (sum(!(cnames %in% rownames(taccmatrix))) > 0)
      {
        wo <- 0.5;
        if (verbose)
        {
          cat("[",length(cnames[!(cnames %in% rownames(taccmatrix))]),"]")
        }        
        cnames <- cnames[cnames %in% rownames(taccmatrix)];        
      }
      mcor <- min(attr(transf,"IDeAEvolution")$Corr)
      cwt <- (mcor - thr)/(1.0 - thr);
      if (cwt < 0) cwt <- 0
      wt <- wo*(1.0 - cwt)^2
      taccmatrix[cnames,cnames] <- taccmatrix[cnames,cnames] + wt*transform[cnames,cnames]
      twts <- twts + wt
      bscore <- attr(transf,"fscore")
      cnames <- names(bscore)
      cnames <- str_remove_all(cnames,"La_");
      names(bscore) <- cnames
      cnames <- cnames[cnames %in% names(fscore)]
      fscore[cnames] <- fscore[cnames] + wt*bscore[cnames];
      if (verbose)
      {         
        cat(sprintf("(r=%3.2f,w=%3.2f)",mcor,wt));
      }
    }
    transform <- NULL;
    taccmatrix <- taccmatrix/twts;
    fscore <- fscore/twts;
    islatent <- apply(1*(taccmatrix != 0),2,sum) > 1;
    ctnames <- colnames(taccmatrix)
    lavariables <- ctnames[islatent] 
    ctnames[islatent] <- paste("La_",ctnames[islatent],sep="")
    snames <- names(fscore)
    snames[snames %in% lavariables] <- paste("La_",snames[snames %in% lavariables],sep="")
    names(fscore) <- snames
    colnames(taccmatrix) <- ctnames
    attr(transf,"UPLTM") <- taccmatrix

    attr(transf,"LatentVariables") <- lavariables;
    transf <- predictDecorrelate(transf,data)
    attr(transf,"UPLTM") <- taccmatrix
    attr(transf,"fscore") <- fscore;
    attr(transf,"TotalAdjustments") <- countf;
    attr(transf,"drivingFeatures") <- AdrivingFeatures;
    attr(transf,"unaltered") <- bfeat;
    attr(transf,"LatentVariables") <- lavariables;
    attr(transf,"useDeCorr") <- useDeCorr;
    attr(transf,"unipvalue") <- adjunipvalue
    attr(transf,"totCorrelated") <- totAtcorr
    attr(transf,"R.critical") <- rcrit
    attr(transf,"IDeAEvolution") <- IDeAEvolution

    if (verbose) cat("\n\r")
  }  
  
  return(transf)

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
  
  totFeatures <- length(varincluded);
  
  
  models <- NULL
  
  if ((useFastCor) && (length(varincluded)>1000))
  {
    cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded]),large=TRUE))
    colnames(cormat) <- varincluded
    rownames(cormat) <- varincluded
  }
  else
  {
    cormat <- abs(cor(refdata[,varincluded],method=method))
  }
  diag(cormat) <- 0;
  
  adjunipvalue <- min(c(unipvalue,3*unipvalue/totFeatures)); ## Adjusting for false association at 3 true associations per feature
  
  ## Min Correlation based on the Pearson distribution
  ndf <- nrow(data)-2
  tvalue <- qt(1.0 - adjunipvalue,ndf)
  rcrit <- tvalue/sqrt(ndf + tvalue^2)
  if (thr < rcrit) thr <- rcrit

  maxcor <- apply(cormat,2,max)

  totAtcorr <- sum(cormat >= thr); ## 
  DeCorrmatrix <- NULL;
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  bfeat <- NULL

  
  ################################ #######################
#  althr <- 0.25;
#  cortoInclude <- max(c(althr*thr,rcrit))
#  
#  varincluded <- names(maxcor)[maxcor >= cortoInclude];
  
#  if (length(varincluded) > 2000) 
#  { 
#    althr <- 0.5;
#    cortoInclude <- max(c(althr*thr,rcrit))
#    varincluded <- names(maxcor)[maxcor >= cortoInclude];
#    varincluded <- names(maxcor)[maxcor >= rcrit];
#  }    

  allFeatAdded <- character()
  decordim <- 0
#  print(c(nrow(cormat),max(cormat),length(varincluded)))

  if (length(varincluded) > 1)
  {
      if (!is.null(Outcome))
      { 
          if (length(drivingFeatures)==1)          
          {
            if (drivingFeatures[1] == Outcome)
            {
              if (inherits(data[,Outcome],"factor") | (length(unique(data[,Outcome])) == 2))
              {
                outcomep <- univariate_correlation(data,Outcome,limit =-1) # the top associated features to the outcome
              }
              else
              {
                outcomep <- univariate_correlation(data,Outcome,method=method,limit= -1) # the top associated features to the outcome
              }
              drivingFeatures <- names(outcomep);
              if (verbose) print(head(outcomep))
              outcomep <- NULL;
            }
          }
      }
      cormat <- cormat[varincluded,varincluded];
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;
      decordim <- length(varincluded);
      
      maxcor <- apply(cormat,2,max)
      ordera <- maxcor;
      if (corRank)
      {       
        axcor <- cormat;
        axcor[axcor < 0.5*(1.0 + rcrit)] <- 0; ## Remove low correlated features
        ordera <- 0.90*maxcor + 0.10*apply(axcor^2,2,mean); ## Add mean of top correlated
        axcor <- NULL
      }
      
      AdrivingFeatures <- varincluded[order(-ordera[varincluded])];
      if (length(drivingFeatures) > 0)
      {    
        AdrivingFeatures <- unique(c(drivingFeatures,AdrivingFeatures))
        AdrivingFeatures <- AdrivingFeatures[AdrivingFeatures %in% varincluded];
      }
      ordera[AdrivingFeatures] <- c(length(AdrivingFeatures):1)
      ordera <- ordera/length(AdrivingFeatures);
      bfeat <- as.character(correlated_Remove(cormat,AdrivingFeatures,thr = thr,isDataCorMatrix=TRUE))
      

      if (verbose) 
      {
        cat ("\n",head(bfeat),"\n");
        print(head(ordera));
        cat ("\n Included:",length(varincluded),
                        ", Uni p:",adjunipvalue,
                        ", Base Size:",length(bfeat),
                        ", Rcrit:",rcrit,
                        "\n")
      }

      thr2 <- thr
      thr <- 1.1*thr
      thrvalues <- c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05);
      nextthr <- 1;
      nextthra <- 0;
      while (((addedlist > 0) || (thr > thr2)) && (lp < maxLoops)) 
      {
        lp = lp + 1;

        addedlist <- 0;
        
        maxcor <- apply(cormat,2,max)
        mxScor <- max(maxcor);
        correlationMeasureEvolution <- c(correlationMeasureEvolution,mxScor);
        sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

        thr <- thr2;
        if ((mxScor >= thr2) && relaxed)
        {
          nextthr <- sum(thrvalues > mxScor) + 1
          thr <- max(thr2,thrvalues[nextthr]);
        }

        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (verbose)  cat("\n\r",lp,sprintf("<R=%5.3f,thr=%5.3f>",mxScor,thr))
        maxcomat <- NULL
        if (length(topfeat)>0)
        {
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          topfeat <- topfeat[order(-ordera[topfeat])]; 
          atopbase <- bfeat[bfeat %in% topfeat];
          topfeat <- unique(c(atopbase,topfeat));
          topfeat <- correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE)
          topfeat <- topfeat[order(-maxcor[topfeat])];
          intopfeat <- character();
          toBeDecorrelated <- length(topfeat)
          if (verbose)  cat(", Top:",toBeDecorrelated);
          featAdded <- character();
          if (toBeDecorrelated > 0)
          {
            if (toBeDecorrelated > 1)
            {
              maxcomat <- apply(cormat[,topfeat],1,max);
            }
            for (feat in topfeat)
            {
                corlist <- cormat[,feat];
                corlist <- corlist[corlist >= thr];
                varlist <- names(corlist)

                if ((toBeDecorrelated > 1) && (length(varlist) > 1))
                {
                  varlist <- names(corlist[corlist >= maxcomat[varlist]])
                }
                
                varlist <- varlist[!(varlist %in% decorrelatedFetureList)]

                if (length(varlist) > 0)
                {
                   if (verbose && (feat==topfeat[1]))  cat("<",length(varlist),">");
                   if (useFastCor)
                   {
                      prebetas <- getAllBetaCoefficients(feat,varlist);
                      varlist <- varlist[prebetas != 0];
                      if (length(varlist) > 0)
                      {
                        betamatrix[feat,varlist] <- -1.0*prebetas[prebetas != 0];
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
                              dataTransformed[,c(feat,varlist)] <- adataframe[,c(feat,varlist)];
                              refdata[,c(feat,varlist)] <- adataframe[refdataids,c(feat,varlist)];
                          }
                      }
                    }
                    if (length(varlist) > 0)
                    {
                        featAdded <- c(featAdded,feat);
                        intopfeat <- unique(c(intopfeat,feat));
                        countf[varlist] <- countf[varlist] + 1;
                        decorrelatedFetureList <- c(decorrelatedFetureList,varlist);
                        fscore[feat] <- fscore[feat] + length(varlist);
                        fscore[varlist] <- fscore[varlist] - 1;
                        if (toBeDecorrelated > 1)
                        {
                          maxcomat[varlist] <- 0;
                        }
                        if (verbose && (length(intopfeat) %% 100 == 99)) cat(".")
                    }
                }
            }
            allFeatAdded <- unique(c(allFeatAdded,featAdded))
            ordera[featAdded] <- ordera[featAdded] + 1;
          }
          addedlist <- length(decorrelatedFetureList);
          if (verbose) cat("[Fa=",length(allFeatAdded),"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
          colused <- character()
          if (addedlist > 0)
          {
             mbetas <- c(intopfeat,decorrelatedFetureList);
             totused <- unique(c(totused,mbetas));
             totalpha <- unique(c(totalpha,intopfeat));
             allused <- varincluded[varincluded %in% totused];
             colused <- varincluded[varincluded %in% decorrelatedFetureList];
             alphaused <- varincluded[varincluded %in% intopfeat];
             totalphaused <- varincluded[varincluded %in% totalpha];
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
          if (verbose) cat("><")
          if (length(colused) > 0)
          {
            colsd <- apply(as.matrix(refdata[,colused]),2,sd,na.rm = TRUE);
            if (sum(colsd==0) > 0)
            {
              zerovar <- colused[colsd==0];
              for (zcheck in zerovar)
              {
                refdata[,zcheck] <- refdata[,zcheck] + rnorm(nrow(refdata),0,1.0e-10);
              }
            }
          }
          topFeatures <- unique(c(topFeatures,intopfeat));
          if ((useFastCor) && (length(colused)>1000))
          {
            cormat <- abs(Rfast::cora(as.matrix(refdata[,varincluded]),large=TRUE))
            colnames(cormat) <- varincluded
            rownames(cormat) <- varincluded
          }
          else
          {
            if (length(colused) > 0)
            {
              cormatS <- abs(cor(refdata[,varincluded],refdata[,colused],method=method))
              cormat[,colused] <- cormatS
              cormat[colused,] <- t(cormatS)
              cormatS <- NULL
            }
          }
          if (verbose) cat(">")
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
        mxScor <- max(cormat);
        nextthr <- sum(thrvalues > mxScor) + 1

        if ((nextthr == nextthra) && (addedlist == 0) && (mxScor > thr))
        {
            lp <- maxLoops
        }
        if ((addedlist == 0) && (mxScor > thr))
        {
          nextthra <- nextthr
        }        
      }
      if (lp >= maxLoops)
      {
        warning("Maximum number of iterations reached: Failed to achieve target correlation\n")
      }
      correlationMeasureEvolution <- c(correlationMeasureEvolution,max(cormat));
      sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

      betamatrix <- NULL;
      tmparincluded <- varincluded;
      varincluded <- tmparincluded[tmparincluded %in% totused];
      bfeat <- bfeat[bfeat %in% varincluded];
      if (useDeCorr)
      {
        DeCorrmatrix <- DeCorrmatrix[varincluded,varincluded]
        abfeat <- varincluded[apply(DeCorrmatrix!=0,2,sum)==1]
        bfeat <- unique(c(bfeat,abfeat));
        dataTransformed <- data
        if (length(varincluded) > 1)
        {
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
             if ((useFastCor) && (length(varincluded)>1000))
             {
              cormat <- abs(Rfast::cora(as.matrix(dataTransformed[,varincluded]),large=TRUE))
              colnames(cormat) <- varincluded
              rownames(cormat) <- varincluded
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
#          fscore <- fscore[varincluded];

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
  
  
  attr(dataTransformed,"UPLTM") <- DeCorrmatrix;
  attr(dataTransformed,"fscore") <- fscore;
  attr(dataTransformed,"TotalAdjustments") <- countf;
  attr(dataTransformed,"drivingFeatures") <- AdrivingFeatures;
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
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"UPLTM")))
  {
    decorMat <- attr(decorrelatedobject,"UPLTM")
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
  UPLTM <- attr(decorrelatedobject,"UPLTM")
  namecol <- colnames(UPLTM);
  if (length(namecol) > 0)
  {
    n=1;
    for (i in 1:ncol(UPLTM))
    {
      associ <- abs(UPLTM[,i])>0;
      if (sum(associ) > 1)
      {
        nonZeronames <- append(nonZeronames,namecol[i]);
        CoeffList[[n]] <- UPLTM[associ,i];
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
  

  UPLTM <- attr(decorrelatedobject,"UPLTM")
  laBetas <- numeric(ncol(UPLTM))
  names(laBetas) <- colnames(UPLTM)
  namesin <- names(latentCoef) %in% colnames(UPLTM)
  lanamesinUPLTM <- names(latentCoef)[namesin]
  lanamesNotinUPLTM <- names(latentCoef)[!namesin]
  
  laBetas[lanamesinUPLTM] <- latentCoef[lanamesinUPLTM]
  obsCoef <-  c(latentCoef[lanamesNotinUPLTM],t(UPLTM %*% laBetas))
  names(obsCoef) <- c(lanamesNotinUPLTM,rownames(UPLTM))
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