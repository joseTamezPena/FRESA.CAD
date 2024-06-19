## ILAA is  a wrapper of IDeA that computes a linear transformation matrix and returns the transformed data set.
# The transformation matrix is an attribute of the returned data frame
# The main motivation behind ILAA is to simplify IDeA calling and to return robust transformation based on small data perturbations.

# ### ILAA Function Interface

# The `ILAA` function is designed to process and analyze a dataset (`data`) using specified thresholds and methods. It allows for bootstrapping to enhance the #robustness of the analysis and produces various metrics and adjustments based on the input parameters. Below is a detailed description of each parameter and #the function's overall operation:

# #### Parameters:

# -data: `data.frame` (default: `NULL`)
# - The input dataset to be analyzed. This should be a data frame containing the variables of interest.

# -thr: `numeric` (default: `0.80`)
# - Threshold value used in the analysis. This typically represents a correlation threshold or a cutoff value for some statistical measure.

# -method: `character` (default: `c("pearson", "spearman")`)
# - The correlation method to be used in the analysis. Options are "pearson" for Pearson correlation and "spearman" for Spearman rank correlation.

# -Outcome: `character` or `NULL` (default: `NULL`)
# - The outcome variable in the dataset. If provided, the analysis will focus on this outcome.

# -drivingFeatures: `character` vector or `NULL` (default: `NULL`)
# - Specific features in the dataset that are considered as driving factors for the analysis. If provided, these features will be given special consideration in the analysis.

# -maxLoops: `integer` (default: `100`)
# - The maximum number of loops or iterations to be performed in the analysis process.

# -verbose: `logical` (default: `FALSE`)
# - If `TRUE`, the function will print detailed output during execution, which is useful for debugging or understanding the analysis steps.

# -bootstrap: `integer` (default: `0`)
# - The number of bootstrap samples to be used. If `bootstrap` is greater than 1, the function will perform bootstrapping to enhance the robustness of the results.

# #### Operation:

# 1. Method Selection:
 # - The function starts by selecting the correlation method based on the input parameter (`method`). It sets the `type` to "LM" for Pearson and "RLM" for Spearman correlations.

# 2. Initial IDeA Call:
 # - The `IDeA` function is called with the specified parameters to perform the initial analysis and transformation of the data.

# 3. Bootstrapping (if applicable):
 # - If `bootstrap` is greater than 1, the function initializes matrices and vectors to store results from multiple bootstrap samples.
 # - It performs bootstrapping by resampling the data and calling the `IDeA` function repeatedly, adjusting and accumulating the results across samples.

# 4. Result Aggregation:
 # - The function aggregates the results from the bootstrapping process, normalizing the accumulated transformation matrices and scores.
 # - It identifies latent variables and adjusts the transformation matrices accordingly.

# 5. Final Adjustments:
 # - The final transformation matrix and other attributes are updated based on the aggregated results.
 # - The function applies decorrelation and calculates variance ratios for the transformed and original data.

# 6. Return:
 # - The function returns the final transformation results (`transf`) with updated attributes, including transformation matrices, scores, latent variables, and variance ratios.

# ---


ILAA <- function(data=NULL,
                thr=0.80,
                method=c("pearson","spearman"),
                Outcome=NULL,
                drivingFeatures=NULL,
                maxLoops=100,
                verbose=FALSE,
                bootstrap=0)
{
# Select method (either "pearson" or "spearman") and set type accordingly
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
  ## First IDeA Call
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
    # Initialize various matrices and vectors for bootstrap
    transform <- attr(transf,"UPLTM")
    ocalnames <- colnames(transform)
    fscore <- attr(transf,"fscore");
    countf <- attr(transf,"TotalAdjustments");
    AdrivingFeatures <- attr(transf,"drivingFeatures");
    bfeat<- attr(transf,"unaltered");
    adjunipvalue <- attr(transf,"unipvalue"); 
    rcrit <- attr(transf,"R.critical");
    IDeAEvolution <- attr(transf,"IDeAEvolution"); 
    useDeCorr <- attr(transf,"useDeCorr")
    

    dsize <- nrow(data);
    lavariables <- names(fscore);
    lavariables <- str_remove_all(lavariables,"La_");
    names(fscore) <- lavariables

    taccmatrix <- diag(length(fscore));
    colnames(taccmatrix) <- lavariables
    rownames(taccmatrix) <- lavariables
    twts <- rep(1,nrow(taccmatrix));
    names(twts) <- lavariables
    colnames(transform) <- str_remove_all(ocalnames,"La_");

    taccmatrix[rownames(transform),colnames(transform)] <- transform
    sgnmatrix <- 1*(taccmatrix > 0) - 1*(taccmatrix < 0)
    tmptransform <- sgnmatrix
    if (verbose) 
    {
      cat("bootstrapping \n")
    }
    
    for (lp in c(1:bootstrap))
    {
      if (verbose)
      {         
        cat(".")
        if (lp %% 40 == 0) cat("\n")
      }
      if (bootstrap > 500)  ## True bootstrap
      {
        smp <- sample(dsize,dsize,replace = TRUE);
      }
      else ## small bootstrap 5% will be duplicated
      {
        smp <- c(sample(dsize,0.95*dsize),sample(dsize,0.05*dsize));
      }
      ## Calling IDeA inside the Loop
      transf <- IDeA(data[smp,],
                   thr=thr,
                   method=method,
                   Outcome=Outcome,
                   drivingFeatures=AdrivingFeatures,
                   maxLoops=maxLoops,
                   type=type)
      ## Norw we will aggregate the ILAA transformations
      transform <- attr(transf,"UPLTM")
      bcolnames <- colnames(transform)     
      colnames(transform) <- str_remove_all(bcolnames,"La_")
      rnames <- rownames(transform)
      cnames <- rnames 
      # wo for weighed aggregation
      wo <- 1.0;
      tmptransform <- sgnmatrix
      tmptransform[rnames,cnames] <- 1*(transform > 0) - 1*(transform < 0)
      sgnmatrix[sgnmatrix == 0] <- tmptransform[sgnmatrix == 0]
      sgchange <- apply(1*(sgnmatrix*tmptransform < 0),2,sum) 
      if ((sum(!(bcolnames %in% ocalnames)) > 0) | (sum(sgchange) > 0))
      {
        if (verbose)
        {
          cat("{",length(bcolnames[!(bcolnames %in% ocalnames)]),",",sum(sgchange),"}")
        }        
        wo <- 0.1;
        cnames <- bcolnames[bcolnames %in% ocalnames];
        cnames <- str_remove_all(cnames,"La_")
        cnames <- cnames[cnames %in% names(sgchange)[sgchange == 0]]
      }
      ## Estimating the weights for the aggregation. It is an ad hoc weights 
      mcor <- min(attr(transf,"IDeAEvolution")$Corr)
      cwt <- (mcor - thr)/(1.0 - thr);
      if (cwt < 0) cwt <- 0
      wt <- wo*(1.0 - cwt)^2
      ## Weighted aggregation. 
      taccmatrix[rnames,cnames] <- taccmatrix[rnames,cnames] + wt*transform[rnames,cnames]
      bscore <- attr(transf,"fscore")
      names(bscore) <- str_remove_all(names(bscore),"La_");
      fscore[cnames] <- fscore[cnames] + wt*bscore[cnames];
      twts[cnames] <- twts[cnames] + wt
      if (verbose)
      {         
        cat(sprintf("(r=%3.2f,w=%3.2f)",mcor,wt));
      }
    }
    ## Updated the transformation and the discovered latent variables
    tmptransform <- NULL
    sgnmatrix <- NULL
    transform <- NULL;
    for (cn in colnames(taccmatrix))
    {
      taccmatrix[,cn] <- taccmatrix[,cn]/twts[cn];
    }

    fscore <- fscore/twts;
    islatent <- apply(1*(taccmatrix != 0),2,sum) > 1;
    ctnames <- colnames(taccmatrix)
    lavariables <- ctnames[islatent] 
    ctnames[islatent] <- paste("La_",ctnames[islatent],sep="")
    snames <- names(fscore)
    snames[snames %in% lavariables] <- paste("La_",snames[snames %in% lavariables],sep="")
    names(fscore) <- snames
    colnames(taccmatrix) <- ctnames
    tokeep <- !((apply(1*(taccmatrix != 0),2,sum) == 1) & (apply(1*(taccmatrix != 0),1,sum) == 1))
    taccmatrix <- taccmatrix[tokeep,tokeep]

    
    attr(transf,"UPLTM") <- taccmatrix
    attr(transf,"LatentVariables") <- lavariables;
    ## Transform the input dataframe
    transf <- predictDecorrelate(transf,data)
    
    VarTrans <- apply(transf[,colnames(taccmatrix)],2,var);
    VarObs <- apply(data[,rownames(taccmatrix)],2,var);
    VarRatio <- VarTrans/VarObs
    VarRatio <- VarRatio[order(-VarRatio)]
    noaltered <- colnames(data)[!(colnames(data) %in% rownames(taccmatrix))]
    if (length(noaltered)>0)
    {
      varone <- rep(1,length(noaltered))
      names(varone) <- noaltered
      VarRatio <-  c(varone,VarRatio)
    }

    ## Setting the attributes of the returned dataframe
    attr(transf,"UPLTM") <- taccmatrix
    attr(transf,"fscore") <- fscore;
    attr(transf,"TotalAdjustments") <- countf;
    attr(transf,"drivingFeatures") <- AdrivingFeatures;
    attr(transf,"unaltered") <- bfeat;
    attr(transf,"LatentVariables") <- lavariables;
    attr(transf,"useDeCorr") <- useDeCorr;
    attr(transf,"unipvalue") <- adjunipvalue
    attr(transf,"R.critical") <- rcrit
    attr(transf,"IDeAEvolution") <- IDeAEvolution
    attr(transf,"VarRatio") <- VarRatio

    if (verbose) 
    {
      cat("\n\r")
      print(head(VarRatio))
    }
  }  
  
  return(transf)

}                                  



## IDeA Computes the linear transformation that mitigates feature collinearity. The transformation is then used to transform the input database.
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

# Ensure 'Rfast' package is installed
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

# Function to calculate linear beta coefficients for a given feature: x_i=b*x_j+c
  getAllBetaCoefficients <- function(feat,varlist=NULL)
  {
    allBetaCoef <- numeric(length(varlist));
    featVector <- refdata[,feat]

    # Fit the feature to another feature and returns the betaCoef
    getBetaCoefficient <- function(x)
    {
      betaCoef <- 0;
      modellm <- try(lm(x~featVector,model = FALSE,na.action=na.exclude))
      if (!inherits(modellm, "try-error"))
      {	
        f <- summary(modellm)$coefficients
        # If not significant associatoin set it to zero
        if (nrow(f)>1)
        {
          betaCoef <- modellm$coef[2]*(f[2,4] < adjunipvalue);
        }
      }
      return (betaCoef)
    }
    # Apply the beta-estimation to all features
    if (length(varlist) > 1)
    {
        allBetaCoef <- apply(refdata[,varlist],2,getBetaCoefficient)
    }
    else
    {
      allBetaCoef <- getBetaCoefficient(refdata[,varlist])
    }
    # If error set to zero
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

  # Initializing variables
  
  DeCorrmatrix <- NULL;
  topFeatures <- character()
  correlatedToBase <- character();
  lavariables <- character();
  lastintopfeat <- character();
  lastdecorrelated <- character();
  totused <- character();
  totalpha <- character();
  bfeat <- NULL



  if (is.null(refdata))
  {
    refdata <- data;
  }
  
  dataids <- rownames(data)
  refdataids <- rownames(refdata);
  dataTransformed <- as.data.frame(rbind(data,refdata[!(refdataids %in% dataids),]));

  
  
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
  # Variables to analyze for collinearity
  varincluded <- colnames(refdata);
  # If Outcome is present remove it from the variable list
  if (!is.null(Outcome))
  {
    varincluded <- varincluded[!(varincluded %in% Outcome)];
  }

  # The character features
  ischaracter <- sapply(refdata[,varincluded],class) == "character";
  # The factors from the list
  isFactor <- sapply(refdata[,varincluded],class) == "factor";
  # Variables with not to many different values
  isContinous <- unlist(lapply(lapply(refdata[,varincluded],unique),length)) >= minUniqueValues; 

  #The final list of variables that may have collinearity issues
  varincluded <- varincluded[!isFactor & isContinous & !ischaracter];
  
  
  fscore <- numeric(length(varincluded));
  names(fscore) <- varincluded;
  
  totFeatures <- length(varincluded);
  
  
  models <- NULL
  
  # Compute the initial correlation matrix. Rfast correlation is used in large matrices.
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
  # Set to zero any NA
  cormat[is.na(cormat)] <- 0
  # Set to zero the diagonal
  diag(cormat) <- 0;
  
  # This may be debatable. But lets use this as an add hoc value for the estimation of the critical threshold 
  adjunipvalue <- min(c(unipvalue,3*unipvalue/totFeatures)); ## Adjusting for false association at 3 true associations per feature
  
  # Calculate critical correlation threshold
  # Any correlation below this threshold is not significant
  # If the specified thr is lower than the critical value thr is set to the critical value  
  ndf <- nrow(data)-2
  tvalue <- qt(1.0 - adjunipvalue,ndf)
  rcrit <- tvalue/sqrt(ndf + tvalue^2)
  if (thr < rcrit) thr <- rcrit

  # Initial values of variables
  maxcor <- apply(cormat,2,max)


  

  allFeatAdded <- character()
  decordim <- 0

  if (length(varincluded) > 1)
  {
    # If an outcome is specified, Let us get the variables associated to the outcome and use them as driving features
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
      # Only continuous variables 
      cormat <- cormat[varincluded,varincluded];
      # Initial shape of the decorrelation matrix
      DeCorrmatrix <- diag(length(varincluded));
      colnames(DeCorrmatrix) <- varincluded;
      rownames(DeCorrmatrix) <- varincluded;
      decordim <- length(varincluded);
      
      maxcor <- apply(cormat,2,max)
      # Ranking the variables accordingly to their maximum association strength
      ordera <- maxcor;
      if (corRank)
      { 
        # If the corRank is set the order now depends on the average correlation of top associated variables
        axcor <- cormat;
        axcor[axcor < 0.5*(1.0 + rcrit)] <- 0; ## Remove low correlated features
        ordera <- 0.90*maxcor + 0.10*apply(axcor^2,2,mean); ## Add mean of top correlated
        axcor <- NULL
      }
      
      AdrivingFeatures <- varincluded[order(-ordera[varincluded])];
      if (length(drivingFeatures) > 0)
      { 
        # The Adrivingfeatures define the order of association exploration. i.e., we want them to be unaltered basis  
        AdrivingFeatures <- unique(c(drivingFeatures,AdrivingFeatures))
        AdrivingFeatures <- AdrivingFeatures[AdrivingFeatures %in% varincluded];
      }
      # We use a small trick to force the AdrivingFeatures to be at the top of the list of unaltered basis candidates.
      ordera[AdrivingFeatures] <- c(length(AdrivingFeatures):1) 
      ordera <- ordera/length(AdrivingFeatures);
      
      # Only features whose association is larger than the critical values will be analyzed
      bfeat <- as.character(correlated_Remove(cormat,AdrivingFeatures,thr = rcrit,isDataCorMatrix=TRUE))
      

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

      thr2 <- thr # Storing the user specified target
      thr <- 1.1*thr
      # The algorithm will may loop using the following threshold values
      thrvalues <- c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05);
      nextthr <- 1;
      nextthra <- 0;
      skv <- 1
      # Main loop for feature selection and decorrelation
      while (((addedlist > 0) || (thr > thr2)) && (lp < maxLoops)) 
      {
        lp = lp + 1;

        addedlist <- 0;
        
        #The maximum correlation
        maxcor <- apply(cormat,2,max)
        mxScor <- max(maxcor); # The current maximum correlation in the current transformed dataset
        mxScorm <- mean(c(maxcor[maxcor>=thr2],mxScor));
        correlationMeasureEvolution <- c(correlationMeasureEvolution,mxScorm);
        sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

        #Adjusting the loop correlation target
        thr <- thr2;
        if ((mxScor >= thr2) && relaxed)
        {
          # Adjust current loop thr if max correlation is still greater than the target thr.
          nextthr <- sum(thrvalues > mxScor) + skv
          thr <- max(thr2,thrvalues[nextthr]);
        }

        # Only features that have a correlation greater than "thr" will be decorrelated
        topfeat <- varincluded;
        names(topfeat) <- topfeat;
        topfeat <- topfeat[maxcor[topfeat] >= thr];
        if (verbose)  cat("\n\r",lp,sprintf("<R=%5.3f,thr=%5.3f>",mxScorm,thr))
        maxcomat <- NULL
        if (length(topfeat)>0)
        {
          # Only if significant associations still exists
          decorrelatedFetureList <- character();
          betamatrix <- diag(length(varincluded));
          colnames(betamatrix) <- varincluded;
          rownames(betamatrix) <- varincluded;
          # Get the minimum subset of features to be decorrelated
          topfeat <- topfeat[order(-ordera[topfeat])]; 
          atopbase <- bfeat[bfeat %in% topfeat];
          topfeat <- unique(c(atopbase,topfeat));
          topfeat <- correlated_Remove(cormat,topfeat,thr = thr,isDataCorMatrix=TRUE)
          topfeat <- topfeat[order(-maxcor[topfeat])]; #topfeat has driving features still associated to other features
          
          intopfeat <- character();
          toBeDecorrelated <- length(topfeat)
          if (verbose)  cat(", Top:",toBeDecorrelated);
          featAdded <- character();
          
          if (toBeDecorrelated > 0)
          {
            # Only if there are features with significant correlation
            if (toBeDecorrelated > 1)
            {
              maxcomat <- apply(cormat[,topfeat],1,max);
            }
            for (feat in topfeat)
            {
              # Loop over all driving features associated with other features
                corlist <- cormat[,feat];
                corlist <- corlist[corlist >= thr];
                # Get the features that are correlated to the driving features
                varlist <- names(corlist)

                if ((toBeDecorrelated > 1) && (length(varlist) > 1))
                {
                  varlist <- names(corlist[corlist >= maxcomat[varlist]])
                }
                
                varlist <- varlist[!(varlist %in% decorrelatedFetureList)]

                if (length(varlist) > 0)
                {
                  # Get the beta coefficients of the driving features to their corresponding associated pair
                   if (verbose && (feat==topfeat[1]))  cat("<",length(varlist),">");
                   if (useFastCor)
                   {
                      #Linear Method
                      prebetas <- getAllBetaCoefficients(feat,varlist);
                      varlist <- varlist[prebetas != 0];
                      if (length(varlist) > 0)
                      {
                        betamatrix[feat,varlist] <- -1.0*prebetas[prebetas != 0];
                      }
                   }
                   else
                   {
                    #Robust methods
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
                      # If any variable was adjusted, we do inloop  bookkeeping
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
            # Store added features
            allFeatAdded <- unique(c(allFeatAdded,featAdded))
            ordera[featAdded] <- ordera[featAdded] + 1;
          }
          addedlist <- length(decorrelatedFetureList);
          if (verbose) cat("[Fa=",length(allFeatAdded),"](",length(intopfeat),",",addedlist,",",length(totalpha),"),<")
          colused <- character()
          if (addedlist > 0)
          {
            # If features were decorrelated we do the new estimation of the transformation matrix
             mbetas <- c(intopfeat,decorrelatedFetureList);
             totused <- unique(c(totused,mbetas));
             totalpha <- unique(c(totalpha,intopfeat));
             allused <- varincluded[varincluded %in% totused];
             colused <- varincluded[varincluded %in% decorrelatedFetureList];
             alphaused <- varincluded[varincluded %in% intopfeat];
             totalphaused <- varincluded[varincluded %in% totalpha];
              if ((length(colused) == 1) && (length(allused) == 1))
              {
                # Updating the transformation matrix for a single variable change
                DeCorrmatrix[totalphaused,colused] <-  as.numeric(DeCorrmatrix[totalphaused,allused]) *
                                                               as.numeric(betamatrix[allused,colused]);
              }
              else
              {
                # Updating the transformation matrix for more than one change
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
                # Update the transformed data set
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
              # We don´t like zero residuals. i.e., we add some minor noise
              zerovar <- colused[colsd==0];
              for (zcheck in zerovar)
              {
                refdata[,zcheck] <- refdata[,zcheck] + rnorm(nrow(refdata),0,1.0e-10);
              }
            }
          }
          # Recompute correlation matrix
          topFeatures <- unique(c(topFeatures,intopfeat));
          if ((useFastCor) && (length(colused)>1000))
          {
            # Fast correlation for large matrices and linear fitting
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
          cormat[is.na(cormat)] <- 0;
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
#
        if ((nextthr == nextthra) && (addedlist == 0) && (mxScor > thr))
        {
          # loop faster if there were no changes
            lp <- lp + 2
        }
        if ((addedlist == 0) && (mxScor > thr))
        {
          # Add to check the next threshold value
          skv <- skv + 1 
        }
        nextthra <- sum(thrvalues > mxScor) + skv
      
      }
      
# The transformation has been completed. Next we do the final estimations 
      if (lp >= maxLoops)
      {
        warning("Maximum number of iterations reached: Failed to achieve target correlation\n")
      }
      # Do the final bookkeeping: Get the final transformation matrix and transform the input data 
      correlationMeasureEvolution <- c(correlationMeasureEvolution,max(cormat));
      sparcityEvolution <- c(sparcityEvolution,sum(DeCorrmatrix == 0));

      betamatrix <- NULL;
      tmparincluded <- varincluded;
      # Only variables involved in transformation
      varincluded <- tmparincluded[tmparincluded %in% totused];
      bfeat <- bfeat[bfeat %in% varincluded];
      if (useDeCorr)
      {
        # If estimating transformations let do the final transformation of the observed data

        DeCorrmatrix <- DeCorrmatrix[varincluded,varincluded] # Remove variables that were not used in the transformation
        abfeat <- varincluded[apply(DeCorrmatrix!=0,2,sum)==1]
        bfeat <- unique(c(bfeat,abfeat));
        dataTransformed <- data
        if (length(varincluded) > 1)
        {
          # Transform the observed data
          dataTransformed[,varincluded] <- Rfast::mat.mult(as.matrix(data[,varincluded]),DeCorrmatrix);
          
          
          # Observing internal behavior
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
            # Recompute correlation matrix only in verbose mode
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
              cormat[is.na(cormat)] <- 0;
              cat ("\n\r [",lp,"],",max(cormat),"Decor Dimension:",length(varincluded),"Nused:",length(totused),". Cor to Base:",length(correlatedToBase),", ABase:",length(AdrivingFeatures),", Outcome Base:",length(drivingFeatures),"\n\r")
          }
          
          # One last step: Update the column names to reflect that they are latent variables
          changed <- (apply(DeCorrmatrix!=0,2,sum) > 1);
          lavariables <- colnames(DeCorrmatrix)[changed]
          newnames <- colnames(dataTransformed)
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(dataTransformed) <- newnames
          newnames <- colnames(DeCorrmatrix)
          newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
          colnames(DeCorrmatrix) <- newnames

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
  # Store evolution data
  IDeAEvolution <-list(Corr=correlationMeasureEvolution,Spar=sparcityEvolution/(decordim^2));
  
  # The VarRatio stores the amount of explained variance
  VarRatio <- apply(dataTransformed[,colnames(DeCorrmatrix)],2,var);
  VarRatio <- VarRatio/apply(data[,rownames(DeCorrmatrix)],2,var);
  
  VarRatio <- VarRatio[order(-VarRatio)]
  noaltered <- colnames(data)[!(colnames(data) %in% rownames(DeCorrmatrix))]
  if (length(noaltered)>0)
  {
    varone <- rep(1,length(noaltered))
    names(varone) <- noaltered
    VarRatio <-  c(varone,VarRatio)
  }

# Return attributes
  attr(dataTransformed,"UPLTM") <- DeCorrmatrix;
  attr(dataTransformed,"fscore") <- fscore;
  attr(dataTransformed,"TotalAdjustments") <- countf;
  attr(dataTransformed,"drivingFeatures") <- AdrivingFeatures;
  attr(dataTransformed,"unaltered") <- bfeat;
  attr(dataTransformed,"LatentVariables") <- lavariables;
  attr(dataTransformed,"useDeCorr") <- useDeCorr;
  attr(dataTransformed,"unipvalue") <- adjunipvalue;
  attr(dataTransformed,"R.critical") <- rcrit;
  attr(dataTransformed,"IDeAEvolution") <- IDeAEvolution;
  attr(dataTransformed,"VarRatio") <- VarRatio;
  return(dataTransformed)
}

# transform new data
predictDecorrelate <- function(decorrelatedobject,testData)
{
  if (attr(decorrelatedobject,"useDeCorr") && !is.null(attr(decorrelatedobject,"UPLTM")))
  {
    decorMat <- attr(decorrelatedobject,"UPLTM")
    testData[,rownames(decorMat)] <- Rfast::mat.mult(as.matrix(testData[,rownames(decorMat)]),decorMat);
    newnames <- colnames(testData)
    lavariables <- attr(decorrelatedobject,"LatentVariables");
    newnames[newnames %in% lavariables] <- paste("La_",newnames[newnames %in% lavariables],sep="")
    colnames(testData) <- newnames
  }
  else
  {
    warning("The object does not have a decorrelation Matrix\n")
  }
  return (testData)
}

# Get the latent variable formulas as characters 
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

# Get the transformation coefficients for each latent variable
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

## Get the coefficients associated with each one of the observed features
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