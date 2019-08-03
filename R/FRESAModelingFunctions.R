CVsignature <- function(formula = formula, data=NULL, ...)
{
	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	if (length(usedFeatures)<5) #if less than 5 features, just use a glm fit.
	{
		warning("Less than five features. Returning a glm model");
		result <- glm(formula,data=data,na.action=na.exclude,family=binomial(link=logit));
	}
	else
	{
		parameters <- list(...);
		target="All";
		CVFolds=0;
		repeats=9;
		distanceFunction=signatureDistance;
		method="pearson";
		if (!is.null(parameters$target)) target=parameters$target;
		if (!is.null(parameters$CVFolds)) CVFolds=parameters$CVFolds;
		if (!is.null(parameters$repeats)) repeats=parameters$repeats;
		if (!is.null(parameters$distanceFunction)) distanceFunction=parameters$distanceFunction;
		if (!is.null(parameters$method)) method=parameters$method;
		cvsig <- getSignature(data=data,varlist=usedFeatures,Outcome=baseformula[2],target,CVFolds,repeats,distanceFunction,method);
		variable.importance <- 1:length(cvsig$featureList);
		names(variable.importance) <- cvsig$featureList;
		result <- list(fit=cvsig,method=method,variable.importance=variable.importance);
		class(result) <- "FRESAsignature";
	}
	return (result);
}

predict.FRESAsignature <- function(object, ...) 
{
  	parameters <- list(...);
  	testframe <- parameters[[1]];
	method <- object$method;
	if (!is.null(parameters$method)) method=parameters$method;
	controlDistances <- signatureDistance(object$fit$controlTemplate,testframe,method);
	caseDistances <- signatureDistance(object$fit$caseTamplate,testframe,method);
	distance <- controlDistances-caseDistances;
	return (distance);
}

KNN_method <- function(formula = formula, data=NULL, ...)
{
  	parameters <- list(...);
	if (is.null(parameters$kn))
	{
		kn <- as.integer(sqrt(nrow(data))+0.5);
	}
	else
	{
		kn <- parameters$kn;
	}

	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	tlabs <- attr(terms(formula,data=data),"term.labels");
	if (length(tlabs) > 0)
	{
		usedFeatures <- tlabs;
	}
#	print(usedFeatures);


	scaledData <- as.data.frame(data[,usedFeatures]);
	scaleMethod <- "None";
	if (is.null(parameters$scaleMethod))
	{
		scaleMethod <- "Norm";
	}
	else
	{
		scaleMethod  <- parameters$scaleMethod;
	}
	mean_vec <- NULL;
	disp_vec <- NULL;
	if (scaleMethod != "None") 
	{
		scaledParam <- FRESAScale(scaledData,scaledData,scaleMethod);
		mean_vec <- scaledParam$refMean;
		disp_vec <- scaledParam$refDisp;
		scaledData <- scaledParam$scaledData;
	}
	
	result <- list(trainData=as.data.frame(data[,usedFeatures]),scaledData=scaledData,classData=data[,baseformula[2]],outcome=baseformula[2],usedFeatures=usedFeatures,mean_col=mean_vec,disp_col=disp_vec,kn=kn,scaleMethod=scaleMethod);
	class(result) <- "FRESAKNN"
	return(result);
}


predict.FRESAKNN <- function(object, ...) 
{
  	parameters <- list(...);
  	testframe <- parameters[[1]];

	testframe <- as.data.frame(testframe[,object$usedFeatures]);
	trainframe <- object$scaledData
	
	if (object$scaleMethod != "None")
	{
		testframe <- FRESAScale(testframe,object$trainData,object$scaleMethod,object$mean_col,object$disp_col)$scaledData;
	}

	knnclass <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn,prob=TRUE))
	if (inherits(knnclass, "try-error")) knnclass <- numeric(nrow(testframe));
	if (length(table(knnclass))==2)
	{
		prop <- attributes(knnclass);
		knnclass <- abs(prop$prob-1*(knnclass=="0"))
	}
    return(knnclass);
}


LASSO_MIN <- function(formula = formula, data=NULL, ...)
{
if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 
	s <- "lambda.min";
	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	if (length(usedFeatures)<5) #if less than 5 features, just use a lm fit.
	{
		warning("Less than five features. Returning a lm model");
		result <- lm(formula,data);
	}
	else
	{
		if (sum(str_count(baseformula,"Surv")) > 0)
    {
      featuresOnSurvival <- baseformula[2]
      featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
      featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
      featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
      featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
      
      usedFeatures <- colnames(data)[!(colnames(data) %in%  featuresOnSurvivalObject[[1]])]
      x <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
      y <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
      baseformula <- gsub(featuresOnSurvivalObject[[1]][1],"x",baseformula)
      baseformula <- gsub(featuresOnSurvivalObject[[1]][2],"y",baseformula)
      result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),survival::Surv(x,y),family = "cox"),s = s,formula = formula,usedFeatures=usedFeatures);
      class(result) <- "FRESA_LASSO"
      
      coef(result$fit,s)
    }
    else
    {
		result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]),...),s = s,formula = formula,outcome = baseformula[2],usedFeatures = usedFeatures)
		class(result) <- "FRESA_LASSO"
		}
	}
  return(result);
}

LASSO_1SE <- function(formula = formula, data=NULL, ...)
{
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    install.packages("glmnet", dependencies = TRUE)
  } 
  s <- "lambda.1se";
  baseformula <- as.character(formula);
  usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
  if (length(usedFeatures)<5) #if less than 5 features, just use a lm fit.
  {
    warning("Less than five features. Returning a lm model");
    result <- lm(formula,data);
  }
  else
  {
    if (sum(str_count(baseformula,"Surv")) > 0)
    {
      featuresOnSurvival <- baseformula[2]
      featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
      featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
      featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
      featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
      
      usedFeatures <- colnames(data)[!(colnames(data) %in%  featuresOnSurvivalObject[[1]])]
      x <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
      y <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
      baseformula <- gsub(featuresOnSurvivalObject[[1]][1],"x",baseformula)
      baseformula <- gsub(featuresOnSurvivalObject[[1]][2],"y",baseformula)
      result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),survival::Surv(x,y),family = "cox"),s = s,formula = formula,usedFeatures=usedFeatures);
      class(result) <- "FRESA_LASSO"
      
      coef(result$fit,s)
    }
    else
    {
      result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]),...),s = s,formula = formula,outcome = baseformula[2],usedFeatures = usedFeatures)
      class(result) <- "FRESA_LASSO"
    }
  }
  return(result);
}

predict.FRESA_LASSO <- function(object,...) 
{
  	parameters <- list(...);
  	testData <- parameters[[1]];
    pLS <- predict(object$fit,as.matrix(testData[,object$usedFeatures]), s = object$s);
    return(pLS);
}

BESS <- function(formula = formula, data=NULL, ...)
{
  if (!requireNamespace("BeSS", quietly = TRUE)) {
    install.packages("BeSS", dependencies = TRUE)
  } 
  
  baseformula <- as.character(formula);
  featuresOnSurvival <- baseformula[2]
  featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
  featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
  
  usedFeatures <- colnames(data)[!(colnames(data) %in%  featuresOnSurvivalObject[[1]])]
  x <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
  y <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
  baseformula <- gsub(featuresOnSurvivalObject[[1]][1],"x",baseformula)
  baseformula <- gsub(featuresOnSurvivalObject[[1]][2],"y",baseformula)
  result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]), survival::Surv(x, y), s.min=1, family = "cox"),formula = formula,usedFeatures=usedFeatures);
  class(result) <- "bess"
  
  return(result);
}


TUNED_SVM <- function(formula = formula, data=NULL,...)
{
  if (!requireNamespace("e1071", quietly = TRUE)) {
	  install.packages("e1071", dependencies = TRUE)
	  }
	obj <- e1071::tune.svm(formula, data=data,gamma = 2^(2*(-10:0)), cost = 2^(2*(-5:2)));
	fit <- e1071::svm(formula, data=data,gamma=obj$best.parameters$gamma,cost=obj$best.parameters$cost,...);
	result <- list(fit = fit,tuneSVM=obj);
	class(result) <- "FRESA_SVM"
	return(result);
}

predict.FRESA_SVM <- function(object,...) 
{
  	parameters <- list(...);
  	testData <- parameters[[1]];
    pLS <- predict(object$fit,...);
    return(pLS);
}

NAIVE_BAYES <- function(formula = formula, data=NULL, ...)
{
if (!requireNamespace("naivebayes", quietly = TRUE)) {
   install.packages("naivebayes", dependencies = TRUE)
} 
	baseformula <- as.character(formula);
	if (class(data[,baseformula[2]]) != "factor") data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	result <- list(fit = naivebayes::naive_bayes(formula,data,...))
	class(result) <- "FRESA_NAIVEBAYES"
	return(result);
}


predict.FRESA_NAIVEBAYES <- function(object,...) 
{
  	parameters <- list(...);
  	testData <- parameters[[1]];
    pLS <- as.numeric(as.character(predict(object$fit,testData)));
	if (length(table(pLS)) == 2)
	{
		prop <- predict(object$fit,testData,type = "prob");
		pLS <- prop[,2];
		pLS[is.nan(pLS)] <- 0.5;
		pLS[is.na(pLS)] <- 0.5;
	}
    return(pLS);
}

LM_RIDGE_MIN <- function(formula = formula, data=NULL, ...)
{
	if (!requireNamespace("MASS", quietly = TRUE)) {
	   install.packages("MASS", dependencies = TRUE)
	}
	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	if (length(usedFeatures)<5) #if less than 5 features, just use a lm fit.
	{
		warning("Less than five features. Returning a lm model");
		fit <- lm(formula,data);
	}
	else
	{
		parameters <- list(...);
		if (is.null(parameters$lambda))
		{
			lambda = seq(0,0.2,0.002);
			fit <- MASS::lm.ridge(formula,data,lambda=lambda,...);
			fit$coef <- fit$coef[,which.min(fit$GCV)];
		}
		else
		{
			fit <- MASS::lm.ridge(formula,data,...);
			if (length(parameters$lambda)>1)
			{
				fit$coef <- fit$coef[,which.min(fit$GCV)];
			}
		}
		class(fit) <- c("FRESA_RIDGE",class(fit))
	}
	return(fit);
}	

predict.FRESA_RIDGE <- function(object,...)
{
 # Predict MASS:lm.ridge is not implemented so I added to FRESA.CAD
  	parameters <- list(...);
  	testData <- parameters[[1]];
	ridgenames <- names(object$xm);
	pr = scale(as.matrix(testData[,ridgenames]),center =  object$xm, scale = object$scales) %*% object$coef + object$ym;
	return(pr)
}


BOOST_BSWiMS <- function(formula = formula, data=NULL, ...)
{
	if (class(formula)=="character")
	{
		formula <- formula(formula);
	}
	else
	{
		baseformula <- as.character(formula);
		baseformula[3] <- str_replace_all(baseformula[3],"[.]","1");
		baseformula <- paste(baseformula[2],"~",baseformula[3]);
		formula <- formula(baseformula);
	}
	varlist <- attr(terms(formula),"variables")
	dependent <- as.character(varlist[[2]])
	Outcome = dependent[1];
	if (length(dependent)==3)
	{
		Outcome = dependent[3];
	}
	outcomedata <- data[,Outcome];
	outcomeTable <- table(outcomedata);

	thr <- 0.5
	thr2 <- 1.0 - thr
	orgModel <- BSWiMS.model(formula,data,...);
	orgPredict <- predict(orgModel,data)
	incorrectSet <- ((orgPredict >= thr) & outcomedata == 0) | ((orgPredict < thr2) & outcomedata == 1);
	alternativeModel <- NULL;
	classsModel <- NULL;
	if (sum(incorrectSet) > 20)
	{
		if (min(table(data[incorrectSet,Outcome])) > 10)
		{
			alternativeModel <- BSWiMS.model(formula,data[incorrectSet,])
			classData <- trainSet[,!(colnames(data) %in% Outcome) ]
			classData$label <- 1*incorrectSet;
			classsModel <- BSWiMS.model(label~1,classData);
		}
	}
	result <- list(original = orgModel,alternativeModel = alternativeModel,classModel = classsModel )
	class(result) <- "FRESA_BOOST"
	return(result);
}


predict.FRESA_BOOST <- function(object,...) 
{
  	parameters <- list(...);
  	testData <- parameters[[1]];
	thr <- 0.55;
	if (length(parameters) > 1)
	{
		thr <- parameters[[2]];
	}
    pLS <- predict(object$original,testData);
	if (!is.null(object$alternativeModel))
	{
		palt <- predict(object$alternativeModel,testData);
		classPred <- predict(object$classsModel,testData);
		pLS[classPred > thr] <- palt[classPred > thr];
	}
    return(pLS);
}

