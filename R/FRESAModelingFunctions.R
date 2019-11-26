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

	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	tlabs <- attr(terms(formula,data=data),"term.labels");
	if (length(tlabs) > 0)
	{
		usedFeatures <- tlabs;
	}
	if (is.null(parameters$kn))
	{
		tb <- table(data[,baseformula[2]]);
		kn <- as.integer(sqrt(min(tb))+0.5);
	}
	else
	{
		kn <- parameters$kn;
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
		scaleMethod	<- parameters$scaleMethod;
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
	result$selectedfeatures <- usedFeatures;
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
		knnclass <- abs(prop$prob-1*(knnclass == "0"))
		if (object$kn > 2)
		{
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn-1,prob=TRUE))
			prop_1 <- attributes(knnclass_1);
			knnclass_2 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn+1,prob=TRUE))
			prop_2 <- attributes(knnclass_2);
			knnclass <- (knnclass + 0.5*abs(prop_1$prob-1*(knnclass_1 == "0")) + 0.5*abs(prop_2$prob-1*(knnclass_2 == "0")))/2.0;
		}
	}
	else
	{
		prop <- attributes(knnclass)$prob;
		if (object$kn > 2)
		{
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn-1,prob=TRUE))
			prop_1 <- attributes(knnclass_1)$prob;
			testclass <- as.character(knnclass_1) == as.character(knnclass)
			prop[testclass] <- 0.75*prop[testclass] + 0.25*prop_1[testclass];
			prop[!testclass] <- 0.75*prop[!testclass] + 0.25*(1.0-prop_1[!testclass]);
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn+2,prob=TRUE))
			prop_1 <- attributes(knnclass_1)$prob;
			testclass <- as.character(knnclass_1) == as.character(knnclass)
			prop[testclass] <- 0.75*prop[testclass] + 0.25*prop_1[testclass];
			prop[!testclass] <- 0.75*prop[!testclass] + 0.25*(1.0-prop_1[!testclass]);
			attr(knnclass,"prob") <- prop;
		}
	}
	return(knnclass);
}


GLMNET <- function(formula = formula, data=NULL,coef.thr=0.001,s="lambda.min",...)
{
	if (!requireNamespace("glmnet", quietly = TRUE)) 
	{
		install.packages("glmnet", dependencies = TRUE)
	} 
	parameters <- list(...);
	isSurv <- FALSE;
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
			isSurv <- TRUE;
			featuresOnSurvival <- baseformula[2]
			featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
			featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
			featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
			featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
			
			usedFeatures <- colnames(data)[!(colnames(data) %in%	featuresOnSurvivalObject[[1]])]
			x <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
			y <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
			baseformula <- gsub(featuresOnSurvivalObject[[1]][1],"x",baseformula)
			baseformula <- gsub(featuresOnSurvivalObject[[1]][2],"y",baseformula)
			result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),survival::Surv(x,y),family = "cox",...),s=s,formula = formula,usedFeatures=usedFeatures);
		}
		else
		{
			result <- list(fit = glmnet::cv.glmnet(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]),...),s=s,formula = formula,outcome = baseformula[2],usedFeatures = usedFeatures)
		}
	}
	coefthr <- numeric(1*(!isSurv)+length(usedFeatures));
	if(!is.null(parameters$alpha))
	{
		if(parameters$alpha < 1)
		{
			coefthr <- apply(data[,usedFeatures],2,sd, na.rm = TRUE);
			coefthr <- coef.thr/coefthr;
			if (!isSurv)
			{
				coefthr <- c(0,coefthr);
			}
		}
	}

	cf <- coef(result$fit,s);
	selectedFeatures <- character();
	lcoef <- numeric();
	if (class(cf) == "list")
	{
		lcoef <- list();
		for (cl in 1:length(cf))
		{
			cenet <- as.matrix(cf[[cl]]);
			if (!is.null(cenet))
			{
				lft <- cenet[as.vector(abs(cenet[,1]) > coefthr),,drop=FALSE];
				lcoef[[cl]] <- as.numeric(lft);
				names(lcoef[[cl]]) <- rownames(lft);
				sF <- rownames(lft);
				selectedFeatures <- sF[rownames(lft) %in% usedFeatures];
			}
		}			
	}
	else
	{
		cenet <- as.matrix(cf);
		if (!is.null(cenet))
		{
			lft <- cenet[as.vector(abs(cenet[,1]) > coefthr),,drop=FALSE];
			lcoef <- as.numeric(lft);
			names(lcoef) <- rownames(lft);
			sF <- rownames(lft);
			selectedFeatures <- sF[rownames(lft) %in% usedFeatures];
		}
	}
	result$selectedfeatures <- unique(selectedFeatures);
	result$coef <- lcoef;
	class(result) <- "FRESA_GLMNET"
	return(result);
}

LASSO_MIN <- function(formula = formula, data=NULL, ...)
{
	result <- GLMNET(formula,data,s = "lambda.min",...);
	return (result);
}


LASSO_1SE <- function(formula = formula, data=NULL,...)
{
	result <- GLMNET(formula,data,s = "lambda.1se",...);
	return (result);
}

GLMNET_RIDGE_MIN <- function(formula = formula, data=NULL, ...)
{
	result <- GLMNET(formula,data,s = "lambda.min", alpha=0,...);
	return (result);
}

GLMNET_ELASTICNET_MIN <- function(formula = formula, data=NULL, ...)
{
	result <- GLMNET(formula,data,s = "lambda.min",alpha=0-95,...);
	return (result);
}

GLMNET_RIDGE_1SE <- function(formula = formula, data=NULL, ...)
{
	result <- GLMNET(formula,data,s = "lambda.1se", alpha=0,...);
	return (result);
}

GLMNET_ELASTICNET_1SE <- function(formula = formula, data=NULL, ...)
{
	result <- GLMNET(formula,data,s = "lambda.1se",alpha=0-95,...);
	return (result);
}

predict.FRESA_GLMNET <- function(object,...) 
{
		parameters <- list(...);
		testData <- parameters[[1]];
		pLS <- predict(object$fit,as.matrix(testData[,object$usedFeatures]), s = object$s);
		return(pLS);
}

BESS <- function(formula = formula, data=NULL, method="sequential", ic.type="BIC",...)
{
	if (!requireNamespace("BeSS", quietly = TRUE)) {
		install.packages("BeSS", dependencies = TRUE)
	} 
	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	
	if (sum(str_count(baseformula,"Surv")) > 0)
	{
		baseformula <- as.character(formula);
		featuresOnSurvival <- baseformula[2]
		featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
		featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
		featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
		featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
		
		usedFeatures <- colnames(data)[!(colnames(data) %in%	featuresOnSurvivalObject[[1]])]
		x <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
		y <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
		baseformula <- gsub(featuresOnSurvivalObject[[1]][1],"x",baseformula)
		baseformula <- gsub(featuresOnSurvivalObject[[1]][2],"y",baseformula)
		result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]), survival::Surv(x, y), method=method, family = "cox",ic.type=ic.type,...),formula = formula,usedFeatures=usedFeatures);
		bessCoefficients <- result$fit$bestmodel$coefficients
	}
	else
	{
		tb <- table(data[,baseformula[2]]);
		if (length(tb)>2)
		{
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "gaussian", epsilon = 1e-12,ic.type=ic.type,...),formula = formula,usedFeatures=usedFeatures);
		}
		else
		{
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "binomial", epsilon = 0,ic.type=ic.type,...), formula = formula,usedFeatures=usedFeatures);

		}
		bessCoefficients <- result$fit$bestmodel$coefficients[-1];
	}
	if (!is.null(bessCoefficients))
	{
		result$selectedfeatures <- gsub("xbest","",names(bessCoefficients));
	}
	result$ic.type <- ic.type;
	class(result) <- "FRESA_BESS"
	
	return(result);
}

BESS_GSECTION <- function(formula = formula, data=NULL, method="gsection", ic.type="NULL",...)
{
	result <- BESS(formula, data, method, ic.type,...);
	return(result);
}

BESS_GIC <- function(formula = formula, data=NULL, ic.type="GIC",...)
{
	result <- BESS(formula = formula, data = data, ic.type=ic.type,...);
	return(result);
}

predict.FRESA_BESS <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pLS <- NULL;
	if (!is.null(parameters$type))
	{
		type <- parameters$type;
		if (type == "response")
		{
			newdata <- as.data.frame(cbind(y=1:nrow(testData),testData[,object$selectedfeatures]))
			colnames(newdata) <- c("y",paste("xbest",object$selectedfeatures,sep=""))
			object$fit$bestmodel$formula <- formula(paste("y~",paste(colnames(newdata)[-1],collapse = " + ")))
			object$fit$bestmodel$terms <- terms(object$fit$bestmodel$formula)
			pLS <- predict(object$fit$bestmodel,newdata,type=type);
		}
	}
	if (is.null(pLS))
	{
		if (object$fit$method == "gsection")
		{
			pLS <- predict(object$fit,testData,type="opt");
		}
		else
		{
			type = object$ic.type;
			if (!is.null(parameters$type))
			{
				type <- parameters$type;
			}
			pLS <- predict(object$fit,testData,type=type);
		}
	}
	return(pLS);
}

TUNED_SVM <- function(formula = formula, data=NULL,gamma = 10^(-5:-1), cost = 10^(-3:1),...)
{
	if (!requireNamespace("e1071", quietly = TRUE)) {
		install.packages("e1071", dependencies = TRUE)
		}
	obj <- e1071::tune.svm(formula, data=data,gamma = gamma, cost = cost);
	fit <- e1071::svm(formula, data=data,gamma=obj$best.parameters$gamma,cost=obj$best.parameters$cost,...);
	
	parameters <- list(...);
	probability <- NULL;
	if 	(!is.null(parameters$probability))
	{
		probability <- parameters$probability
	}

	result <- list(fit = fit,tuneSVM=obj,probability = probability);
		
	class(result) <- "FRESA_SVM"
	return(result);
}

predict.FRESA_SVM <- function(object,...) 
{
		parameters <- list(...);
		testData <- parameters[[1]];
		if (is.null(object$probability))
		{
			pLS <- predict(object$fit,...);
		}
		else
		{
			pLS <- predict(object$fit,testData,probability = object$probability);
			pLS <- attr(pLS,"probabilities")[,"1"];
		}
		return(pLS);
}

NAIVE_BAYES <- function(formula = formula, data=NULL,pca=TRUE,...)
{
if (!requireNamespace("naivebayes", quietly = TRUE)) {
	 install.packages("naivebayes", dependencies = TRUE)
} 
	baseformula <- as.character(formula);
	if (class(data[,baseformula[2]]) != "factor") data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	pcaobj <- NULL;
	scaleparm <- NULL;
	if (pca && (nrow(data) > 2*ncol(data)))
	{
		outcome <- data[,baseformula[2]];
		data <- data[,!(colnames(data) %in% baseformula[2])];
		scaleparm <- FRESAScale(data,method="Order");
		pcaobj <- prcomp(scaleparm$scaledData);
		data <- as.data.frame(cbind(as.numeric(as.character(outcome)),pcaobj$x));
		colnames(data) <- c(baseformula[2],colnames(pcaobj$x));
		data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	}
	result <- list(fit = naivebayes::naive_bayes(formula,data,...),pcaobj=pcaobj,outcome=baseformula[2],scaleparm=scaleparm);
	class(result) <- "FRESA_NAIVEBAYES"
	return(result);
}


predict.FRESA_NAIVEBAYES <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	if (!is.null(object$pcaobj))
	{
		testData <- FRESAScale(testData,refMean=object$scaleparm$refMean,refDisp=object$scaleparm$refDisp)$scaledData;
		testData <- predict(object$pcaobj,testData);
	}
	else
	{
		testData <- testData[,!(colnames(testData) %in% object$outcome)];
	}
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
	pr = scale(as.matrix(testData[,ridgenames]),center =	object$xm, scale = object$scales) %*% object$coef + object$ym;
	return(pr)
}


HCLAS_CLUSTER <- function(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.0,classMethod=KNN_method,classModel.Control=NULL,minsize=10,...)
{
	if (class(formula) == "character")
	{
		formula <- formula(formula);
	}
	varlist <- attr(terms(formula,data=data),"variables")
	dependent <- as.character(varlist[[2]])
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[3];
	}


	alternativeModel <- list();
	correctSet <- list();
	classModel <- NULL;
	selectedfeatures <- colnames(data)[!(colnames(data) %in% Outcome)]
	orgModel <- method(formula,data,...);
	if (!is.null(orgModel$selectedfeatures))
	{
		selectedfeatures <- orgModel$selectedfeatures;
	}
	else
	{
		if (!is.null(orgModel$bagging))
		{
			selectedfeatures <- names(orgModel$bagging$frequencyTable);
		}
	}
	accuracy <- 1.0;
	inserted <- TRUE;
	n=0;
	classData <- data;
	classData[,Outcome] <- rep(0,nrow(classData));
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(orgModel,data);
		if ((min(thePredict) < -0.1) || (max(thePredict) > 1.1 ))
		{
			thePredict <- 1.0/(1.0+exp(-thePredict));
		}
		outcomedata <- data[,Outcome];
		correct <- ((thePredict >= 0.5) == (outcomedata > 0));
		accuracy <- sum(correct)/nrow(data);
		if (sum(!correct) > minsize)
		{
			nextdata <- data;
			while (inserted)
			{
				inserted <- FALSE;
				preData <- nextdata;
				outcomedata <- preData[,Outcome];
				falseP <- (thePredict > (0.5 - hysteresis)) & (outcomedata == 0);
				falseN <- (thePredict < (0.5 + hysteresis)) & (outcomedata == 1);
				if ( sum(falseP | falseN) > (nrow(preData)/2) )
				{
					falseP <- (thePredict > 0.5) & (outcomedata == 0);
					falseN <- (thePredict < 0.5) & (outcomedata == 1);
				}
				if ((sum(falseP) >= (minsize/2)) && (sum(falseN) >= (minsize/2)))
				{
					incorrectSet <- falseP | falseN;
					nextdata <- preData[incorrectSet,];
					alternativeM <- method(formula,nextdata,...);
					nselected <- character();
					if (!is.null(alternativeM$selectedfeatures))
					{
						nselected <- alternativeM$selectedfeatures;
					}
					else
					{
						if (!is.null(alternativeM$bagging))
						{
							nselected <- names(alternativeM$bagging$frequencyTable);
						}
					}
					if (length(nselected)>0)
					{
						inserted <- TRUE;
						n <- n + 1;
						selectedfeatures <- c(selectedfeatures,nselected);
						selectedfeatures <- unique(selectedfeatures);
						thePredict <- rpredict(alternativeM,nextdata);
						if ((min(thePredict) < -0.1) || (max(thePredict) > 1.1))
						{
							thePredict <- 1.0/(1.0+exp(-thePredict));
						}
						correctSet[[n]] <- rownames(preData[!incorrectSet,]);
						cat("[",sum(incorrectSet),"]")
						alternativeModel[[n]] <- alternativeM;
					}
				}
				else
				{
					if ( ((sum(falseP) >= minsize) || (sum(falseN) >= minsize)) )
					{
						incorrectSet <- falseP | falseN;
						inserted <- FALSE;
						n <- n + 1;
						correctSet[[n]] <- rownames(preData[!incorrectSet,]);
						cat("(",sum(incorrectSet),")")
						alternativeModel[[n]] <- sum(falseN)/(sum(falseP)+sum(falseN));
					}
				}
			}
			if (n>0)
			{
				allCorrectSet <- correctSet[[1]];
				classData[,Outcome] <- rep(n,nrow(classData));
				for (i in 1:n)
				{
					classData[correctSet[[i]],Outcome] <- i - 1;
					allCorrectSet <- c(allCorrectSet,correctSet[[i]]);
				}
				if (is.null(classModel.Control))
				{
					classModel <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
				}
				else
				{
					classData[,Outcome] <- as.factor(classData[,Outcome]);
					classModel <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
				}
			}
		}
	}
	result <- list(original = orgModel,
					alternativeModel = alternativeModel,
					classModel = classModel,
					accuracy=accuracy,
					selectedfeatures = selectedfeatures,
					hysteresis=hysteresis,
					classSet=classData[,Outcome]
					)
	class(result) <- "FRESA_HCLAS"
	return(result);
}

HCLAS_EM_CLUSTER <- function(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.0,classMethod=KNN_method,classModel.Control=NULL,minsize=10,...)
{
	if (class(formula) == "character")
	{
		formula <- formula(formula);
	}
	varlist <- attr(terms(formula,data=data),"variables")
	dependent <- as.character(varlist[[2]])
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[3];
	}


	alternativeModel <- list();
	correctSet <- NULL;
	classModel <- NULL;
	selectedfeatures <- colnames(data)[!(colnames(data) %in% Outcome)]
	baseModel <- method(formula,data,...);
	orgModel <- baseModel;
	if (!is.null(baseModel$selectedfeatures))
	{
		selectedfeatures <- baseModel$selectedfeatures;
	}
	else
	{
		if (!is.null(baseModel$bagging))
		{
			selectedfeatures <- names(baseModel$bagging$frequencyTable);
		}
	}
	accuracy <- 1.0;
	changes <- 1;
	n=0;
	classData <- data;
	classData[,Outcome] <- rep(0,nrow(classData));
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(baseModel,data);
		if ((min(thePredict) < -0.1) || (max(thePredict) > 1.1 ))
		{
			thePredict <- 1.0/(1.0+exp(-thePredict));
		}
		outcomedata <- data[,Outcome];
		correct <- ((thePredict >= 0.5) == (outcomedata > 0));
		accuracy <- sum(correct)/nrow(data);
		if (sum(!correct) > minsize)
		{
			outcomedata <- data[,Outcome];
			falseP <- (thePredict > (0.5 - hysteresis)) & (outcomedata == 0);
			falseN <- (thePredict < (0.5 + hysteresis)) & (outcomedata == 1);
			if ( sum(falseP | falseN) > (nrow(data)/2) )
			{
				falseP <- (thePredict > 0.5) & (outcomedata == 0);
				falseN <- (thePredict < 0.5) & (outcomedata == 1);
			}
			if ((sum(falseP) >= (minsize/2)) && (sum(falseN) >= (minsize/2)))
			{
				secondSet <- falseP | falseN;
				firstSet <- !secondSet;
				loops <- 0;
				firstModel <- NULL;
				secondModel <- NULL;
				while ((changes > 0) && (loops < 10))
				{
					loops <- loops + 1;
					n <- 0;
					changes <- 0;
					firstdata <- data[firstSet,];
					seconddata <- data[secondSet,];
					if ( (min(table(firstdata[,Outcome])) > (minsize/2)) && (min(table(seconddata[,Outcome])) > (minsize/2)) )
					{
						firstModel <- method(formula,firstdata,...);
						secondModel <- method(formula,seconddata,...);
						nselected <- character();
						if (!is.null(secondModel$selectedfeatures))
						{
							nselected <- secondModel$selectedfeatures;
						}
						else
						{
							if (!is.null(secondModel$bagging))
							{
								nselected <- names(secondModel$bagging$frequencyTable);
							}
						}
						if (length(nselected)>0)
						{
							n <- 1;
							firstPredict <- rpredict(firstModel,data);
							if ((min(firstPredict) < -0.1) || (max(firstPredict) > 1.1))
							{
								firstPredict <- 1.0/(1.0+exp(-firstPredict));
							}
							secondPredict <- rpredict(secondModel,data);
							if ((min(secondPredict) < -0.1) || (max(secondPredict) > 1.1))
							{
								secondPredict <- 1.0/(1.0+exp(-secondPredict));
							}
							d1 <-  abs(firstPredict-outcomedata);
							d2 <-  abs(secondPredict-outcomedata);
							nfirstSet <- (d1 < (d2 + hysteresis));
							changes <- sum(nfirstSet != firstSet);
							firstSet <- nfirstSet;
							secondSet <- (d2 < (d1 + hysteresis));
						}
					}
					cat("(",changes,")");
				}
				cat("[",sum(secondSet),"]")
				if (n > 0)
				{
					d1 <-  abs(firstPredict-outcomedata);
					d2 <-  abs(secondPredict-outcomedata);
					firstSet <- (d1 < d2);
				}
				correctSet <- rownames(data[firstSet,]);
				orgModel <- firstModel;
				alternativeModel[[1]] <- secondModel;
			}
			if (n>0)
			{
				nselected <- character();
				if (!is.null(firstModel$selectedfeatures))
				{
					nselected <- firstModel$selectedfeatures;
				}
				else
				{
					if (!is.null(firstModel$bagging))
					{
						nselected <- names(firstModel$bagging$frequencyTable);
					}
				}
				selectedfeatures <- c(selectedfeatures,nselected);
				selectedfeatures <- unique(selectedfeatures);
				nselected <- character();
				if (!is.null(secondModel$selectedfeatures))
				{
					nselected <- secondModel$selectedfeatures;
				}
				else
				{
					if (!is.null(secondModel$bagging))
					{
						nselected <- names(secondModel$bagging$frequencyTable);
					}
				}
				selectedfeatures <- c(selectedfeatures,nselected);
				selectedfeatures <- unique(selectedfeatures);
				classData[,Outcome] <- rep(1,nrow(classData));
				classData[correctSet,Outcome] <- 0;
				if (is.null(classModel.Control))
				{
					classModel <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
				}
				else
				{
					classData[,Outcome] <- as.factor(classData[,Outcome]);
					classModel <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
				}
			}
		}
	}
	result <- list(original = orgModel,
					alternativeModel = alternativeModel,
					classModel = classModel,
					accuracy=accuracy,
					selectedfeatures = selectedfeatures,
					hysteresis=hysteresis,
					classSet=classData[,Outcome],
					baseModel = baseModel
					)
	class(result) <- "FRESA_HCLAS"
	return(result);
}

predict.FRESA_HCLAS <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pLS <- rpredict(object$original,testData);
	if ((min(pLS) < -0.1) || (max(pLS) > 1.1))
	{
		pLS <- 1.0/(1.0+exp(-pLS));
	}
	if (!is.null(object$classModel))
	{
		tb <- table(object$classSet)/length(object$classSet);
		if (class(object$classModel)[[1]] == "FRESAKNN")
		{
			classPred <- predict(object$classModel,testData);
			if (length(object$alternativeModel) == 1)
			{
				classPred <- 1.0 - classPred;
				palt <- rpredict(object$alternativeModel[[1]],testData);
				if ((min(palt) < -0.1) || (max(palt) > 1.1))
				{
					palt <- 1.0/(1.0 + exp(-palt));
				}
				pLS <- classPred*pLS + (1.0 - classPred)*(palt + tb[2]*(1.0 - pLS))/(1.0 + tb[2]);
			}
			else
			{
				pmodel <- pLS;
				prbclas <- attributes(classPred)$prob;
				classPred <- as.numeric(as.character(classPred)) + 1;
				nm <- length(object$alternativeModel)
				for (n in 1:nm)
				{
					ptmp <- rpredict(object$alternativeModel[[n]],testData)
					if ((min(ptmp) < -0.1) || (max(ptmp) > 1.1))
					{
						ptmp <- 1.0/(1.0 + exp(-ptmp));
					}
					pmodel <- cbind(pmodel,ptmp);
				}
				itotclas <- 1.0/nm; 
				for (i in 1:length(pLS))
				{
					wts <- prbclas[i];
					pLS[i] <- pmodel[i,classPred[i]]*wts;
					for (n in 1:(nm + 1))
					{
						if (n != classPred[i])
						{
							wt <- (1.0 - prbclas[i])*tb[n]*itotclas;
							pLS[i] <- pLS[i] + pmodel[i,n]*wt;
							wts <- wts + wt;
						}
					}
					if (is.null(object$baseModel))
					{
						if (classPred[i] > 1)
						{
							for (n in 1:(classPred[i] - 1))
							{
								wt <- prbclas[i]*tb[n + 1];
								pLS[i] <- pLS[i] + wt*(1.0 - pmodel[i,n]);
								wts <- wts + wt;
							}
						}
					}
					pLS[i] <- pLS[i]/wts;
				}
			}
		}
		else
		{
			if (class(object$classModel)[[1]] == "svm.formula")
			{
				classPred <- predict(object$classModel,testData,probability = TRUE);
				pclase <- attributes(classPred)$probabilities;
				pmodel <- pLS;
				nm <- length(object$alternativeModel)
				for (n in 1:nm)
				{
					ptmp <- rpredict(object$alternativeModel[[n]],testData)
					if ((min(ptmp) < -0.1) || (max(ptmp) > 1.1))
					{
						ptmp <- 1.0/(1.0 + exp(-ptmp));
					}
					pmodel <- cbind(pmodel,ptmp);
				}
				classPred <- as.numeric(as.character(classPred));
				for (i in 1:length(pLS))
				{
					wm <- classPred[i];
					wts <- 0;
					pLS[i] <- 0;
					for (n in 0:nm)
					{	
						wt <-  pclase[i,as.character(n)];
						pLS[i] <- pLS[i]+wt*pmodel[i,(n + 1)];
						wts <- wts + wt;
					}
					if (is.null(object$baseModel))
					{
						if (wm > 0)
						{
							for (n in 0:(wm - 1))
							{	
								wt <-  pclase[i,wm+1]*tb[n + 2];
								pLS[i] <- pLS[i]+wt*(1.0 - pmodel[i,(n + 1)]);
								wts <- wts + wt;
							}
						}
					}
					pLS[i] <- pLS[i]/wts;
				}
			}
		}
	}
	return(pLS);
}


filteredFit <- function(formula = formula, data=NULL, filtermethod=univariate_Wilcoxon, classmethod=e1071::svm,filtermethod.control=list(pvalue=0.05,limit=0.1),...)
{
	if (class(formula) == "character")
	{
		formula <- formula(formula);
	}
	varlist <- attr(terms(formula,data=data),"variables")
	dependent <- as.character(varlist[[2]])
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[3];
	}
	fm <- NULL
	
	if (is.null(filtermethod.control))
	{
		fm <- filtermethod(data,Outcome);
	}
	else
	{
		fm <- do.call(filtermethod,c(list(data,Outcome),filtermethod.control));
	}
	usedFeatures <-  c(Outcome,names(fm));
	fit <- classmethod(formula,data[,usedFeatures],...);
	parameters <- list(...);
	result <- list(fit=fit,filter=fm,selectedfeatures = names(fm),usedFeatures = usedFeatures,parameters=parameters,asFactor=(class(data[,Outcome])=="factor"),classLen=length(table(data[,Outcome])));
	class(result) <- "FRESA_FILTERFIT";
	return (result)	
}

predict.FRESA_FILTERFIT <- function(object,...)
{
	parameters <- list(...);
	testData <- parameters[[1]];
	probability <- FALSE;
	if (!is.null(object$parameters$probability))
	{
		probability <- object$parameters$probability;
	}
	pLS <- rpredict(object$fit,testData[,object$usedFeatures],asFactor=object$asFactor,classLen=object$classLen,probability=probability,...);
	return (pLS);
}

ClustClass <- function(formula = formula, data=NULL, filtermethod=univariate_Wilcoxon, clustermethod=NULL, classmethod=LASSO_1SE,filtermethod.control=list(pvalue=0.05,limit=0.1),clustermethod.control=NULL,classmethod.control=list(family = "binomial"))
{
	if (class(formula) == "character")
	{
		formula <- formula(formula);
	}
	varlist <- attr(terms(formula,data=data),"variables")
	dependent <- as.character(varlist[[2]])
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[3];
	}

	outcomedata <- data[,Outcome];
	totsamples <- nrow(data);
	minSamples <- max(5,0.05*totsamples);
	clus <- NULL
	fm <- NULL
	
	if (is.null(filtermethod.control))
	{
		fm <- filtermethod(data,Outcome);
	}
	else
	{
		fm <- do.call(filtermethod,c(list(data,Outcome),filtermethod.control));
	}
	if (is.null(clustermethod.control))
	{
		clus <- clustermethod(data[,names(fm)]);
	}
	else
	{
		clus <- do.call(clustermethod,c(list(data[,names(fm)]),clustermethod.control));
	}
	selectedfeatures <- names(fm);
	tb <- table(clus$classification);
	classlabels <- as.numeric(names(tb));
	models <- list();
	if (length(classlabels) > 1)
	{
			tb <- table(clus$classification,outcomedata);
			for (i	in 1:nrow(tb))
			{
				if (min(tb[i,]) > minSamples)
				{
					if (is.null(classmethod.control))
					{
						models[[i]] <- classmethod(formula,subset(data,clus$classification == classlabels[i]));
					}
					else
					{
						models[[i]] <- do.call(classmethod,c(list(formula,subset(data,clus$classification == classlabels[i])),classmethod.control));
					}
				}
				else
				{
					models[[i]] <- as.numeric(colnames(tb)[which.max(tb[i,])]);
				}
			}
	}
	else
	{
		if (is.null(classmethod.control))
		{
			models[[1]] <- classmethod(formula,data);
		}
		else
		{
			models[[1]] <- do.call(classmethod,c(list(formula,data),classmethod.control));
		}
	}
	result <- list(features = fm,cluster = clus,models = models);
	result$selectedfeatures <- selectedfeatures;

	class(result) <- "CLUSTER_CLASS"
	return(result);
}

predict.CLUSTER_CLASS <- function(object,...)
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pLS <- predict(object$cluster,testData[,names(object$features)])$classification;
	tb <- table(pLS);
	index <- as.numeric(names(tb));
	for (i in 1:nrow(tb))
	{
		predeictset <- (pLS == index[i]);
		if (class(object$models[[index[i]]]) == "numeric")
		{
			pLS[predeictset] <- object$models[[index[i]]];
		}
		else
		{
			pLS[predeictset] <- predict(object$models[[index[i]]],testData[predeictset,])
		}
	}
	return (pLS);
}

GMVEBSWiMS <- function(formula = formula, data=NULL, GMVE.control = list(p.threshold = 0.85,p.samplingthreshold = 0.5), ...)
{
	if (class(formula) == "character")
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
	if (length(dependent) == 3)
	{
		Outcome = dependent[3];
	}
	outcomedata <- data[,Outcome];
	totsamples <- nrow(data);
	minSamples <- max(5,0.05*totsamples);
	clus <- NULL;
	fm <- NULL;
	baseClass <- BSWiMS.model(formula,data,...);
#	barplot(baseClass$bagging$frequencyTable);
	error <- sum(1*(baseClass$bagging$bagged.model$linear.predictors > 0.5) != outcomedata)/totsamples;

	models <- list();
	selectedfeatures <- names(baseClass$bagging$frequencyTable);
	
	if (length(baseClass$bagging$frequencyTable) > 0)
	{
		if (length(baseClass$BSWiMS.model$back.model$coefficients) >= 2)
		{
			fm <- names(baseClass$BSWiMS.model$back.model$coefficients)[-1];
		}
		else
		{
			fm <- unique(c(selectedfeatures,as.character(baseClass$univariate$Name)[1:2]));
			fm <- correlated_Remove(data,fm,thr=0.85);
		}
		
		if (length(fm) > (totsamples/10)) # we will keep the number of selected features small
		{
			fm <- fm[1:(totsamples/10)];
		}
#		print(fm)
		if (error > 0.025) # more than 2.5% of error
		{
			if (is.null(GMVE.control))
			{
				clus <- GMVECluster(as.data.frame(data[,fm]));
			}
			else
			{
				clus <- do.call(GMVECluster,c(list(as.data.frame(data[,fm])),GMVE.control));
			}
			tb <- table(clus$cluster);
			classlabels <- as.numeric(names(tb));
			if (nrow(tb) > 1)
			{
				tb <- table(clus$cluster,outcomedata);
				for (i	in 1:nrow(tb))
				{
					if (min(tb[i,]) > minSamples)
					{
							models[[i]] <- BSWiMS.model(formula,subset(data,clus$cluster == classlabels[i]),...);
					}
					else
					{
						models[[i]] <- as.numeric(colnames(tb)[which.max(tb[i,])]);
					}
				}
			}
			else
			{
					models[[1]] <- baseClass;
			}
		}
		else
		{
			models[[1]] <- baseClass;
		}
	}
	else
	{
		models[[1]] <- baseClass;
	}
	result <- list(features = fm,cluster = clus,models = models, baseModel = baseClass);
	result$selectedfeatures <- selectedfeatures;
	class(result) <- "GMVE_BSWiMS"
	return(result);
}

predict.GMVE_BSWiMS <- function(object,...)
{
	parameters <- list(...);
	testData <- parameters[[1]];
	if (!is.null(object$cluster))
	{
		if (length(object$models) > 1)
		{
			pLS <- nearestCentroid(testData[,object$features],object$cluster$centers,object$cluster$covariances,0)
			tb <- table(pLS);
			index <- as.numeric(names(tb));
			for (i in 1:nrow(tb))
			{
				predeictset <- (pLS == index[i]);
				if (class(object$models[[index[i]]])[1] == "numeric")
				{
					pLS[predeictset] <- object$models[[index[i]]];
				}
				else
				{
					pLS[predeictset] <- predict(object$models[[index[i]]],testData[predeictset,])
				}
			}
		}
		else
		{
			pLS <- predict(object$models[[1]],testData);
		}
	}
	else
	{
		pLS <- predict(object$models[[1]],testData);
	}
	return(pLS);
}
