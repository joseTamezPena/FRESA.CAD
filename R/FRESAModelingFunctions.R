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
		kn <- 2*(as.integer(sqrt(min(tb))/2 + 0.5)) + 1;
	}
	else
	{
		kn <- parameters$kn;
	}
#	print(usedFeatures);


	scaledData <- as.data.frame(data[,usedFeatures]);
	colnames(scaledData) <- usedFeatures;
	rownames(scaledData) <- rownames(data);
	traindata <- scaledData;

	scaleMethod <- "None";
	if (is.null(parameters$scaleMethod))
	{
		scaleMethod <- "OrderLogit";
	}
	else
	{
		scaleMethod	<- parameters$scaleMethod;
	}
	mean_vec <- NULL;
	disp_vec <- NULL;
	if (scaleMethod != "None") 
	{
		scaledParam <- FRESAScale(scaledData,method=scaleMethod);
		mean_vec <- scaledParam$refMean;
		disp_vec <- scaledParam$refDisp;
		scaledData <- scaledParam$scaledData;
	}
	
	result <- list(trainData=traindata,scaledData=scaledData,classData=data[,baseformula[2]],outcome=baseformula[2],usedFeatures=usedFeatures,mean_col=mean_vec,disp_col=disp_vec,kn=kn,scaleMethod=scaleMethod);
	result$selectedfeatures <- usedFeatures;
	class(result) <- "FRESAKNN"
	return(result);
}


predict.FRESAKNN <- function(object, ...) 
{
	parameters <- list(...);
	testframe <- parameters[[1]];

	testframe <- as.data.frame(testframe[,object$usedFeatures]);
	colnames(testframe) <- object$usedFeatures;
	trainframe <- object$scaledData
	
	if (object$scaleMethod != "None")
	{
		testframe <- FRESAScale(testframe,object$trainData,object$scaleMethod,object$mean_col,object$disp_col)$scaledData;
	}

	knnclass <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn,prob=TRUE))
	if (inherits(knnclass, "try-error")) knnclass <- numeric(nrow(testframe));
	if (length(table(object$classData))==2)
	{
		classk <- knnclass;
		prop <- attributes(knnclass);
		knnclass <- abs(prop$prob - 1*(knnclass == "0"))
		if (object$kn > 3)
		{
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn - 2,prob=TRUE))
			prop_1 <- attributes(knnclass_1);
			knnclass_2 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn + 2,prob=TRUE))
			prop_2 <- attributes(knnclass_2);
			knnclass <- (knnclass + 0.5*abs(prop_1$prob - 1*(knnclass_1 == "0")) + 0.5*abs(prop_2$prob - 1*(knnclass_2 == "0")))/2.0;
		}
		attr(knnclass,"class") <- as.character(classk);
	}
	else
	{
		prop <- attributes(knnclass)$prob;
		if (object$kn > 3)
		{
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn - 2,prob=TRUE))
			prop_1 <- attributes(knnclass_1)$prob;
			testclass <- as.character(knnclass_1) == as.character(knnclass)
			prop[testclass] <- 0.75*prop[testclass] + 0.25*prop_1[testclass];
#			prop[!testclass] <- 0.75*prop[!testclass] + 0.25*(1.0 - prop_1[!testclass]);
			knnclass_1 <- try(class::knn(trainframe,testframe,factor(object$classData),object$kn + 2,prob=TRUE))
			prop_2 <- attributes(knnclass_1)$prob;
			testclass_2 <- as.character(knnclass_1) == as.character(knnclass)
			prop[testclass_2] <- 0.75*prop[testclass_2] + 0.25*prop_2[testclass_2];
#			prop[!testclass_2] <- 0.75*prop[!testclass_2] + 0.25*(1.0 - prop_2[!testclass_2]);
		}
		attr(knnclass,"prob") <- prop;
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
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "gaussian", epsilon = 1e-12,...),ic.type=ic.type,formula = formula,usedFeatures=usedFeatures);
		}
		else
		{
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "binomial", epsilon = 0,...),ic.type=ic.type, formula = formula,usedFeatures=usedFeatures);

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
			if (class(object$fit$bestmodel) != "coxph")
			{
				pLS <- predict(object$fit$bestmodel,newdata,type=type);
			}
			else
			{
				pLS <- predict(object$fit$bestmodel,newdata);
			}
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

NAIVE_BAYES <- function(formula = formula, data=NULL,pca=TRUE,normalize=TRUE,...)
{
if (!requireNamespace("naivebayes", quietly = TRUE)) {
	 install.packages("naivebayes", dependencies = TRUE)
} 
	baseformula <- as.character(formula);
	if (class(data[,baseformula[2]]) != "factor") data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	pcaobj <- NULL;
	scaleparm <- NULL;
	numclases <- length(table(data[,baseformula[2]]))
	if (pca && (nrow(data) > 2*ncol(data)) && (ncol(data) > 3))
	{
		outcome <- data[,baseformula[2]];
		data <- as.data.frame(data[,!(colnames(data) %in% baseformula[2])]);
		if (normalize)
		{
			scaleparm <- FRESAScale(data,method="OrderLogit");
			pcaobj <- prcomp(scaleparm$scaledData);
		}
		else
		{
			pcaobj <- prcomp(data);
		}
		data <- as.data.frame(cbind(as.numeric(as.character(outcome)),pcaobj$x));
		colnames(data) <- c(baseformula[2],colnames(pcaobj$x));
		data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	}
	if (length(list(...)) == 0)
	{
		result <- list(fit = naivebayes::naive_bayes(formula,data,usekernel = TRUE,bw="SJ",adjust=0.5),pcaobj=pcaobj,outcome=baseformula[2],scaleparm=scaleparm,numClases=numclases);
	}
	else
	{
		result <- list(fit = naivebayes::naive_bayes(formula,data,...),pcaobj=pcaobj,outcome=baseformula[2],scaleparm=scaleparm,numClases=numclases);
	}
	class(result) <- "FRESA_NAIVEBAYES"
	return(result);
}


predict.FRESA_NAIVEBAYES <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	if (!is.null(object$pcaobj))
	{
		if (!is.null(object$scaleparm))
		{
			testData <- FRESAScale(testData,method=object$scaleparm$method,refMean=object$scaleparm$refMean,refDisp=object$scaleparm$refDisp)$scaledData;
		}
		testData <- predict(object$pcaobj,testData);
	}
	else
	{
		testData <- as.data.frame(testData[,!(colnames(testData) %in% object$outcome)]);
	}
	pLS <- as.numeric(as.character(predict(object$fit,testData)));
	if (is.null(parameters$probability))
	{
		if (object$numClases == 2)
		{
			prop <- predict(object$fit,testData,type = "prob");
			pLS <- prop[,"1"];
			pLS[is.nan(pLS)] <- 0.5;
			pLS[is.na(pLS)] <- 0.5;
		}
		else
		{
			attr(pLS,"prob") <- predict(object$fit,testData,type = "prob");
		}
	}
	else
	{
		attr(pLS,"probabilities") <- predict(object$fit,testData,type = "prob");
	}
	if (!is.null(object$pcaobj))
	{
		attr(pLS,"PCAData") <- testData;
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


HLCM <- function(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.1,classMethod=KNN_method,classModel.Control=NULL,minsize=10,...)
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
	classModel <- list();
	selectedfeatures <- colnames(data)[!(colnames(data) %in% Outcome)]
	orgModel <- try(method(formula,data,...));
	if (inherits(orgModel, "try-error"))
	{
		orgModel <- 0.5;
	}
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
	classData[,Outcome] <- numeric(nrow(classData));
	classKey <- numeric(nrow(classData));
	names(classKey) <- rownames(data);
	errorfreq <- numeric();
	classfreq <- numeric();
	baseClass <- numeric();
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(orgModel,data);
		outcomedata <- data[,Outcome];
		correct <- ((thePredict >= 0.5) == (outcomedata > 0));
		accuracy <- sum(correct)/nrow(data);
		toterror <- sum(!correct);
		baseClass <- c(baseClass,0);
		correctSet[[n+1]] <- rownames(data[correct,]);
		classfreq <- c(classfreq,length(correctSet[[n+1]]));
		errorfreq <- c(errorfreq,toterror);
		if (toterror > minsize)
		{
			nextdata <- data;
			while (inserted && (toterror > minsize) && (toterror < (nrow(nextdata) - 1)))
			{
				inserted <- FALSE;
				preData <- nextdata;
				outcomedata <- preData[,Outcome];
				falseP <- (thePredict >= (0.5 - hysteresis)) & (outcomedata == 0);
				falseN <- (thePredict <= (0.5 + hysteresis)) & (outcomedata == 1);
				incorrectSet <- falseP | falseN;
				if (sum(incorrectSet) > nrow(nextdata)/2)
				{
					falseP <- (thePredict >= 0.5) & (outcomedata == 0);
					falseN <- (thePredict <= 0.5) & (outcomedata == 1);
					incorrectSet <- falseP | falseN;
				}
				if ((sum(falseP) >= (minsize/2)) && (sum(falseN) >= (minsize/2)))
				{
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
						n <- n + 1;
						baseClass <- c(baseClass,n+1);
						selectedfeatures <- c(selectedfeatures,nselected);
						selectedfeatures <- unique(selectedfeatures);
						cat("[",sum(incorrectSet),"]");
						alternativeModel[[n]] <- alternativeM;
						thePredict <- rpredict(alternativeM,nextdata);
						correct <- ((thePredict >= 0.5) == (nextdata[,Outcome] > 0));
						correctSet[[n+1]] <- rownames(nextdata[correct,]);
						classfreq <- c(classfreq,length(correctSet[[n+1]]));
						toterror <- sum(abs(thePredict - nextdata[,Outcome]) > 0.5 );
						errorfreq <- c(errorfreq,toterror);
						inserted <- TRUE;
					}	
				}
			}
			if ( !inserted )
			{
				cat("<",sum(incorrectSet),">")
				if (sum(incorrectSet) > minsize) 
				{
					n <- n + 1;
					baseClass <- c(baseClass,n);
					cat("(",sum(incorrectSet),")")
					alternativeModel[[n]] <- (sum(preData[incorrectSet,Outcome])/nrow(preData[incorrectSet,]));
					correctSet[[n+1]] <- rownames(preData[incorrectSet,]);
					classfreq <- c(classfreq,length(correctSet[[n+1]]));
					errorfreq <- c(errorfreq,0);
				}
			}
			if (n > 0)
			{
				for (i in 1:(n+1))
				{
					classData[,Outcome] <- data[,Outcome];
					classData[correctSet[[i]],Outcome] <- classData[correctSet[[i]],Outcome] + 2;
					if (is.null(classModel.Control))
					{
						classModel[[i]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
					}
					else
					{
						classData[,Outcome] <- as.factor(classData[,Outcome]);
						classModel[[i]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
					}
					classKey[correctSet[[i]]] <- classKey[correctSet[[i]]] + 2^(i-1);
				}
			}
		}
	}
	errorfreq <- errorfreq/nrow(data);
	classfreq <- classfreq/nrow(data);
	result <- list(original = orgModel,
					alternativeModel = alternativeModel,
					classModel = classModel,
					accuracy = accuracy,
					selectedfeatures = selectedfeatures,
					hysteresis = hysteresis,
					classSet = classKey,
					errorfreq = errorfreq,
					classfreq = classfreq,
					baseClass = baseClass
					)
	class(result) <- "FRESA_HLCM"
	return(result);
}

HLCM_EM <- function(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.1,classMethod=KNN_method,classModel.Control=NULL,minsize=10,...)
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

	errorfreq <- numeric();
	classfreq <- numeric();
	alternativeModel <- list();
	correctSet <- NULL;
	firstSet <- NULL;
	secondSet <- NULL;
	originalSet <- NULL;
	classModel <- list();
	lastModel <- 0.5;
	firstModel <- 0.5;
	secondModel <- 0.5;
	selectedfeatures <- colnames(data)[!(colnames(data) %in% Outcome)]
	orgModel <- try(method(formula,data,...));
	if (inherits(orgModel, "try-error"))
	{
		warning("Error. Setting prediction to 0.5\n")
		orgModel <- 0.5;
	}
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
	changes <- 1;
	n=0;
	classData <- data;
	classData[,Outcome] <- rep(0,nrow(classData));
	sselectedfeatures <- colnames(data);
	baseClass <- numeric();
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(orgModel,data);
		outcomedata <- data[,Outcome];
		correct <- ((thePredict >= 0.5) == (outcomedata > 0));
		originalSet <- correct;
		accuracy <- sum(correct)/nrow(data);
		if (sum(1*(!correct),na.rm=TRUE) > minsize)
		{
			outcomedata <- data[,Outcome];
			falseP <- (thePredict >= (0.5 - hysteresis)) & (outcomedata == 0);
			falseN <- (thePredict <= (0.5 + hysteresis)) & (outcomedata == 1);
			secondSet <- falseP | falseN;
			if ( sum(secondSet) > (2*nrow(data)/3) )
			{
				falseP <- (thePredict >= 0.5) & (outcomedata == 0);
				falseN <- (thePredict <= 0.5) & (outcomedata == 1);
			}
			trueP <- (thePredict >= (0.5 - hysteresis)) & (outcomedata == 1);
			trueN <- (thePredict <= (0.5 + hysteresis)) & (outcomedata == 0);
			firstSet <- trueP | trueN;
			if ((sum(falseP) >= (minsize/2)) && (sum(falseN) >= (minsize/2)))
			{
				loops <- 0;
				firstPredict <- thePredict;
				while ((changes > 0) && (loops < 10))
				{
					loops <- loops + 1;
					n <- 1;
					changes <- 0;
					firstdata <- data[firstSet,];
					seconddata <- data[secondSet,];
					tb1 <- table(firstdata[,Outcome]);
					tb2 <- table(seconddata[,Outcome]);
					if ((length(tb1) > 1) && (length(tb2) > 1) && (min(tb1) > minsize) && (min(tb2) > minsize))
					{
						firstModel <- method(formula,firstdata,...);
						alternativeModel[[1]] <- firstModel;
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
						if (length(nselected) > 0)
						{
							n <- 2;
							firstPredict <- rpredict(firstModel,data);
							secondPredict <- rpredict(secondModel,data);
							d1 <-  abs(firstPredict - outcomedata);
							d2 <-  abs(secondPredict - outcomedata);
							nfirstSet <- (d1 <= (d2 + hysteresis));
							changes <- sum(nfirstSet != firstSet);
							firstSet <- nfirstSet;
							nsecondSet <- (d2 <= (d1 + hysteresis));
							changes <- changes + sum(nsecondSet != secondSet);
							secondSet <- nsecondSet;
						}
					}
					else
					{
						if ((sum(secondSet) > minsize))
						{
							n <- 2;
							secondModel <- sum(seconddata[,Outcome])/nrow(seconddata);
							secondPredict <- rpredict(secondModel,data);
							cat("{",sum(nrow(seconddata)),"}")
						}
						else
						{
							secondModel <- 0.5;
						}
					}
					if (sum(secondSet) == 0)
					{
						changes <- 0;
						secondModel <- 0.5;
					}
					cat("(",changes,")");
				}
#				cat("[",sum(secondSet),"]")
				if (n > 0)
				{
					firstPredict <- rpredict(firstModel,data);
					secondPredict <- rpredict(secondModel,data);
					d0 <-  abs(thePredict - outcomedata);
					d1 <-  abs(firstPredict - outcomedata);
					d2 <-  abs(secondPredict - outcomedata);
					firstSet <- (d1 < 0.5);
					secondSet <- (d2 < 0.5);
					originalSet <- (d0 < 0.5);

#					cat("[",sum(!firstSet),"]")
					alternativeModel[[1]] <- firstModel;
					alternativeModel[[2]] <- secondModel;
					errorSet <- (d1 >= 0.5) & (d2 >= 0.5);
					errorfreq <- c(sum(!originalSet),sum(!firstSet),sum(!secondSet),0);
					classfreq <- c(sum(originalSet),sum(firstSet),sum(secondSet),sum(errorSet));
					baseClass <- c(0,0,0,0);
					cat("<",sum(originalSet),",",sum(firstSet),",",sum(secondSet),",",sum(errorSet),">") 
					nselected <- character();
					if (sum(errorSet) > minsize)
					{
						cat("%",sum(errorSet),"%")
						errordata <- data[errorSet,c(Outcome,selectedfeatures)];
						tb <- table(errordata[,Outcome]);
						if ((length(tb) > 1) && (min(tb) > (minsize/2)))
						{
							lastModel <- method(formula,errordata,...);
							if (!is.null(lastModel$selectedfeatures))
							{
								nselected <- lastModel$selectedfeatures;
							}
							else
							{
								if (!is.null(lastModel$bagging))
								{
									nselected <- names(lastModel$bagging$frequencyTable);
								}
							}
							if (length(nselected) > 0)
							{
								n = 3;
								alternativeModel[[3]] <- lastModel;
							}
						}
						else
						{
							n = 3;
							alternativeModel[[3]] <- sum(errordata[,Outcome])/nrow(errordata);
						}
					}
				}
			}
			else
			{
				sumsecond <- sum(secondSet)
				if (sumsecond > minsize)
				{
					n <- 2;
					alternativeModel[[1]] <- firstModel;
					alternativeModel[[2]] <- sum(data[secondSet,Outcome])/sumsecond;
					cat("{",sumsecond,"}")
					firstPredict <- rpredict(firstModel,data);
					secondPredict <- rpredict(secondModel,data);
					d0 <-  abs(thePredict - outcomedata);
					d1 <-  abs(firstPredict - outcomedata);
					d2 <-  abs(secondPredict - outcomedata);
					firstSet <- (d1 < 0.5);
					secondSet <- (d2 < 0.5);
					originalSet <- (d0 < 0.5);
					errorSet <- (d1 >= 0.5) & (d2 >= 0.5);
					errorfreq <- c(sum(!originalSet),sum(!firstSet),sum(!secondSet),0);
					classfreq <- c(sum(originalSet),sum(firstSet),sum(secondSet),sum(errorSet));
					baseClass <- c(0,0,0,0);
					cat("<<",sum(originalSet),",",sum(firstSet),",",sum(secondSet),",",sum(errorSet),">>") 
				}
			}
			if (n > 0)
			{
				nselected <- character();
				if (class(firstModel)[1] != "numeric")
				{
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
				}
				nselected <- character();
				if (class(secondModel)[1] != "numeric")
				{
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
				}
				selectedfeatures <- unique(selectedfeatures);
				classData[,Outcome] <- 2*originalSet + outcomedata;
				classData[,Outcome] <- as.factor(classData[,Outcome]);
				if (is.null(classModel.Control))
				{
					classModel[[1]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
				}
				else
				{
					classModel[[1]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
				}
				classData[,Outcome] <- 2*firstSet + outcomedata;
				classData[,Outcome] <- as.factor(classData[,Outcome]);
				if (is.null(classModel.Control))
				{
					classModel[[2]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
				}
				else
				{
					classModel[[2]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
				}
				if (n > 1)
				{
					classData[,Outcome] <- 2*secondSet + outcomedata;
					classData[,Outcome] <- as.factor(classData[,Outcome]);
					if (is.null(classModel.Control))
					{
						classModel[[3]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
					}
					else
					{
						classModel[[3]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
					}
					if ( n > 2)
					{
						classData[,Outcome] <- 2*errorSet + outcomedata;
						classData[,Outcome] <- as.factor(classData[,Outcome]);
						if (is.null(classModel.Control))
						{
							classModel[[4]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]);
						}
						else
						{
							classModel[[4]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,selectedfeatures)]),classModel.Control));
						}
					}
				}
				classData[,Outcome] <- numeric(nrow(classData));
				classData[originalSet,Outcome] <- 1;
				classData[firstSet,Outcome] <- classData[firstSet,Outcome] + 2;
				classData[secondSet,Outcome] <- classData[secondSet,Outcome] + 4;
				classData[errorSet,Outcome] <- classData[errorSet,Outcome] + 8;
			}
			else
			{
				alternativeModel <- list();
			}
		}
	}
	errorfreq <- errorfreq/nrow(data);
	classfreq <- classfreq/nrow(data);
	result <- list(original = orgModel,
					alternativeModel = alternativeModel,
					classModel = classModel,
					accuracy = accuracy,
					selectedfeatures = selectedfeatures,
					hysteresis = hysteresis,
					classSet = classData[,Outcome],
					errorfreq = errorfreq,
					classfreq = classfreq,
					baseClass = baseClass
					)
	class(result) <- "FRESA_HLCM"
	return(result);
}

predict.FRESA_HLCM <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pLS <- rpredict(object$original,testData);

	if (length(object$classModel) > 0)
	{
		prbclas <- matrix(0,nrow=nrow(testData),ncol=length(object$classModel));
		for (n in 1:length(object$classModel))
		{
			if (class(object$classModel[[n]])[1] == "FRESAKNN")
			{
				classPred <- predict(object$classModel[[n]],testData);
				if (class(classPred) == "factor")
				{
#					print(table(object$classModel[[n]]$classData))
					object$classModel[[n]]$classData <- as.integer(as.numeric(as.character(object$classModel[[n]]$classData))/2);
#					print(table(object$classModel[[n]]$classData))
					classPred2 <- as.numeric(predict(object$classModel[[n]],testData));
#					print(table(classPred2))
					nclass <- as.numeric(as.character(classPred));
					prbclas[,n] <- attributes(classPred)$prob*(nclass > 1) + (1.0 - attributes(classPred)$prob)*(nclass < 2);
					condone <- (classPred2 >= 0.5) & (classPred2 > prbclas[,n]);
					prbclas[condone,n] <- classPred2[condone];
					condtwo <- (classPred2 < 0.5) & (classPred2 < prbclas[,n]);
					prbclas[condtwo,n] <- classPred2[condtwo];
				}
				else
				{
					nclass <- as.numeric(attributes(classPred)$class);
					prbclas[,n] <- classPred*(nclass == 3) + (1.0 - classPred)*(nclass == 2);
				}
			}
			else
			{
				classPred <- predict(object$classModel[[n]],testData,probability = TRUE);
				classnames <- colnames(attributes(classPred)$probabilities)
				prbclas[,n] <- 0;
				if ("2" %in% classnames)
				{
					prbclas[,n] <- attributes(classPred)$probabilities[,"2"];
				}
				if ("3" %in% classnames)
				{
					prbclas[,n] <- prbclas[,n] + attributes(classPred)$probabilities[,"3"];
				}
			}
		}
		pmodel <- pLS;
		nm <- length(object$alternativeModel);
		for (n in 1:nm)
		{
			ptmp <- rpredict(object$alternativeModel[[n]],testData);
			pmodel <- cbind(pmodel,ptmp);
		}
		nm <- length(object$classModel);
		mpclas <- 1.0*(apply(prbclas,1,max) > 0.75);
		for (i in 1:length(pLS))
		{
			nwt <- object$classfreq[1]*(1.0 - prbclas[i,1])*mpclas[i];
			wts <- prbclas[i,1] + nwt;
			pLS[i] <- prbclas[i,1]*pmodel[i,1] + nwt*(1.0-pmodel[i,1]);
			if (nm > 1)
			{
				for (n in 2:nm)
				{
					if ((object$baseClass[n] == 0) || (object$baseClass[n] > nm))
					{
						nwt <- object$classfreq[n]*(1.0 - prbclas[i,n])*mpclas[i];
					}
					else
					{
						nwt <- object$classfreq[n]*prbclas[i,object$baseClass[n]]*mpclas[i];
					}
					wts <- wts + prbclas[i,n]*mpclas[i] + nwt;
					pLS[i] <- pLS[i] + prbclas[i,n]*pmodel[i,n]*mpclas[i] + nwt*(1.0-pmodel[i,n]);
				}
			}
			if (wts > 0) pLS[i] <- pLS[i]/wts;
		}
		attr(pLS,"probabilities") <- prbclas;
	}
	return(pLS);
}


filteredFit <- function(formula = formula, data=NULL, filtermethod=univariate_Wilcoxon, fitmethod=e1071::svm,filtermethod.control=list(pvalue=0.10,limit=0.1),...)
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
	fit <- try(fitmethod(formula,data[,usedFeatures],...));
	parameters <- list(...);
	result <- list(fit=fit,filter=fm,selectedfeatures = names(fm),usedFeatures = usedFeatures,parameters=parameters,asFactor=(class(data[,Outcome])=="factor"),classLen=length(table(data[,Outcome])));
	class(result) <- c("FRESA_FILTERFIT");
	if (inherits(fit, "try-error"))
	{
		warning("Fit error\n")
		class(result) <- c("FRESA_FILTERFIT","try-error");
	}
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

ClustClass <- function(formula = formula, data=NULL, filtermethod=univariate_Wilcoxon, clustermethod=GMVECluster, classmethod=LASSO_1SE,filtermethod.control=list(pvalue=0.1,limit=10),clustermethod.control=list(p.threshold = 0.95,p.samplingthreshold = 0.5),classmethod.control=list(family = "binomial"))
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

GMVEBSWiMS <- function(formula = formula, data=NULL, GMVE.control = list(p.threshold = 0.95,p.samplingthreshold = 0.5), ...)
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
	error <- sum(1*(baseClass$bagging$bagged.model$linear.predictors > 0.0) != outcomedata)/totsamples;
#	cat(error)

	models <- list();
	selectedfeatures <- names(baseClass$bagging$frequencyTable);
	fm <- selectedfeatures
#			print(fm)
	if (length(fm) > 0)
	{		
		if (error > 0.025) # more than 2.5% of error
		{
			fm <- names(univariate_Wilcoxon(data,Outcome,pvalue=0.05,limit=10));
			selectedfeatures <- unique(fm,selectedfeatures);
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
							selectedfeatures <- unique(selectedfeatures,names(models[[i]]$bagging$frequencyTable));

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
			pLS <- predict(object$cluster,testData)$classification
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
