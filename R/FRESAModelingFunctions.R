CVsignature <- function(formula = formula, data=NULL, ...)
{
	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	if (length(usedFeatures)<3) #if less than 3 features, just use a glm fit.
	{
#		print(usedFeatures)
		warning("Less than five features. Returning a glm model");
		result <- glm(formula,data=data,na.action=na.exclude,family=binomial(link=logit));
	}
	else
	{
		parameters <- list(...);
		target="All";
		CVFolds <- 0;
		repeats <- 9;
		distanceFunction=signatureDistance;
		method="pearson";
		if (!is.null(parameters$target)) target=parameters$target;
		if (!is.null(parameters$CVFolds)) CVFolds=parameters$CVFolds;
		if (!is.null(parameters$repeats)) repeats=parameters$repeats;
		if (!is.null(parameters$distanceFunction)) distanceFunction=parameters$distanceFunction;
		if (!is.null(parameters$method)) method=parameters$method;
		
		cvsig <- getSignature(data=data,varlist=usedFeatures,Outcome=baseformula[2],target,CVFolds,repeats,distanceFunction,method);
#		variable.importance <- 1:length(cvsig$featureList);
#		names(variable.importance) <- cvsig$featureList;
		
		misam <- min(cvsig$caseTamplate$samples,cvsig$controlTemplate$samples);
		lowtthr <- qt(0.40,misam-1,lower.tail = FALSE); # The distance has to greater than t values with a chance of 0.6 success
		uptthr  <- qt(0.01,misam-1,lower.tail = FALSE); # The upper distance of the t value of weights 
		ca_tmp <- cvsig$caseTamplate$template;
		co_tmp <- cvsig$controlTemplate$template;
		mvs <- as.integer((nrow(ca_tmp)+ 1.0)/2);
		wts <- cvsig$caseTamplate$meanv - cvsig$controlTemplate$meanv;
		sddCa <- (ca_tmp[mvs,]-ca_tmp[mvs-2,])/pt(1,cvsig$caseTamplate$samples-1);
		sddCo <- (co_tmp[mvs+2,]-co_tmp[mvs,])/pt(1,cvsig$controlTemplate$samples-1);
		sddP <- pmin(sddCa,sddCo);
		sddCa <- (ca_tmp[mvs+2,]-ca_tmp[mvs,])/pt(1,cvsig$caseTamplate$samples-1);
		sddCo <- (co_tmp[mvs,]-co_tmp[mvs-2,])/pt(1,cvsig$controlTemplate$samples-1);
		sddN <- pmin(sddCa,sddCo);
		sdd <- sddP;
		sdd[wts < 0] <- sddN[wts < 0];
		sdd[sdd == 0] <- 0.5*(sddP[sdd == 0] + sddN[sdd == 0]);
		sdd[sdd == 0] <- 0.5*(cvsig$controlTemplate$sdv[sdd == 0]/pt(1,cvsig$controlTemplate$samples) + cvsig$caseTamplate$sdv[sdd == 0]/pt(1,cvsig$caseTamplate$samples));
		sdd[sdd == 0] <- 0.25;
		wts <- (abs(wts)/sdd);
		variable.importance <- wts[order(-wts)];
		wts[wts > uptthr] <- uptthr;
		wts[wts < lowtthr] <- 0.01*wts[wts < lowtthr]; 
		cmat <- abs(cor(data[,colnames(cvsig$caseTamplate$template)],method="spearman"));
		cmat[cmat < 0.75] <- 0;
		cmat <- cmat*cmat;
		cwts <- sqrt(1.0/apply(cmat,1,sum));
#		print(cwts);

		wts <- wts*cwts;

		result <- list(fit=cvsig,method=method,variable.importance=variable.importance,wts=wts,cwts=cwts);
		class(result) <- "FRESAsignature";
	}
	return (result);
}

predict.FRESAsignature <- function(object, ...) 
{
	parameters <- list(...);
	wts <- NULL;
	testframe <- parameters[[1]];
	method <- object$method;
	if (!is.null(parameters$method)) method <- parameters$method;
	if (!is.null(parameters$wts)) 
	{
		wts <- parameters$wts;
	}
	else
	{
		if (inherits(object$fit$caseTamplate,"list"))
		{
			wts <- object$wts;
		}
	}
	

	controlDistances <- signatureDistance(object$fit$controlTemplate,testframe,method,wts);
	caseDistances <- signatureDistance(object$fit$caseTamplate,testframe,method,wts);
	distancep <- pnorm(controlDistances-caseDistances);
	attr(distancep,"controlDistances") <- controlDistances;
	attr(distancep,"caseDistances") <- caseDistances;
	attr(distancep,"wts") <- wts;
	return (distancep);
}

KNN_method <- function(formula = formula, data=NULL, ...)
{
	parameters <- list(...);

	baseformula <- as.character(formula);
	usedFeatures <- colnames(data)[!(colnames(data) %in% baseformula[2])]
	tlabs <- attr(terms(formula,data=data),"term.labels");
#	tlabs <- attr(terms(formula),"term.labels");
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
#	if (length(usedFeatures)<5) #if less than 5 features, just use a lm fit.
#	{
#		warning("Less than five features. Returning a lm model");
#		result <- lm(formula,data);
#	}
#	else
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
	if (inherits(cf,"list"))
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
		result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]), survival::Surv(x, y), method=method, family = "cox",...),formula = formula,usedFeatures=usedFeatures);
		bessCoefficients <- result$fit$bestmodel$coefficients
	}
	else
	{
		tb <- table(data[,baseformula[2]]);
		if (length(tb)>2)
		{
#			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "gaussian", epsilon = 1e-12,...),ic.type=ic.type,formula = formula,usedFeatures=usedFeatures);
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "gaussian",...),formula = formula,usedFeatures=usedFeatures);
		}
		else
		{
			result <- list(fit=BeSS::bess(as.matrix(data[,usedFeatures]),as.vector(data[,baseformula[2]]), method=method, family = "binomial",...), formula = formula,usedFeatures=usedFeatures);

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

BESS_EBIC <- function(formula = formula, data=NULL, ic.type="EBIC",...)
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
		if ((type == "response") || (type == "link"))
		{
			newdata <- as.data.frame(cbind(y=1:nrow(testData),testData[,object$selectedfeatures]))
			colnames(newdata) <- c("y",paste("xbest",object$selectedfeatures,sep=""))
			object$fit$bestmodel$formula <- formula(paste("y~",paste(colnames(newdata)[-1],collapse = " + ")))
			object$fit$bestmodel$terms <- terms(object$fit$bestmodel$formula)
			if (!inherits(object$fit$bestmodel,"coxph"))
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
	obj <- e1071::tune.svm(formula, data=data,gamma = gamma, cost = cost,...);
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
	if (!inherits(data[,baseformula[2]],"factor")) data[,baseformula[2]] <- as.factor(data[,baseformula[2]])
	pcaobj <- NULL;
	scaleparm <- NULL;
	numclases <- length(table(data[,baseformula[2]]))
	if (pca && (ncol(data) > 3))
	{
		outcome <- data[,baseformula[2]];
		data <- as.data.frame(data[,!(colnames(data) %in% baseformula[2])]);
		if (normalize)
		{
			scaleparm <- FRESAScale(data,method="OrderLogit");
			pcaobj <- prcomp(scaleparm$scaledData,center = FALSE,tol=0.025);
			scaleparm$scaledData <- NULL
		}
		else
		{
			pcaobj <- prcomp(data,center = FALSE,tol=0.025);
		}
		data <- as.data.frame(cbind(outcome,pcaobj$x));
		colnames(data) <- c(baseformula[2],colnames(pcaobj$x));
		data[,baseformula[2]] <- as.factor(outcome)
	}
	if (length(list(...)) == 0)
	{
		laplace = 1.0/ncol(data);
#		prior = rep(1.0/numclases,numclases);
		prior = NULL;
#		print(prior)
        fit <- try (naivebayes::naive_bayes(formula,data,prior=prior,
		laplace = laplace,
		usekernel = TRUE,
		bw="SJ",adjust=1.10,window="optcosine"), silent=TRUE)
		if (inherits(fit, "try-error"))
		{
#			print("Error. Try again")
			for (fn in colnames(data)[!(colnames(data) %in% baseformula[2])])
			{
				if (length(table(data[,fn])) < 5)
				{
					data[,fn] <- data[,fn] + rnorm(nrow(data),0,0.01*sd(data[,fn]));
				}
			}
	        fit <- try(naivebayes::naive_bayes(formula,data,prior=prior,
			laplace = laplace,
			usekernel = TRUE), silent=TRUE);
			if (inherits(fit, "try-error")) 
			{
				fit <- naivebayes::naive_bayes(formula,data,prior=prior,laplace = laplace);
			}

		}
		result <- list(fit = fit,
		pcaobj=pcaobj,
		outcome=baseformula[2],
		scaleparm=scaleparm,
		numClases=numclases);
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
		usedcolumns <- colnames(testData)[!(colnames(testData) %in% object$outcome)];
		testData <- as.data.frame(as.matrix(testData[,usedcolumns]));
		colnames(testData) <- usedcolumns;
	}
	pLS <- predict(object$fit,testData);
#	pLS <- as.numeric(as.character(predict(object$fit,testData)));
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
	if (inherits(formula, "character"))
	{
		formula <- formula(formula);
	}
	dependent <- all.vars(formula)
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[2];
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
	allClassFeatures <- character();
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(orgModel,data);
		outcomedata <- data[,Outcome];
		correct <- (abs(thePredict - outcomedata) <= (0.5 - hysteresis) );
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
						correct <- (abs(thePredict - nextdata[,Outcome])  < (0.5 - hysteresis)) ;
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
					classData[,Outcome] <- numeric(nrow(classData));
					classData[correctSet[[i]],Outcome] <- 1;
					classData[,Outcome] <- as.factor(classData[,Outcome]);
					classFeatures <- names(univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.025,thr=0.9));
					allClassFeatures <- unique(c(allClassFeatures,classFeatures))

					if (is.null(classModel.Control))
					{
						classModel[[i]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]);
					}
					else
					{
						classModel[[i]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]),classModel.Control));
					}
					classKey[correctSet[[i]]] <-  classKey[correctSet[[i]]] + 2^(i-1);
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
					baseClass = baseClass,
					allClassFeatures = allClassFeatures
					)
	class(result) <- "FRESA_HLCM"
	return(result);
}

HLCM_EM <- function(formula = formula, data=NULL,method=BSWiMS.model,hysteresis = 0.10,classMethod=KNN_method,classModel.Control=NULL,minsize=10,...)
{
	if (inherits(formula, "character"))
	{
		formula <- formula(formula);
	}
	dependent <- all.vars(formula)
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[2];
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
	outcomeFeatures <- c(Outcome,names(univariate_BinEnsemble(data=data,Outcome=Outcome,pvalue=0.01,limit = -1)));
	nselected <- selectedfeatures;

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
	fselected <- selectedfeatures;
	outcomeFeatures <- unique(c(outcomeFeatures,selectedfeatures))
	accuracy <- 1.0;
	changes <- minsize;
	n=0;
	classData <- data;
	classData[,Outcome] <- rep(0,nrow(classData));
	sselectedfeatures <- colnames(data);
	baseClass <- numeric();
	allClassFeatures <- character();
	truelimit <- 0.5 + hysteresis;
	falselimit <- 0.5 - hysteresis;
#	cat(length(outcomeFeatures),"\n");
	if (length(selectedfeatures) > 0)
	{
		thePredict <- rpredict(orgModel,data);
		outcomedata <- data[,Outcome];
		d0 <-  abs(thePredict - outcomedata);
		originalSet <- (d0 <= 0.5 );
		correct <- originalSet;
		accuracy <- sum(correct)/nrow(data);
		originalSet <- (d0 < falselimit);
		if (sum(1*(!correct),na.rm=TRUE) > minsize)
		{
			outcomedata <- data[,Outcome];
			secondSet <- (d0 >= 0.9*falselimit);
			firstSet <- (d0 <= truelimit);
			nfirstSet <- firstSet;
			nsecondSet <- secondSet;
			falseP <-  (thePredict >= 0.5) & (data[,Outcome] == 0);
			falseN <-  (thePredict <= 0.5) & (data[,Outcome] == 1);
			if ((sum(falseP) >= (minsize/2)) && (sum(falseN) >= (minsize/2)))
			{
				loops <- 0;
				firstPredict <- thePredict;
				while ((changes > (minsize/2)) && (loops < 10))
				{
					loops <- loops + 1;
					changes <- 0;
					firstSet <- nfirstSet;
					secondSet <- nsecondSet;
					firstdata <- data[firstSet,c(Outcome,selectedfeatures)];
					seconddata <- data[secondSet,outcomeFeatures];
					tb1 <- table(firstdata[,Outcome]);
					tb2 <- table(seconddata[,Outcome]);
					if ((length(tb1) > 1) && (length(tb2) > 1) && (min(tb1) >= minsize) && (min(tb2) >= minsize))
					{
						n <- 2;
						firstModel <- method(formula,firstdata,...);
						secondModel <- method(formula,seconddata,...);
					}
					else
					{
						if ((sum(secondSet) > minsize))
						{
							n <- 2;
							if (inherits(secondModel,"numeric"))
							{
								secondModel <- sum(seconddata[,Outcome])/nrow(seconddata);
#								cat("{",sum(nrow(seconddata)),"}")
							}
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
					else
					{
						firstPredict <- rpredict(firstModel,data);
						secondPredict <- rpredict(secondModel,data);
						d1 <-  abs(firstPredict - outcomedata);
						d2 <-  abs(secondPredict - outcomedata);
						nfirstSet <-  (d1 <= d2) | (d1 < falselimit);
						changes <- sum(nfirstSet != firstSet) ;
						nsecondSet <- (d2 < (d1 + hysteresis)) | (d2 < falselimit);
						changes <- changes + sum(nsecondSet != secondSet);
					}
#					cat("[(",changes,")<",length(firstModel$selectedfeatures),",",length(secondModel$selectedfeatures),">]");
				}
#				cat("[",sum(secondSet),"]")
				if (n > 0)
				{
					firstPredict <- rpredict(firstModel,data);
					secondPredict <- rpredict(secondModel,data);
					d1 <-  abs(firstPredict - outcomedata);
					d2 <-  abs(secondPredict - outcomedata);

#					cat("[",sum(!firstSet),"]")
					alternativeModel[[1]] <- firstModel;
					alternativeModel[[2]] <- secondModel;
					errorSet <- (d1 >= 0.5) & (d2 >= 0.5) & (d0 >= 0.5);
					errorfreq <- c(sum(!originalSet),sum(!firstSet),sum(!secondSet),0);
					classfreq <- c(sum(originalSet),sum(firstSet),sum(secondSet),sum(errorSet));
					baseClass <- c(0,0,0,0);
					if (!is.null(firstModel$selectedfeatures))
					{
						nselected <- secondModel$selectedfeatures;
						fselected <- firstModel$selectedfeatures;
					}
					else
					{
						if (!is.null(firstModel$bagging))
						{
							nselected <- names(secondModel$bagging$frequencyTable);
							fselected <- names(firstModel$bagging$frequencyTable);
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
#					cat("{",sumsecond,"}")
					firstPredict <- rpredict(firstModel,data);
					secondPredict <- rpredict(secondModel,data);
					d1 <-  abs(firstPredict - outcomedata);
					d2 <-  abs(secondPredict - outcomedata);
					errorSet <- (d1 >= 0.5) & (d2 >= 0.5) & (d0 >= 0.5);
					errorfreq <- c(sum(!originalSet),sum(!firstSet),sum(!secondSet),0);
					classfreq <- c(sum(originalSet),sum(firstSet),sum(secondSet),sum(errorSet));
					baseClass <- c(0,0,0,0);
#					cat("<<",sum(originalSet),",",sum(firstSet),",",sum(secondSet),",",sum(errorSet),">>") 
				}
			}
			if (n > 1)
			{
				originalSet <- (d0 <= falselimit);
				firstSet <- (d1 < 0.5);
				secondSet <- (d2 < 0.5);
				if (sum(secondSet) > minsize)
				{

	#Check for statistical significance of second set. If not do not, there is no latent class
					classData[,Outcome] <- as.factor(1*firstSet);
					classpFeatures <- univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.20,thr=0.90);
					classData[,Outcome] <- as.factor(1*secondSet);
					classpFeatures <- c(classpFeatures,univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.20,thr=0.90));
					classpFeatures <- classpFeatures[order(classpFeatures)]
					
					if (classpFeatures[1] <= 0.05) # If significant then go ahead
					{
						errorSet <- (d0 >= 0.5) & (d1 >= 0.5) & (d2 >= 0.5);
#						cat("[p=",classpFeatures[1],",",length(classpFeatures),"]<",sum(originalSet),",",sum(firstSet),",",sum(secondSet),",",sum(errorSet),">") 
#						cat("<",sum(originalSet),",",sum(firstSet),",",sum(secondSet),",",sum(errorSet),",",length(nselected),">") 
						classData[,Outcome] <- as.factor(1*originalSet);
						classpFeatures <- univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.20,thr=0.90);
						classFeatures <- names(classpFeatures);
						allClassFeatures <- unique(c(allClassFeatures,classFeatures))
						if (is.null(classModel.Control))
						{
							classModel[[1]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]);
						}
						else
						{
							classModel[[1]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]),classModel.Control));
						}
						classData[,Outcome] <- as.factor(1*firstSet);
						classpFeatures <- univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.20,thr=0.90);
						classFeatures <- names(classpFeatures);
						allClassFeatures <- unique(c(allClassFeatures,classFeatures))
						if (is.null(classModel.Control))
						{
							classModel[[2]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]);
						}
						else
						{
							classModel[[2]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]),classModel.Control));
						}
						classData[,Outcome] <- as.factor(1*secondSet);
						classpFeatures <- univariate_KS(data=classData, Outcome=Outcome,pvalue = 0.20,thr=0.90);
						classFeatures <- names(classpFeatures);
						allClassFeatures <- unique(c(allClassFeatures,classFeatures))
						if (is.null(classModel.Control))
						{
							classModel[[3]] <- classMethod(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]);
						}
						else
						{
							classModel[[3]] <- do.call(classMethod,c(list(formula(paste(Outcome,"~.")),classData[,c(Outcome,classFeatures)]),classModel.Control));
						}
						classData[,Outcome] <- 1*originalSet + 2*firstSet + 4*secondSet;
						selectedfeatures <- unique(c(selectedfeatures,nselected,fselected,allClassFeatures));
					}
					else
					{
						alternativeModel <- list();
					}
				}
				else
				{
					alternativeModel <- list();
				}
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
					baseClass = baseClass,
					allClassFeatures = allClassFeatures
					)
	class(result) <- "FRESA_HLCM"
	return(result);
}

predict.FRESA_HLCM <- function(object,...) 
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pLS <- rpredict(object$original,testData);
	pLSOr <- pLS;

	nm <- length(object$classModel);
	if (nm > 0)
	{
		prbclas <- matrix(0,nrow=nrow(testData),ncol=nm);
		for (n in 1:nm)
		{
			if (!is.null(object$classModel[[n]]))
			{
				if (class(object$classModel[[n]])[1] == "FRESAKNN")
				{
					classPred <- predict(object$classModel[[n]],testData);
					prbclas[,n] <- classPred;
				}
				else
				{
					classPred <- predict(object$classModel[[n]],testData);
					if (inherits(classPred,"numeric"))
					{
						prbclas[,n] <- classPred;
					}
					else
					{
						classPred <- predict(object$classModel[[n]],testData,probability = TRUE);
						prbclas[,n] <- attributes(classPred)$probabilities[,"1"];
						if (is.null(prbclas[,n]))
						{
							prbclas[,n] <- classPred;
						}
					}
				}
			}
		}
		prbclas[prbclas < 1.0e-10] <- 1.0e-10;
		pmodel <- pLS;
		nm <- length(object$alternativeModel);
		for (n in 1:nm)
		{
			if (!is.null(object$alternativeModel[[n]]))
			{
				ptmp <- rpredict(object$alternativeModel[[n]],testData);
				pmodel <- cbind(pmodel,ptmp);
			}
		}
		lmodels <- c(2:length(object$classModel));
		if (length(lmodels)>1)
		{
			pLS <- apply(pmodel[,lmodels]*prbclas[,lmodels],1,sum)/apply(prbclas[,lmodels],1,sum);
		}
		else
		{
			pLS <- pmodel[,lmodels];
		}
		wt <- prbclas[,1]*object$classfreq[1];
		pLS <- wt*pLSOr + (1.0-wt)*pLS;
		attr(pLS,"probabilities.class") <- prbclas;
		attr(pLS,"probabilities.model") <- pmodel;
	}
	return(pLS);
}


filteredFit <- function(formula = formula, data=NULL, 
							filtermethod=univariate_KS, 
							fitmethod=e1071::svm,
							filtermethod.control=list(pvalue=0.10,limit=0),
							Scale="none",
							PCA=FALSE,
							WHITE=c("none","CCA"),
							DECOR=FALSE,
							DECOR.control=list(thr=0.80,method="fast",type="NZLM"),
							...
						)
{

	WHITE <- match.arg(WHITE);
	if (inherits(formula, "character"))
	{
		formula <- formula(formula);
	}
	dependent <- all.vars(formula)
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[2];
	}
	fm <- colnames(data)
	fm <- fm[!(fm %in% dependent)]
	usedFeatures <-  c(Outcome,fm);

	scaleparm <- NULL;
	UPSTM <- NULL;
	pcaobj <- NULL;
	ccaobj <- NULL;
	transColnames <- NULL;

	if (DECOR && (length(fm) > 1))
	{
		if (is.null(DECOR.control))
		{
			data <- IDeA(data);
		}
		else
		{
			data <- do.call(IDeA,c(list(data),DECOR.control));
		}
		UPSTM <- attr(data,"UPSTM")
		attr(data,"UPSTM") <- NULL
		transColnames <- colnames(data);
		fm <- colnames(data)
		fm <- fm[!(fm %in% dependent)]
		usedFeatures <-  c(Outcome,fm);
	}
	if (!PCA && (WHITE=="CCA") && (Scale == "none"))
	{
		Scale <- "Norm"
	}
	if ((Scale != "none") && (length(fm) > 1) )
	{
		scaleparm <- FRESAScale(as.data.frame(data[,fm]),method=Scale);
		data[,fm] <- as.data.frame(scaleparm$scaledData);
		scaleparm$scaledData <- NULL;
	}
	filtout <- filtermethod
	if (!is.null(filtermethod))
	{
		
		if (is.null(filtermethod.control))
		{
			fm <- filtermethod(data,formula);
		}
		else
		{
			fm <- do.call(filtermethod,c(list(data,formula),filtermethod.control));
		}
		filtout <- fm;
		if (length(fm) > 1)
		{
			selpvalue <- sqrt(median(fm)*(min(fm)+1.0e-12));
			maxpvalue <- selpvalue*1.0e3 + 1.0e-6;
			fm <- fm[fm <= maxpvalue];
		}
		fm <- names(fm)
		usedFeatures <-  c(Outcome,fm);
	}
	data <- data[,usedFeatures]


	binOutcome <- length(table(data[,Outcome])) == 2
	isFactor <- inherits(data[,Outcome], "factor")
	if (PCA && (length(fm) > 1))
	{
		if (binOutcome)
		{
			controlSet <- subset(data,get(Outcome) == 0)
			if (length(usedFeatures) > 2)
			{
				pcaobj <- prcomp(controlSet[,fm],center = (Scale == "none"), scale.= (Scale == "none"),tol=0.025);
				data <- as.data.frame(cbind(data[,Outcome],as.data.frame(predict(pcaobj,data[,fm]))));
				colnames(data) <- c(Outcome,colnames(pcaobj$x));
				if (isFactor)
				{
					data[,Outcome] <-as.factor(data[,Outcome])
				}
			}
		}
		else
		{
			pcaobj <- prcomp(data[,fm],center = (Scale == "none"), scale.= (Scale == "none"),tol=0.025);
			data <- as.data.frame(cbind(data[,Outcome],as.data.frame(pcaobj$x)));
			colnames(data) <- c(Outcome,colnames(pcaobj$x));
			if (isFactor)
			{
				data[,Outcome] <-as.factor(data[,Outcome])
			}
		}
	}
	if (!PCA && (WHITE=="CCA") && (length(fm) > 1))
	{
		if (!requireNamespace("whitening", quietly = TRUE)) {
			 install.packages("whitening", dependencies = TRUE)
		}
		mx <- as.matrix(data[,fm]);
		ccaobj <- whitening::scca(mx, mx,verbose=FALSE);
		ccaobj$WY <- NULL
		CCAX <- as.data.frame(tcrossprod( mx, ccaobj$WX ))
		data <- as.data.frame(cbind(data[,Outcome],CCAX));
		colnames(data) <-  c(Outcome,colnames(CCAX));
		if (isFactor)
		{
			data[,Outcome] <-as.factor(data[,Outcome])
		}
#		cat(colnames(data));
	}
	
	fit <- try(fitmethod(formula,data,...));
	selectedfeatures <- fm;
	if (!is.null(fit$selectedfeatures))
	{
		selectedfeatures <- fit$selectedfeatures;
	}
	else
	{
		 vs <- NULL;
          if (!is.null(fit$importance))
          {
            vs <- fit$importance[,1];
          }
          if (!is.null(fit$variable.importance))
          {
            vs <- fit$variable.importance;
          }
          if (!is.null(vs))
          {
            #						print(vs)
            if (length(vs)>1)
            {
              if (length(vs) > 2) vs <- vs[order(-vs)];
              if (sum(vs > 0.0) > 1)
              {
                vs <- vs[vs > 0.0];							
              }
              selectedfeatures <- names(vs);
            }
          }
	}

	
	parameters <- list(...);
	result <- list(fit=fit,
					filter=filtout,
					processedFeatures = fm,
					selectedfeatures = selectedfeatures,
					usedFeatures = usedFeatures,
					parameters=parameters,
					asFactor=(class(data[,Outcome])=="factor"),
					classLen=length(table(data[,Outcome])),
					Scale = scaleparm,
					binOutcome = binOutcome,
					Outcome = Outcome,
					pcaobj = pcaobj,
					ccaobj = ccaobj,
					UPSTM = UPSTM,
					transColnames = transColnames
					);
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
	if (!is.null(object$UPSTM))
	{
	    testData[,rownames(object$UPSTM)] <- Rfast::mat.mult(as.matrix(testData[,rownames(object$UPSTM)]),object$UPSTM);
		colnames(testData) <- object$transColnames;
	}
	if (!is.null(object$Scale))
	{
		testData <- FRESAScale(as.data.frame(testData),
								method=object$Scale$method,
								refMean=object$Scale$refMean,
								refDisp=object$Scale$refDisp
							  )$scaledData;
	}
	testData <- as.data.frame(testData[,object$usedFeatures])
	if (!is.null(object$pcaobj))
	{
		pcapred <- predict(object$pcaobj,testData[,object$processedFeatures]);
		testData <- as.data.frame(cbind(testData[,object$usedFeatures[1]],pcapred));
		colnames(testData) <-  c(object$usedFeatures[1],colnames(pcapred));
	}
	if (!is.null(object$ccaobj))
	{
		mx <- as.matrix(testData[,object$processedFeatures]);
		CCAX <- as.data.frame(tcrossprod( mx, object$ccaobj$WX ))
		testData <- as.data.frame(cbind(testData[,object$usedFeatures[1]],CCAX));
		colnames(testData) <-  c(object$usedFeatures[1],colnames(CCAX));
#		cat(colnames(testData));
	}
	
	probability <- FALSE;
	if (!is.null(object$parameters$probability))
	{
		probability <- object$parameters$probability;
	}
	pLS <- rpredict(object$fit,testData,asFactor=object$asFactor,classLen=object$classLen,probability=probability,...);
	return (pLS);
}

ClustClass <- function(formula = formula, data=NULL, filtermethod=univariate_KS, clustermethod=GMVECluster, classmethod=LASSO_1SE,filtermethod.control=list(pvalue=0.1,limit=21),clustermethod.control=list(p.threshold = 0.95,p.samplingthreshold = 0.5),classmethod.control=list(family = "binomial"),pca=TRUE,normalize=TRUE)
{
	if (inherits(formula, "character"))
	{
		formula <- formula(formula);
	}
	dependent <- all.vars(formula)
	Outcome = dependent[1];
	if (length(dependent) == 3)
	{
		Outcome = dependent[2];
	}

	outcomedata <- data[,Outcome];
	totsamples <- nrow(data);
	minSamples <- max(5,0.025*totsamples);
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
	selectedfeatures <- names(fm);
	datapca <- as.data.frame(data[,selectedfeatures]);
	pcaobj <- NULL;
	scaleparm <- NULL;
	if (pca &&  (ncol(datapca) > 3))
	{
		rank. = max(3,as.integer(length(selectedfeatures)/3+0.5));
		rank. = min(rank.,length(selectedfeatures))
		if (normalize)
		{
			scaleparm <- FRESAScale(datapca,method="OrderLogit");
			pcaobj <- prcomp(scaleparm$scaledData,center = FALSE,rank.=rank.,tol=0.05);
			scaleparm$scaledData <- NULL
		}
		else
		{
			pcaobj <- prcomp(datapca,center = TRUE,rank.=rank.,tol=0.05);
		}
		datapca <- as.data.frame(pcaobj$x);
	}

	if (is.null(clustermethod.control))
	{
		clus <- clustermethod(datapca);
	}
	else
	{
		clus <- do.call(clustermethod,c(list(datapca),clustermethod.control));
	}
#	print(selectedfeatures)
	tb <- table(clus$classification);
	classlabels <- as.numeric(names(tb));
	models <- list();
#	data <- data[,c(Outcome,selectedfeatures)];
	allfeatures <- selectedfeatures;
	if (length(classlabels) > 1)
	{
			tb <- table(clus$classification,outcomedata);
#			print(tb)
			for (i	in 1:nrow(tb))
			{
#				cat(i,":")
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
					allfeatures <- c(allfeatures,models[[i]]$selectedFeatures);
				}
				else
				{
					models[[i]] <- as.numeric(colnames(tb)[which.max(tb[i,])]);
				}
#				cat(tb[i,],"<")
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
		allfeatures <- c(allfeatures,models[[1]]$selectedFeatures);
	}
	result <- list(features = selectedfeatures,cluster = clus,models = models,pcaobj = pcaobj,scaleparm=scaleparm );
	result$selectedfeatures <- allfeatures;

	class(result) <- "CLUSTER_CLASS"
	return(result);
}

predict.CLUSTER_CLASS <- function(object,...)
{
	parameters <- list(...);
	testData <- parameters[[1]];
	pcaData <- testData[,object$features];
	if (!is.null(object$pcaobj))
	{
		if (!is.null(object$scaleparm))
		{
			pcaData <- FRESAScale(pcaData,method=object$scaleparm$method,refMean=object$scaleparm$refMean,refDisp=object$scaleparm$refDisp)$scaledData;
		}
		pcaData <- as.data.frame(predict(object$pcaobj,pcaData));
	}
	pLS <- predict(object$cluster,pcaData)$classification;
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
	return (pLS);
}

GMVEBSWiMS <- function(formula = formula, data=NULL, GMVE.control = list(p.threshold = 0.95,p.samplingthreshold = 0.5), ...)
{
	if (inherits(formula, "character"))
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
	selectedfeatures <- baseClass$selectedfeatures;
#	print(selectedfeatures);
	fm <- names(baseClass$BSWiMS.model$at.opt.model$coefficients)[-1]
#			print(fm)
	if (length(fm) > 0)
	{		
		if (error > 0.025) # more than 2.5% of error
		{
			fm <- unique(c(fm,names(univariate_KS(data,Outcome,pvalue=0.05,limit=10,thr=0.8))));
			selectedfeatures <- unique(c(fm,selectedfeatures));
#			print(selectedfeatures);
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
							selectedfeatures <- unique(c(selectedfeatures,names(models[[i]]$bagging$frequencyTable)));

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
