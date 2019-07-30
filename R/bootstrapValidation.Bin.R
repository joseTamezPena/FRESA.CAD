bootstrapValidation_Bin <-
function(fraction=1.00,loops=200,model.formula,Outcome,data,type=c("LM","LOGIT","COX"),plots=FALSE,best.model.formula=NULL)
{

	if (class(model.formula) == "character")
	{
		model.formula = formula(model.formula);
	}
	varsList <- all.vars(model.formula)
	basemodel <- modelFitting(model.formula,data,type,TRUE);
	

	bootmodel <- basemodel;
	if (is.null(best.model.formula))
	{
		best.model.formula = model.formula;
	}
	else
	{
		if (class(best.model.formula) == "character")
		{
			best.model.formula = formula(best.model.formula);
		}
	}
	
	varsList2 <- all.vars(best.model.formula)

    svar <- 1+1*(type=="COX")
	epts <- c(svar:length(basemodel$coefficients));
	
	if ((length(varsList)>svar)&&(length(varsList2)>svar))
	{

		type <- match.arg(type);
		casesample = subset(data,get(Outcome)== 1);
		controlsample = subset(data,get(Outcome) == 0);
		sizecases = nrow(casesample);
		sizecontrol = nrow(controlsample);

		modelMat <- model.matrix(model.formula,data);
		sdmod <- apply(modelMat,2,sd, na.rm = TRUE);
		sdmod[1] <- 1.0;
		sdmod <- as.vector(1.0/sdmod);

		bestMat <- model.matrix(best.model.formula,data);
		if (type=="COX")
		{
			response <- data[,all.vars(model.formula)[1:2]];
			sdmod <- sdmod[-1];
		}
		else
		{
			response <- data[,Outcome];
		}

		output<-.Call("bootstrapValidationBinCpp", fraction, loops, modelMat,type, data.matrix(response),bestMat);

#		colnames(output$bcoef) <- names(basemodel$coefficients);
		basepredict <- predict.fitFRESA(basemodel,testData=data,predictType = 'linear');
		basemodel$linear.predictors <- basepredict;
		framesize = nrow(data);
		acc = 0.0;
		sen = 0.0;
		spe = 0.0;
		for (i in 1:framesize)
		{
			if ((data [i,Outcome] > 0) && (basepredict[i] > 0) ) 
			{
				acc = acc + 1.0; 
				sen = sen + 1.0;
			}
			if ((data [i,Outcome] == 0) && (basepredict[i] < 0) ) 
			{
				acc = acc + 1.0; 
				spe = spe + 1.0;
			}
		}
		baseAcc = acc/framesize;
		baseSen = sen/sizecases;
		baseSpe = spe/sizecontrol;
		bootmodel <- basemodel;
		environment(bootmodel$formula) <- NULL
		environment(bootmodel$terms) <- NULL
#		bootmodel$coefficients[epts] <-  apply(output$bcoef,2,median, na.rm = TRUE);
#		print(sdmod)
#		cat(nrow(output$bcoef),":",ncol(output$bcoef),"\n");
		wts <- 2.0/(1.0+(abs(output$bcoef) %*% sdmod));
#		print(wts)
		bootmodel$coefficients[epts] <-  apply(output$bcoef,2,weighted.mean,w=wts, na.rm = TRUE);
		bootmodel$estimations[epts] <-  bootmodel$coefficients[epts];
#		print(apply(output$bcoef,2,weighted.mean, na.rm = TRUE));
#		print(bootmodel$coefficients[epts])

		pr <- predict.fitFRESA(bootmodel,testData=data,predictType = 'linear');
		bootmodel$linear.predictors <- pr;
		p <- predict.fitFRESA(bootmodel,testData=data,predictType = 'prob');
		bootmodel$fitted.values <- p;
		sen <- sum( 1*((data[,Outcome] > 0)*( bootmodel$linear.predictors >= 0 )) , na.rm = TRUE)
		spe <- sum( 1*((data[,Outcome] == 0)*( bootmodel$linear.predictors < 0 )) , na.rm = TRUE)
		acc <- sen+spe;
		acc = acc/framesize;
		sen = sen/sizecases;
		spe = spe/sizecontrol;
		blidRoc = 0.5;
		bootRoc = NULL;
		if (plots && (length(output$testoutcome)>0))
		{
			par(mfrow=c(2,2))
			pROC::roc( as.vector(data[,Outcome]), basemodel$linear.predictors, col="blue",plot=TRUE,smooth=FALSE,progress= 'none',quiet = TRUE);
			par(new=TRUE)
			blidRoc <- pROC::roc(as.vector(output$testoutcome),output$testprediction,col="red",auc=TRUE,print.auc=TRUE,plot=TRUE,smooth=FALSE,progress= 'none',quiet = TRUE)
			par(new=TRUE)
			bootRoc <- pROC::roc( as.vector(data[,Outcome]), bootmodel$linear.predictors,plot=TRUE,ci=plots,auc=TRUE,of='se',specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),boot.n=200,smooth=FALSE,progress= 'none',quiet = TRUE);
				plot(ecdf(output$taccuracy),main="Accuracy",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindAccuracy,col = "red");
				abline(v=acc,col = "blue");
				plot(ecdf(output$tsensitivity),main="Sensitivity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindSensitivity,col = "red");
				abline(v=sen,col = "blue");
				plot(ecdf(output$tspecificity),main="Specificity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=output$BlindSpecificity,col = "red");
				abline(v=spe,col = "blue");
			par(mfrow=c(1,1))
		}
		else
		{
			blidRoc <- llist(auc=output$resul$blindAUC);
			if ((sum(response) == 0) || sum(response) == nrow(data))
			{
			  warning ("Constant Outcome");
			}
			else
			{
				if (length(pr) == nrow(data))
				{
					pr[is.na(pr)] <- 0;
					pr[pr == -Inf] <- -100000.0;
					pr[pr == Inf] <- 100000.0;
					bootRoc <- pROC::roc( as.vector(data[,Outcome]), pr,plot=FALSE,auc=TRUE,smooth=FALSE,progress= 'none',quiet = TRUE);
				}
				else
				{
					warning ("No prediction");
					bootRoc <- pROC::roc( as.vector(data[,Outcome]), numeric(nrow(data)),plot=FALSE,auc=TRUE,smooth=FALSE,progress= 'none',quiet = TRUE);
				}
			}
		}
		colnames(output$IDI) <- attr(terms(model.formula),'term.labels');
		colnames(output$zIDI) <- attr(terms(model.formula),'term.labels');
		colnames(output$NRI) <- attr(terms(model.formula),'term.labels');
		colnames(output$zNRI) <- attr(terms(model.formula),'term.labels');
		colnames(output$bcoef) <- names(basemodel$coefficients)[epts];

		tem=output$resul;

		result <- structure(llist(
		data=data,
		outcome=data[,Outcome],
		blind.accuracy=output$BlindAccuracy,
		blind.sensitivity=output$BlindSensitivity,
		blind.specificity=output$BlindSpecificity,
		train.ROCAUC=output$trainRoc,
		blind.ROCAUC= blidRoc,
		boot.ROCAUC=bootRoc, 
		fraction=fraction,
		loops=loops,
		base.Accuracy=baseAcc,
		base.sensitivity=baseSen,
		base.specificity=baseSpe,
		accuracy=output$accuracy,
		sensitivity=output$sensitivity,
		specificity=output$specificity,
		train.accuracy=output$taccuracy,
		train.sensitivity=output$tsensitivity,
		train.specificity=output$tspecificity,
		s.coef=output$bcoef,
		boot.model=bootmodel,
		boot.accuracy=acc,
		boot.sensitivity=sen,
		boot.specificity=spe,
		z.NRIs=output$zNRI,
		z.IDIs=output$zIDI,
		test.z.NRIs=tem$test_zNRI,
		test.z.IDIs=tem$test_zIDI,
		startBetas=tem$rawbetas,
		NRIs=output$NRI,
		IDIs=output$IDI,
		testOutcome=output$testoutcome,
		testPrediction=output$testprediction,
		labels = FALSE), class = "bootstrapValidation_Bin");
	}
	else
	{
		result <- structure(llist(
			data=data,
			outcome=data[,Outcome],
			boot.model=bootmodel
		),class = "bootstrapValidation_Bin_error")			
	}
	return (result);
}

