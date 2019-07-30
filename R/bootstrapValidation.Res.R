bootstrapValidation_Res <-
function(fraction=1.00,loops=200,model.formula,Outcome,data,type=c("LM","LOGIT","COX"),plots=FALSE,bestmodel.formula=NULL)
{

    type <- match.arg(type);
	basemodel <- modelFitting(model.formula,data,type,TRUE)
	bootmodel <- basemodel;
	if (class(model.formula) != "formula")
	{
		model.formula <- formula(model.formula);
	}
	if (is.null(bestmodel.formula))
	{
		bestmodel.formula <- model.formula;
	}
	else
	{
		if (class(bestmodel.formula) != "formula")
		{
#			cat("NRI with :",bestmodel.formula,"\n");
			bestmodel.formula <- formula(bestmodel.formula);
		}
	}
#	cat("Loops :",loops,"\n");


	varsList <- all.vars(model.formula)
	svar <- 1+1*(type=="COX")
	epts <- c(svar:length(basemodel$coefficients));
	if (length(varsList)>svar)
	{
		outn <- length(table(data[,Outcome]))
#		print(summary(basemodel))
		
		if ( inherits(basemodel, "try-error"))
		{
			print(model.formula);
			if (type == "LM")
			{
				data[,Outcome] = data[,Outcome] + rnorm(nrow(data),0,1e-10);
			}
			else
			{
				cp <- sample(1:nrow(data),0.05*nrow(data)+1,replace=FALSE);
				data[cp,Outcome] = 1*(data[cp,Outcome]==0);
			}
			basemodel <- modelFitting(model.formula,data,type,TRUE)
			if ( inherits(basemodel, "try-error"))
			{
				stop("BootstapValidation: Initial model Error\n")
			}
		}

		modelMat <- model.matrix(model.formula,data);
		bestmodelMat <- model.matrix(bestmodel.formula,data);
		sdmod <- apply(modelMat,2,sd, na.rm = TRUE);
		sdmod[1] <- 1.0;
		sdmod <- as.vector(1.0/sdmod);

		if (type=="COX")
		{
			response <- data[,all.vars(model.formula)[1:2]];
			sdmod <- sdmod[-1];
		}
		else
		{
			response <- data[,Outcome];
		}

		 output<-.Call("bootstrapValidationResCpp", fraction, loops, modelMat,type,data.matrix(response),bestmodelMat);

		modelMat <- NULL;
		bestmodelMat <- NULL;
		bootmodel <- basemodel;
		environment(bootmodel$formula) <- NULL
		environment(bootmodel$terms) <- NULL

		trainRMSE <- sqrt(mean(output$trainResiduals^2));
		bootvar <- mean(output$testResiduals^2);
		testRMSE <- sqrt(bootvar);
		
#		wt <- ((1.0e-10+2.0*output$trainSampledRMSE)/(1.0e-10+output$trainSampledRMSE+output$testSampledRMSE)); # the pre weights an F variance ratio
		wts <- 2.0/(1.0+(abs(output$bcoef) %*% sdmod));
#		wts <- wt/colSums(t(abs(output$bcoef)));

		pro <- predict.fitFRESA(bootmodel,testData=data,predictType = 'linear');
		bootResRMSE<- sqrt(mean((pro-data[,Outcome])^2))


		bootmodel$coefficients[epts] <-  apply(output$bcoef,2,weighted.mean,w = wts, na.rm = TRUE)
		bootmodel$estimations[epts] <-  bootmodel$coefficients[epts]
#		print(apply(output$bcoef,2,weighted.mean, na.rm = TRUE));
#		print(bootmodel$coefficients[epts])

		pr <- predict.fitFRESA(bootmodel,testData=data,predictType = 'linear');
		bootmodel$linear.predictors <- pr;
		bootResRMSE<- sqrt(mean((pr-data[,Outcome])^2))


		colnames(output$NeRi) <- attr(terms(model.formula),'term.labels');
		colnames(output$tpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$wpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$spvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$Fpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_tpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_wpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_spvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_Fpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$bcoef) <- names(basemodel$coefficients)[epts];


		if (plots)
		{
			if (outn>2) 
			{	
				par(mfrow=c(1,3))
				plot(output$testPrediction ~ output$testOutcome,main="Prediction~Outcome");
				plot(output$testResiduals ~ output$testOutcome,main="Residuals");
			}
			else
			{
				par(mfrow=c(2,2))
				pROC::roc(as.vector(output$testOutcome), output$testPrediction,auc=TRUE,plot=TRUE,smooth=FALSE,print.auc=TRUE,quiet = TRUE);
				boxplot(output$testPrediction ~ output$testOutcome,main="Prediction~Outcome");
				boxplot(output$testResiduals ~ output$testOutcome,main="Residuals");
			}
			plot(ecdf(output$trainSampledRMSE),main="RMSE",col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
			abline(v=testRMSE,col = "red");
			par(mfrow=c(1,1))
			
		}
#		cat("End :",loops,"\n");

		result <- structure(llist( data = data,
			outcome=data[,Outcome],
			s.coef=output$bcoef,
			boot.model=bootmodel,
			NeRis = output$NeRi,
			tStudent.pvalues = output$tpvalue,
			wilcox.pvalues = output$wpvalue,
			bin.pvalues = output$spvalue,
			F.pvalues = output$Fpvalue,
			test.tStudent.pvalues = output$test_tpvalue,
			test.wilcox.pvalues = output$test_wpvalue,
			test.bin.pvalues = output$test_spvalue,
			test.F.pvalues = output$test_Fpvalue,
			testPrediction=output$testPrediction,
			testOutcome=output$testOutcome,
			testResiduals=output$testResiduals,
			testRMSE=testRMSE,
			trainRMSE=trainRMSE,
			trainSampledRMSE=output$trainSampledRMSE,
			testSampledRMSE=output$testSampledRMSE,
			outcomeSD=output$OutcomeSD
		), class = "bootstrapValidation_Res");
	}
	else
	{
		result <- structure(llist(
			data=data,
			outcome=data[,Outcome],
			boot.model=bootmodel
		),class = "bootstrapValidation_Res_error")			
	}
	return (result);
}
