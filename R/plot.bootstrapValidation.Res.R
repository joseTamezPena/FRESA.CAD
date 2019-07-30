#' @method plot bootstrapValidation_Res
plot.bootstrapValidation_Res <-
function(x,xlab = "Years", ylab="Survival",...) 
{
  opg <- par(no.readonly = TRUE)

	par(mfrow=c(1,1),pty='m')
	classlen=length(class(x$boot.model))
	
	cobj <- substr(class(x$boot.model)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			varsList <- as.list(attr(terms(x$boot.model),"variables"))
			form = paste(varsList[2]," ~ strata(x$boot.model$linear.predictors >0)")
			cat(form)
			plot(survival::survfit(formula(form),data = x$data),xlab=xlab,ylab=ylab,col=c("blue","red"),conf.int=TRUE,lty = 2:3,...)
			legend(0.1, 0.3, c("Low Risk", "High Risk"), col=c("blue","red"), lty = 2:3)
			title("Kaplan-Meier Curve");
		},
		{
#			par(mfrow=c(2,2))
#			try(plot(x$boot.model,ask=FALSE,...));
#			par(mfrow=c(1,1))
		}
	)
	par(mfrow=c(3,1))
	plot(x$testOutcome ~ x$testPrediction,main="Outcome ~ Test Presiction",xlab="Test Prediction",ylab="Outcome");
	plot(x$testOutcome ~ x$testResiduals,main="Test Residuals",xlab="Test Residuals",ylab="Outcome");
	plot(ecdf(x$trainSampledRMSE),main="RMSE",col="black",lwd = 2,verticals = TRUE, do.points = FALSE,xlab="RMSE",ylab="Probability");
	abline(v=x$testRMSE,col = "red");
	legend("topleft", legend=c("Train CDF","Test RMSE"),
       col=c("black", "red"), lwd=2)
	par(opg)
}
