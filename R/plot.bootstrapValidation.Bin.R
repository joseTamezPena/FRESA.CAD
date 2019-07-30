#' @method plot bootstrapValidation_Bin
plot.bootstrapValidation_Bin <-
function(x,xlab = "Years", ylab="Survival",strata.levels=c(0),main="ROC",cex=1.0,...) 
{

	op <- par(no.readonly=TRUE)


	par(mfrow=c(1,1),pty='m')
	classlen=length(class(x$boot.model))
	
	cobj <- substr(class(x$boot.model)[classlen], 1, 2);
	if (x$boot.model$type=="COX")
	{	
		cobj = "co"
	}
	lrp = NULL;
	switch(cobj,
		co =
		{
			varsList <- as.list(attr(terms(x$boot.model),"variables"))
			if (length(strata.levels)==1)
			{
				form = paste(varsList[2]," ~ strata(x$boot.model$linear.predictors > strata.levels)")
				cat(form)
				lrp <- EmpiricalSurvDiff(x$boot.model$model[,1],x$boot.model$model[,2],1*(x$boot.model$linear.predictors > strata.levels),plots=TRUE,main="SLR Null")
				plot(survival::survfit(formula(form),data = x$data),xlab=xlab,ylab=ylab,col=c("blue","red"),conf.int=TRUE,lty = 2:3,...)
				legend(0.1, 0.3, c("Low Risk", "High Risk"), col=c("blue","red"), lty = 2:3)
			}
			else
			{
				strat <- quantile(x$boot.model$linear.predictors,  probs = strata.levels)
				groups <- 1*(x$boot.model$linear.predictors < strat[1])+2*((x$boot.model$linear.predictors >= strat[1]) & (x$boot.model$linear.predictors <= strat[2]) ) + 3*(x$boot.model$linear.predictors > strat[2])
				form = paste(varsList[2]," ~ strata(groups)")
				cat(form)
				plot(survival::survfit(formula(form),data = x$data),xlab=xlab,ylab=ylab,col=c("blue","green","red"),conf.int=TRUE,lty = 2:4,...)
				legend(0.1, 0.3, c("Low Risk", "Mid Risk", "High Risk"), col=c("blue","green","red"), lty = 2:4)
			}
			title("Kaplan-Meier Curve");
		}
	)
	par(mfrow=c(3,1))
	plot(ecdf(x$train.accuracy),main="Accuracy",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
	abline(v=x$blind.accuracy,col = "red");
	abline(v=x$boot.accuracy,col = "blue");
	legend("topleft", legend=c("CDF","Model", "CV"),
       col=c("black","blue", "red"), lwd=c(2,1,1))
	plot(ecdf(x$train.sensitivity),main="Sensitivity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
	abline(v=x$blind.sensitivity,col = "red");
	abline(v=x$boot.sensitivity,col = "blue");
	legend("topleft", legend=c("CDF","Model", "CV"),
       col=c("black","blue", "red"), lwd=c(2,1,1))
	plot(ecdf(x$train.specificity),main="Specificity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
	abline(v=x$blind.specificity,col = "red");
	abline(v=x$boot.specificity,col = "blue");
	legend("topleft", legend=c("CDF","Model", "CV"),
       col=c("black","blue", "red"), lwd=c(2,1,1))
	par(mfrow=c(1,1))	
	par(pty='s');

	prp <- pROC::roc( as.vector(x$outcome), x$boot.model$linear.predictors,progress= 'none',quiet = TRUE);
	thres=0;
	if (min(x$boot.model$linear.predictors)>=0) thres=0.5;
	dtable <- table(x$boot.model$linear.predictors<thres,1-x$outcome);
	

	ci.sp.obj <- ci.sp(prp , sensitivities=seq(0, 1, .05), boot.n=100,progress= 'none',quiet = TRUE)
	ci.se.obj <- ci.se(prp , specificities=seq(0, 1, .05), boot.n=100,progress= 'none',quiet = TRUE)
#	ci.sp.obj <- ci.sp(prp , sensitivities=seq(0, 1, .05), boot.n=100)
#	ci.se.obj <- ci.se(prp , specificities=seq(0, 1, .05), boot.n=100)
	plot(prp,grid=c(0.1, 0.1),grid.col=c("gray", "gray")) # restart a new plot
	plot(ci.sp.obj, type="s", col="gray")
	plot(ci.se.obj, type="s", col="light gray")

	par(new=TRUE)
	pr <- pROC::roc(as.vector(x$testOutcome),x$testPrediction,col="red",auc=TRUE,print.auc=FALSE,plot=TRUE,smooth=FALSE,main=main,progress= 'none',quiet = TRUE)
#	legend("bottomright", legend=c("Predictor", "Bootstrap CV"),
 #      col=c("black", "red"), lwd=2)

	if (nrow(dtable)==ncol(dtable) & (ncol(dtable)>1))
	{
		colnames(dtable) <- c("O(+)","O(-)")
		rownames(dtable) <- c("T(+)","T(-)")
		Sen=dtable[1,1]/(dtable[1,1]+dtable[2,1])
		Spe=dtable[2,2]/(dtable[1,2]+dtable[2,2])
		enauc = 0.5*(Sen+Spe)	

		lines(c(1,Spe,0),c(0,Sen,1),col="green",lwd=1.0,lty=1);
    
		ley.names <- c(paste("Predictor (",sprintf("%.3f",prp$auc),")"),paste("Bootstrap CV (",sprintf("%.3f",pr$auc),")"),paste(" Test Prediction (",sprintf("%.3f",enauc),")"));
		ley.colors <- c("black", "red","green");
		ley.lty <- c(1,1,1);

		legend(0.7,0.20, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n",cex=0.9*cex)
		Acc=(dtable[1,1]+dtable[2,2])/(dtable[1,1]+dtable[2,2]+dtable[1,2]+dtable[2,1])
		F1=(2*dtable[1,1])/(2*dtable[1,1]+dtable[1,2]+dtable[2,1])
		x=1.025
		y=1.025
		text(x,y,paste("(TPR=",sprintf("%.3f",Sen),",TNR=",sprintf("%.3f",Spe),",ACC=",sprintf("%.3f",Acc),",F1=",sprintf("%.3f",F1),",AUC=",sprintf("%.3f",enauc),")"),adj = c(0,1),cex=0.7*cex,col="dark green")
		par(new=TRUE,plt=c(0.6,0.8,0.37,0.57),pty='s',cex=0.9*cex)
		plot(t(dtable),main="Confusion Matrix",ylab="Test",xlab="Outcome",cex=0.9*cex)
	}
	else
	{
		legend("bottomright", legend=c("Predictor", "Bootstrap CV"),
		col=c("black", "red"), lwd=2)
	}


	par(op)

	result <- list(proc.test=pr,proc.prediction=prp,diagnosticMatrix=dtable,LRP=lrp)
	return (result)
}