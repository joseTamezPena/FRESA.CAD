#' @method summary bootstrapValidation_Bin

summary.bootstrapValidation_Bin <- 
function(object, ...) 
{


	 digits = 3 
	conf.int = 0.95
	cilow = (1.0-conf.int)/2;
	cihig = 1.0-cilow;

	ciAccuracy <- as.vector(quantile(object$train.accuracy, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	ciSensitivity <- as.vector(quantile(object$train.sensitivity, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	ciSpecificity <- as.vector(quantile(object$train.specificity, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));

	citROC <- as.vector(quantile(object$train.ROCAUC, probs = c(cilow,0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
	
		
	bootAuc <- pROC::roc( as.vector(object$outcome), object$boot.model$linear.predictors,plot=FALSE,ci=TRUE,boot.n=object$loops,quiet = TRUE);
	
    cat("\nModel Cross-Validation with Improvement in Predicted Probability\n\n")
    cat("Number of Cases:", sum(object$outcome), "\t Number of Controls", sum(object$outcome==0), "\n\n")
    cat("Number of Bootstraps:", length(object$train.ROCAUC), "\t Sampled Fraction", object$fraction, "\n\n")

	performance <- vector();
	
	cat(l1 <- sprintf("Blind    Accuracy: %8.3f : Bootstrapped    Accuracy: %8.3f (%8.3f to %8.3f) \n",object$blind.accuracy,ciAccuracy[2],ciAccuracy[1], ciAccuracy[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind Sensitivity: %8.3f : Bootstrapped Sensitivity: %8.3f (%8.3f to %8.3f) \n",object$blind.sensitivity,ciSensitivity[2],ciSensitivity[1], ciSensitivity[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind Specificity: %8.3f : Bootstrapped Specificity: %8.3f (%8.3f to %8.3f) \n",object$blind.specificity,ciSpecificity[2],ciSpecificity[1], ciSpecificity[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind     ROC AUC: %8.3f : Bootstrapped     ROC AUC: %8.3f (%8.3f to %8.3f) \n",object$blind.ROCAUC$auc,object$boot.ROCAUC$auc,bootAuc$ci[1],bootAuc$ci[3]));
	performance <- append(performance,l1);
	cat(l1 <- sprintf("Blind     ROC AUC: %8.3f : ModelBootstrap   ROC AUC: %8.3f (%8.3f to %8.3f) \n\n",object$blind.ROCAUC$auc,citROC[2],citROC[1],citROC[3]));
	performance <- append(performance,l1);
	
	performance.table <- rbind(c(object$blind.accuracy,ciAccuracy[2],ciAccuracy[1], ciAccuracy[3]));
	performance.table <- rbind(performance.table,c(object$blind.sensitivity,ciSensitivity[2],ciSensitivity[1], ciSensitivity[3]));
	performance.table <- rbind(performance.table,c(object$blind.specificity,ciSpecificity[2],ciSpecificity[1], ciSpecificity[3]));
	performance.table <- rbind(performance.table,c(object$blind.ROCAUC$auc,object$boot.ROCAUC$auc,bootAuc$ci[1],bootAuc$ci[3]));
	performance.table <- rbind(performance.table,c(object$blind.ROCAUC$auc,citROC[2],citROC[1],citROC[3]));
	colnames(performance.table) <- c("Blind","Train","LCI","UCI")
	rownames(performance.table) <- c("Accuracy","Sensitivity","Specificity","ROCAUC 1","ROCAUC 2")
	
	
#	smry <- summary(object$boot.model, ...);
	meancoef <- colMeans(object$s.coef,na.rm = TRUE);
	
	ncoef <- length(meancoef);
	lowci <- numeric(ncoef);
	topci <- numeric(ncoef);
	zidimedian <- numeric(ncoef);
	idimedian <- numeric(ncoef);
	idilowci <- numeric(ncoef);
	iditopci <- numeric(ncoef);
	znrimedian <- numeric(ncoef);
	nrimedian <- numeric(ncoef);
	nrilowci <- numeric(ncoef);
	nritopci <- numeric(ncoef);
	soffset <- ncoef-ncol(object$IDIs);
	for ( i in 1:ncoef)
	{
		ci <- as.vector(quantile(object$s.coef[,i], probs = c(cilow, cihig), na.rm = TRUE,names = FALSE, type = 7));
		lowci[i] <- ci[1];
		topci[i] <- ci[2];
		j <- i-soffset;
		if (j>0)
		{
			ci <- as.vector(quantile(object$IDIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
			idimedian[i] <- ci[2];
			idilowci[i] <- ci[1];
			iditopci[i] <- ci[3];		
			ci <- as.vector(quantile(object$z.IDIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
			zidimedian[i] <- ci[2];

			ci <- as.vector(quantile(object$NRIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
			nrimedian[i] <- ci[2];
			nrilowci[i] <- ci[1];
			nritopci[i] <- ci[3];		
			ci <- as.vector(quantile(object$z.NRIs[,j], probs = c(cilow, 0.5, cihig), na.rm = TRUE,names = FALSE, type = 7));
			znrimedian[i] <- ci[2];
		}
	}

	p <- meancoef;
	p <- cbind(p,lowci);
	p <- cbind(p,topci);
	p <- cbind(p,idimedian);
	p <- cbind(p,idilowci);
	p <- cbind(p,iditopci);
	p <- cbind(p,zidimedian);

	p <- cbind(p,nrimedian);
	p <- cbind(p,nrilowci);
	p <- cbind(p,nritopci);
	p <- cbind(p,znrimedian);
	
	print(p)
	print(colnames(object$s.coef))

	cnames <- c("Coef","Low CI","High CI","Median IDI","Low IDI","High IDI","z IDI","Median NRI","Low NRI","High NRI","z NRI");
	colnames (p) <- cnames;
	rownames(p) <- colnames(object$s.coef);
	print(p, digits = digits);
	
	result <- list(performance = performance,
	coef = p, 
	performance.table = performance.table);
	
	return (result);
        
}
