#' @method update fitFRESA
summary.fitFRESA <- function(object,type=c("Improvement","Residual"),ci=c(0.025,0.975),data=NULL,...)
{
	coefficients <- NULL;
	type <- match.arg(type);
	coffset <- 1.0*(object$type == "COX") 
	outcome <- as.character(all.vars(object$formula)[1+coffset]);
	FRESAsummary <- NULL;
	bv <- NULL;
	classlen=length(class(object))
	cobj <- substr(class(object)[classlen], 1, 2);
	uthr <- 0;
	if (cobj == "or")
	{
		FRESAsummary <- list();
		for (s in 1:length(object$theScores[-1]))
		{
			FRESAsummary[[s]] <- summary.fitFRESA(object$theBaggedModels[[s]]$bagged.model);
		}
		return (FRESAsummary);
	}
	if (cobj == "BS")
	{
		if (is.null(object$oridinalModels))
		{
			return (summary.fitFRESA(object$bagging$bagged.model));
		}
		else
		{
			return (summary.fitFRESA(object$oridinalModels));
		}
	}
	if (length(object$coefficients)>1)
	{
		zval1 <- qnorm(ci[1]);
		zval2 <- qnorm(ci[2]);
		if (is.null(data)) 
		{
			data <- object$model;
		}
		if ((object$type=="LM") || (type=="Residual"))
		{
			coefficients <- object$coefficients[-1];
			ncoef <- length(coefficients);
			cis <- matrix(rep(NA,ncoef*3),ncoef,3);
			cinames <- c("lower","mean","upper");
			if (hasName(object,"baggingAnalysis"))
			{
				vres <- list();
				vres$unitrainMSE <- object$baggingAnalysis$uMS_values;
				vres$redtrainMSE <- object$baggingAnalysis$rMS_values;
				vres$NeRIs <- object$baggingAnalysis$NeRI_values;
				vres$tP.value <- object$baggingAnalysis$pt_values;
				vres$WilcoxP.value <- object$baggingAnalysis$pWilcox_values;
				vres$FP.value <- object$baggingAnalysis$pF_values;
				vres$BinP.value <- object$baggingAnalysis$pBin_values;
				vres$FullTrainMSE <- object$baggingAnalysis$mMSE_values;
				vres$RelativeFrequency <- as.vector(object$baggingAnalysis$RelativeFrequency[names(coefficients)]);
#				print(vres$RelativeFrequency)
#				print(coefficients);
				lowci <- object$baggingAnalysis$coefficients;
				hici <- object$baggingAnalysis$coefficients;
				for (n in 1:length(object$baggingAnalysis$coefficients))
				{
					# if (object$baggingAnalysis$coeff_n_samples[n]>10)
					# {
						# tvalue1 <- qt(ci[1],object$baggingAnalysis$coeff_n_samples[n]-1);
						# tvalue2 <- qt(ci[2],object$baggingAnalysis$coeff_n_samples[n]-1);
						# lowci[n] <- tvalue1*object$baggingAnalysis$coefstd[n];
						# hici[n] <- tvalue2*object$baggingAnalysis$coefstd[n];
						# lowci[n] <- max(lowci[n],zval1*abs(object$baggingAnalysis$coefficients[n]/object$baggingAnalysis$avgZvalues[n]));
						# hici[n] <- max(hici[n],zval2*abs(object$baggingAnalysis$coefficients[n]/object$baggingAnalysis$avgZvalues[n]));
					# }
					# else
					{
						azval <- qnorm(exp(-object$baggingAnalysis$avgLogPvalues[n]))
						lowci[n] <- zval1*abs(lowci[n]/azval);
						hici[n] <- zval2*abs(hici[n]/azval);
					}
				}
				cis = cbind(object$baggingAnalysis$coefficients+lowci,object$baggingAnalysis$coefficients,object$baggingAnalysis$coefficients+hici);
			}
			else
			{
				vres <- getVar.Res(object,data=data,Outcome=outcome,type=object$type);
				bv <- bootstrapValidation_Res(model.formula=object$formula,Outcome=outcome,data=data,type=object$type,...)
				coffset <- 1.0*(object$type != "COX") 
				for (i in 1:ncoef)
				{
					cis[i,] = as.vector(quantile(bv$s.coef[,i+coffset], probs = c(ci[1],0.5,ci[2]), na.rm = TRUE,names = FALSE, type = 7));
				}
				cinames <- c("lower","median","upper");
				vres$RelativeFrequency <- rep(1,ncoef);
			}
			coefficients=as.data.frame(cbind(coefficients,cis,vres$unitrainMSE,vres$redtrainMSE,vres$FullTrainMSE,vres$NeRIs,vres$FP.value,vres$tP.value,vres$BinP.value,vres$WilcoxP.value,vres$RelativeFrequency));
			colnames(coefficients) <- c("Estimate",cinames,"u.MSE","r.MSE","model.MSE","NeRI","F.pvalue","t.pvalue","Sign.pvalue","Wilcox.pvalue","Frequency");
			coefficients <- coefficients[order(vres$FP.value),];

			residaulsMSE <- mean((object$response[,2]-predict(object))^2);
			Rsquare <- var(object$response[,2]);
			Rsquare <- (Rsquare-residaulsMSE)/Rsquare;
			FRESAsummary <- list(coefficients=coefficients,MSE=residaulsMSE,R2=Rsquare,bootstrap=bv);
		}
		else
		{
			coefficients <- object$coefficients[-1];
			ncoef <- length(coefficients);
			cis <- matrix(rep(NA,ncoef*3),ncoef,3);
			cinames <- c("lower","mean","upper");
			if (hasName(object,"baggingAnalysis"))
			{
				vres <- list();
				uthr <- 0.5;
				vres$uniTrainAccuracy <- object$baggingAnalysis$uAcc_values;
				vres$redtrainAccuracy <- object$baggingAnalysis$rAcc_values;
				vres$uniTrainAUC <- object$baggingAnalysis$uAUC_values;
				vres$redtrainAUC <- object$baggingAnalysis$rAUC_values;
				vres$IDIs <- object$baggingAnalysis$idi_values;
				vres$NRIs <- object$baggingAnalysis$nri_values;
				vres$z.IDIs <- object$baggingAnalysis$zidi_values;
				vres$z.NRIs <- object$baggingAnalysis$znri_values;
				vres$fullTrainAccuracy <- object$baggingAnalysis$mACC_values;
				vres$fullTrainAUC <- object$baggingAnalysis$mAUC_values;
				vres$RelativeFrequency <- as.vector(object$baggingAnalysis$RelativeFrequency[names(coefficients)]);
#				print(vres$RelativeFrequency);
#				print(coefficients);

				lowci <- object$baggingAnalysis$coefficients;
				hici <- object$baggingAnalysis$coefficients;
				for (n in 1:length(object$baggingAnalysis$coefficients))
				{
					# if (object$baggingAnalysis$coeff_n_samples[n]>10)
					# {
						# tvalue1 <- qt(ci[1],object$baggingAnalysis$coeff_n_samples[n]-1);
						# tvalue2 <- qt(ci[2],object$baggingAnalysis$coeff_n_samples[n]-1);
						# lowci[n] <- tvalue1*object$baggingAnalysis$coefstd[n];
						# hici[n] <- tvalue2*object$baggingAnalysis$coefstd[n];
						# lowci[n] <- max(lowci[n],zval1*abs(object$baggingAnalysis$coefficients[n]/object$baggingAnalysis$avgZvalues[n]));
						# hici[n] <- max(hici[n],zval2*abs(object$baggingAnalysis$coefficients[n]/object$baggingAnalysis$avgZvalues[n]));
					# }
					# else
					{
						azval <- qnorm(exp(-object$baggingAnalysis$avgLogPvalues[n]))
						lowci[n] <- zval1*abs(lowci[n]/azval);
						hici[n] <- zval2*abs(hici[n]/azval);
					}
				}
				cis = cbind(object$baggingAnalysis$coefficients+lowci,object$baggingAnalysis$coefficients,object$baggingAnalysis$coefficients+hici);
				cis = exp(cis);
				cinames <- c("lower","OR","upper");
			}
			else
			{
				vres <- getVar.Bin(object,data=data,Outcome=outcome,type=object$type);
				bv <- bootstrapValidation_Bin(model.formula=object$formula,Outcome=outcome,data=data,type=object$type,...)
				coffset <- 1.0*(object$type != "COX") 
				for (i in 1:ncoef)
				{
					cis[i,] = as.vector(quantile(bv$s.coef[,i+coffset], probs = c(ci[1],0.5,ci[2]), na.rm = TRUE,names = FALSE, type = 7));
				}
				vres$RelativeFrequency <- rep(1,ncoef);
			}
			coefficients=as.data.frame(cbind(coefficients,cis,vres$uniTrainAccuracy,vres$redtrainAccuracy,vres$fullTrainAccuracy,vres$uniTrainAUC,vres$redtrainAUC,vres$fullTrainAUC,vres$IDIs,vres$NRIs,vres$z.IDIs,vres$z.NRIs,vres$RelativeFrequency));
			colnames(coefficients) <- c("Estimate",cinames,"u.Accuracy","r.Accuracy","full.Accuracy","u.AUC","r.AUC","full.AUC","IDI","NRI","z.IDI","z.NRI","Frequency");
			coefficients <- coefficients[order(-vres$z.IDIs),];
			
			pred <- 1*(predict(object)>uthr);
			Accuracy <- mean(1.0*(object$response[,2]==pred));
			sensitivity <- sum((object$response[,2]==pred)*object$response[,2])/sum(object$response[,2]); 
			specificity <- sum((object$response[,2]==pred)*(object$response[,2]==0))/sum(object$response[,2]==0);
			tAUC <- 0.5*(sensitivity+specificity);
			FRESAsummary <- list(coefficients=coefficients,Accuracy=Accuracy,tAUC=tAUC,sensitivity=sensitivity,specificity=specificity,bootstrap=bv);
		}
	}
	return (FRESAsummary);
}
