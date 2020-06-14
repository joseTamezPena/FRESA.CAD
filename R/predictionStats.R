predictionStats_survival <-  function(predictions, plotname="",...)
{
	if(!sum(predictions[,4] == 0) == length(predictions[,4]))
	{
		if (!requireNamespace("survminer", quietly = TRUE)) {
			install.packages("survminer", dependencies = TRUE)
		}

		data <- data.frame(times=predictions[,1],preds=predictions[,5])
		CIFollowUp <- concordance95ci(datatest = data, followUp = TRUE)
		data <- data.frame(times=predictions[,1],preds=-predictions[,6])
		CIRisk <- concordance95ci(datatest = data, followUp = FALSE)
		newData <- data.frame(times=predictions[,1],status=predictions[,2],preds=predictions[,6]);
		Curves <- survival::survfit(survival::Surv(times, status) ~ preds>=median(preds),newData)
		
		function (fit, scale, k = 2, ...) 
		{
			edf <- sum(!is.na(fit$coefficients))
			loglik <- fit$loglik[length(fit$loglik)]
			c(edf, -2 * loglik + k * edf)
		}
		
		#Curves <- survival::survfitBSWiMS <- survival::survfit(Surv(predictions[,1], predictions[,2]) ~ predictions[,6]>=median(predictions[,6]))
		if(plotname!="")
		{
			 graph <- survminer::ggsurvplot(Curves, data=newData, conf.int = TRUE, legend.labs=c("Low Risk", "High Risk"),
                        palette = c("#00bbff", "#ff0000"),
                        ggtheme = ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 30)),
                        title = plotname,
                        risk.table = TRUE,
                        tables.height = 0.2,
                        tables.theme = survminer::theme_cleantable())
      		print(graph)
			# graph <- try(survminer::ggsurvplot(Curves, data=newData, conf.int = TRUE, legend.labs=c("Low Risk", "High Risk"),
			# 					palette = c("#00bbff", "#ff0000"),
			# 					ggtheme = ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 30)),
			# 					title = plotname,
			# 					risk.table = TRUE,
			# 					tables.height = 0.2,
			# 					tables.theme = survminer::theme_cleantable()))
			# # if (!inherits(graph, "try-error"))
			# {
			# 	print(graph)
			# }
			# else{
			# 	warning("Survplot failed");
			# }
			#plot(Curves,main=plotname)
		}
		
		LogRank <- survival::survdiff(Surv(predictions[,1], predictions[,2]) ~ predictions[,6]>=median(predictions[,6]))
		LogRank <- cbind(LogRank$chisq,  1 - pchisq(LogRank$chisq, length(LogRank$n) - 1));
		colnames(LogRank) <- cbind("chisq","pvalue");
		return( list(CIFollowUp=CIFollowUp, CIRisk = CIRisk, LogRank = LogRank, Curves = Curves) );
	}
	else{
		return( list(CIFollowUp=rep(0,nrow(predictions)), CIRisk = rep(0,nrow(predictions)), LogRank = rep(0,nrow(predictions)), Curves = NULL) );
	}
}

concordance95ci <- function(datatest,nss=1000,followUp=FALSE)
{
  sz <- nrow(datatest)
  sesci <- c(0,0,0);
  if (sz>2)
  {
    ses <- numeric(nss);
    for (i in 1:nss)
    {
      bootsample <- datatest[sample(sz,sz,replace=TRUE),];
      if(followUp)
      {
        ses[i] <- rcorr.cens(bootsample[,1], bootsample[,2], outx = TRUE)[1]
      }
      else
      {
        ses[i] <- rcorr.cens(bootsample[,1], bootsample[,2])[1]
      }
    }
    sesci <- quantile(ses, probs = c(0.5,0.025,0.975),na.rm = TRUE);
    sesci[1]<-mean(ses)
    names(sesci)<-c("median","lower","upper");
  }
  return (sesci);
}

predictionStats_ordinal <-  function(predictions,plotname="",...)
{
    cat(plotname,"\n")
	dpoints <- nrow(predictions)
	tint <- qt(0.975,dpoints - 1)/sqrt(dpoints)
    ScoresOutcome <- predictions[,1]
    if (nchar(plotname) > 0) 
    {
      boxplot(predictions[,2] ~ ScoresOutcome,main = plotname,xlab = "Outcome" ,ylab = "Prediction",...)
    }
    ct <-  DescTools::KendallTauB(ScoresOutcome,as.integer(predictions[,2]+0.5),conf.level = 0.95)
    bias <- mean(predictions[,2] - ScoresOutcome)
    rstd <- sqrt(mean((predictions[,2] - ScoresOutcome)^2) - bias^2)
    Bias <- c(bias,bias - tint*rstd,bias + tint*rstd)
	theScores <- as.numeric(names(table(ScoresOutcome)))
	BMAE <- NULL;
	for (s in theScores)
	{
		BMAE <- rbind(BMAE,MAE95ci(predictions[ScoresOutcome==s,])); 
	}
	BMAE <- colMeans(BMAE);
	class95ci <- ClassMetric95ci(predictions);
    kp <- irr::kappa2(cbind(as.integer(predictions[,2] + 0.5),ScoresOutcome),"squared")
    zdis <- 2*kp$value/kp$statistic;
    Kapp <- c(kp$value, kp$value - zdis, kp$value + zdis)
    results <- list(Kendall = ct,
					Bias = Bias,
					BMAE = BMAE,
					Kapp = Kapp,
					class95ci = class95ci,
					KendallTauB = ct,
					Kappa.analysis = kp
					);
    return(results);
}

metric95ci <- function(metric,nss=1000,ssize=0)
{
	sz <- length(metric);
	if (ssize == 0)
	{
		ssize <- sz;
	}
	ssize <- min(sz,ssize);
	metricci <- c(0,0,0);
	if (sz>1)
	{
		meanMetric <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- metric[sample(sz,ssize,replace=TRUE)];
		  meanMetric[i] <- mean(bootsample);
		}
		metricci <- quantile(meanMetric, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (metricci);
}

corcen95ci <- function(dataTable,nss=1000)
{
	sz <- nrow(dataTable);
	metricci <- c(0.5,0.5,0.5);
	if (sz>1)
	{
		meanMetric <- numeric(nss);
		for (i in 1:nss)
		{
		  bootsample <- dataTable[sample(sz,sz,replace=TRUE),];
		  meanMetric[i] <- rcorr.cens(bootsample[,2],bootsample[,1], outx=FALSE)[1];
		}
		metricci <- quantile(meanMetric, probs = c(0.5,0.025, 0.975),na.rm = TRUE);
	}
	return (metricci);
}

predictionStats_binary <-  function(predictions, plotname="", center=FALSE,...)
{
    cat(plotname,"\n")
	cstat <- NULL;
	cstatCI <- c(0.5,0.5,0.5);
	medianTest <- NULL;
	if (ncol(predictions)>2)
	{
		numberOfModels <- table(predictions[,2]);
		numberOfModels <- as.integer(names(numberOfModels));
		cstat <- numeric(length(numberOfModels))
		modSize <- 0;
		for (mi in numberOfModels)
		{
			  mtest <- predictions[,2] == mi;
			  modSize <- modSize+sum(mtest);
			  cstat[mi] <- rcorr.cens(predictions[mtest,3],predictions[mtest,1], outx=FALSE)[1];
		}
		modSize <- modSize/length(numberOfModels);
		boxstaTest <- try(boxplot(as.numeric(as.character(predictions[,3]))~rownames(predictions),plot = FALSE));
		if (!inherits(boxstaTest, "try-error"))
		{
			medianTest <- cbind(predictions[boxstaTest$names,1],boxstaTest$stats[3,]);
			rownames(medianTest) <- boxstaTest$names;
		}
		else
		{
			warning("boxplot test failed");
			medianTest <- cbind(predictions[,1],rep(0,nrow(predictions)));
			rownames(medianTest) <- rownames(predictions);
		}
		colnames(medianTest) <- c("Outcome","Median");
		smpCI <- as.integer(nrow(medianTest)/modSize+0.5);
#		cat("Avg size:", modSize,"Samp size:",smpCI,"\n");
		if (length(cstat) > 1) 
		{
			cstatCI <- metric95ci(cstat,ssize=smpCI);
		}
		else 
		{
			cstatCI <- c(cstat[1],0.0,1.0);
		}
		predictions	<- medianTest;
	}
	else
	{
		cstatCI <- corcen95ci(predictions,200 + 800*(nrow(predictions) < 1000) );
	}
	if (any(is.na(predictions[,2])))
	{	
		predictions <- predictions[!is.na(predictions[,2]),]
	}	
    if (center) predictions[,2] <- predictions[,2] - 0.5;
	if (min(predictions[,2]) >= 0)
	{
		predictions[,2] <- predictions[,2] - 0.5;
	}
	
    pm <- NULL;
	citest <- NULL;
    if (nchar(plotname) > 1)
    {
      pm <- plotModels.ROC(predictions,main = plotname,...);
      cis <- ci.auc(pm$roc.predictor)
    }
    else
    {
      pm <- pROC::roc(as.vector(predictions[,1]),predictions[,2],quiet = TRUE);
      cis <- ci.auc(pm);
      pm$predictionTable <- table(predictions[,2] < 0,1 - predictions[,1]);
		if (nrow(pm$predictionTable) == 1)
		{
			if ((rownames(pm$predictionTable) == "0") || (rownames(pm$predictionTable) == "FALSE"))
			{
				pm$predictionTable <- rbind(c(0,0),pm$predictionTable);
			}
			else
			{
				pm$predictionTable <- rbind(pm$predictionTable,c(0,0));
			}
			rownames(pm$predictionTable) <- c("0","1")
		}
    }
	class95ci <- ClassMetric95ci(cbind(predictions[,1],predictions[,2] >= 0),200 + 800*(nrow(predictions) < 1000) );

#    print(pm$predictionTable)
    if (length(pm$predictionTable) > 2 )
    {
      ci <- epiR::epi.tests(pm$predictionTable);
      accc <- ci$elements$diag.acc;
      berror <- class95ci$berci;
      sensitivity <- ci$elements$sensitivity;
      specificity <- ci$elements$specificity;
    }
    else
    {
      accc <- c(0.5,0.5,0.5);
      cIndexSet <- c(0.5,0.5,0.5);
      cstatCI <- c(0.5,0.5,0.5);
      berror <- c(0.5,0.5,0.5);
      sensitivity <- c(0.0,0.0,0.0);
      specificity <- c(0.0,0.0,0.0);
      names(accc) <- c("est","lower","upper")
      names(berror) <- c("est","lower","upper")
    }
    aucs <- cis[c(2,1,3)];
	names(aucs) <- c("est","lower","upper")
    results <- list(accc = accc,berror = berror,
					aucs = aucs,sensitivity = sensitivity,
					specificity = specificity,
					ROC.analysis = pm,
					CM.analysis = ci,
					ClassMetrics = class95ci,
					cIndexSet = cstat,
					cIndexCI = cstatCI,
					medianTest = medianTest
					);
    return(results);
}

predictionStats_regression <-  function(predictions, plotname="",...)
{
      cat(plotname,"\n")
	  dpoints <- nrow(predictions)
	  chsqup <- sqrt(dpoints/qchisq(0.025, df = dpoints))
	  chsqdown <- sqrt(dpoints/qchisq(0.975, df = dpoints))
	  tint <- qt(0.975,dpoints - 1)/sqrt(dpoints)

	  if (nchar(plotname) > 0) 
      {
         plot(predictions[,2] ~ predictions[,1],main = plotname,xlab = "Outcome", ylab = "Prediction",...)
      }
	  ct <- cor.test(predictions[,1],predictions[,2],method = "pearson");
	  corci <- c(ct$estimate,ct$conf.int);
	  bias <- mean(predictions[,2] - predictions[,1]);
	  rmse <- sqrt(mean((predictions[,2] - predictions[,1])^2));
	  rstd <- sqrt(rmse^2 - bias^2);
	  biasci <- c(bias,bias - tint*rstd,bias + tint*rstd);
	  RMSEci <- c(rmse,chsqdown*rmse,chsqup*rmse);
	  spearmanci <- sperman95ci(predictions);
	  MAEci <- MAE95ci(predictions);
	  results <- list(corci = corci, 
						biasci= biasci,
						RMSEci=RMSEci,
						spearmanci=spearmanci,
						MAEci=MAEci,
						pearson=ct
						);
	  return(results);
}
