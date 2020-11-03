featureAdjustment <-
function(variableList,baseModel,strata=NA,data,referenceframe,type=c("LM","GLS","RLM","RQNT","SPLINE"),pvalue=0.05,correlationGroup = "ID") 
{

if (!requireNamespace("nlme", quietly = TRUE)) {
   install.packages("nlme", dependencies = TRUE)
} 
if (!requireNamespace("MASS", quietly = TRUE)) {
	install.packages("MASS", dependencies = TRUE)
}
	type <- match.arg(type);
	## the reference frame will be used to predict a variable from the basemodel. At output the residuals are returned.
	## strata is a numeric column varname in the data frame from 0 to S, where S is the maximum number of strata
	colnamesList <- as.vector(variableList[,1]);
	size = length(colnamesList);
	if (!is.na(strata)) 
	{
		maxStrata = max(referenceframe[,strata],na.rm = TRUE);
		minStrata = min(referenceframe[,strata],na.rm = TRUE);
	}
	else 
	{
		maxStrata=1;
		minStrata=1;
	}
#	cat ("Min Strata:",minStrata,"Max Strata:",maxStrata,"\n");

#	SortedCtr = referenceframe;
#	nrowsCtr =  nrow(referenceframe);
	created = 0;
	models <- list();
	idx <- 1; 
	tbbaseModel <- NULL;
	for (sta in minStrata:maxStrata)
	{
		if (!is.na(strata))
		{
			stracondition = paste (strata,paste('==',sta));
			strastatement = paste ("subset(referenceframe,",paste(stracondition,")"));
#			cat ("Strata:",stracondition,"\n");
			cstrataref <- eval(parse(text=strastatement));
			strastatement = paste ("subset(data,",paste(stracondition,")"));
			cstrata <- eval(parse(text=strastatement));
#			cat ("Rows:",nrow(cstrataref),"Rows 2",nrow(cstrata)," \n");
		}
		else
		{
			cstrataref = referenceframe;
			cstrata = data;
		}
		if ((nrow(cstrata)>1) && ( nrow(cstrataref)>1))
		{
			if (type=="SPLINE")
			{
				tbbaseModel <- table(cstrataref[,baseModel])
			}
			for (i in 1:size)
			{ 
				tb <- table(cstrataref[,colnamesList[i]])
				if (length(tb) > 4)
				{
					ftm1 <- paste(colnamesList[i],paste(" ~ ",baseModel));
					ftmp <- formula(ftm1);

					switch(type,
    					SPLINE =
						{
#							cat(baseModel,": Variable:\t ",colnamesList[i]);
							rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)
							if (length(tbbaseModel) > 9)
							{
								model <- smooth.spline(cstrataref[,baseModel], 
											  y = cstrataref[,colnamesList[i]],
											  nknots = 5,
#											  df = 3,
											  keep.data = FALSE)						
								dgf = nrow(cstrataref) - model$df;
								pred <- predict(model,cstrataref[,baseModel]);
								ress <- cstrataref[,colnamesList[i]] - pred$y;
								rss2 <- var(ress,na.rm = TRUE);
							}
							else
							{
								model <- lm(ftmp,data=cstrataref,na.action=na.exclude);						
								dgf = nrow(cstrataref) - 1;
								rss2 <- var(model$residuals,na.rm = TRUE);
							}
							f1 = rss1/rss2;
							p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
#							cat(colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},						
						RQNT = 
						{ 
							model <- lm(ftmp,data=cstrataref,na.action=na.exclude);
							iqrref <- as.vector(quantile(cstrataref[,baseModel], probs = c(0.1, 0.5, 0.9), na.rm = TRUE,names = FALSE, type = 7));
							medref <- iqrref[2];
							iqrref <- iqrref[3] - iqrref[1];
							if (iqrref == 0) iqrref <- sd(cstrataref[,baseModel],na.rm = TRUE);
							if (iqrref == 0) iqrref <- 1.0;
							iqrtarget <- as.vector(quantile(cstrataref[,colnamesList[i]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE,names = FALSE, type = 7));
							medtarget <- iqrtarget[2];
							iqrtarget <- iqrtarget[3] - iqrtarget[1];
							if (iqrtarget == 0) iqrtarget <- sd(cstrataref[,colnamesList[i]],na.rm = TRUE);

							slope = iqrtarget/iqrref;

							dpp <- (cstrataref[,colnamesList[i]]-medtarget) - (cstrataref[,baseModel]-medref)*slope;
							iqrdiff <- as.vector(quantile(dpp, probs = c(0.1, 0.5, 0.9), na.rm = TRUE,names = FALSE, type = 7));
							iqrdiff <- iqrdiff[3] - iqrdiff[1];

							model$coefficients[2] <- (iqrtarget-iqrdiff)/iqrref;
							model$coefficients[1] <- medtarget-medref*model$coefficients[2];
							
							if (iqrdiff < 0.00001*iqrtarget) iqrdiff <- 0.00001*iqrtarget;

							f1 = iqrtarget/iqrdiff;
							
							dgf = nrow(cstrataref)-1;
							p <- 1.0-pf(dgf*(iqrtarget/iqrdiff)^2-dgf,1,dgf);
							
#							print(summary(model)$coef)
#							cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},
						LM = 
						{ 
							model <- lm(ftmp,data=cstrataref,na.action=na.exclude)
							f <- summary(model)$fstatistic
							f1 = f[1];
							p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},
						RLM = 
						{ 
							model <- MASS::rlm(ftmp,data=cstrataref,na.action=na.exclude,method = "MM")
							sw <- sum(model$w);
							dgf = sw-length(model$coef)+1;
							m1 <- sum(model$w*cstrataref[,colnamesList[i]],na.rm = TRUE)/sw
							rss1 <- sum(model$w*(cstrataref[,colnamesList[i]]^2),na.rm = TRUE)/sw-m1*m1
							m2 <- sum(model$w*model$residuals,na.rm = TRUE)/sw
							rss2 <- sum(model$w*(model$residuals^2),na.rm = TRUE)/sw-m2*m2
							f1 = rss1/rss2;
							p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},
						GLS =
						{
							model <- eval(parse(text=paste("try(nlme::gls(formula(",ftm1,"),cstrataref,na.action=na.exclude,correlation = nlme::corAR1(0.9,form = ~ 1 | ",correlationGroup,")))")))
							dgf = nrow(cstrataref)-length(model$coef)+1;
							rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)
							rss2 <- var(model$residuals,na.rm = TRUE);
							f1 = rss1/rss2;
							p1 <- 1.0-pf(dgf*rss1/rss2-dgf,1,dgf);
							reg <- summary(model);						
							p <- min(p1,reg$tTable[-1,4])
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p," ",p1,"\n");
						},
						{ 
							model <- lm(ftmp,data=cstrataref,na.action=na.exclude)
							f <- summary(model)$fstatistic
							f1 = f[1];
							p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						}
					)
					avg <- mean(cstrataref[,colnamesList[i]],na.rm = TRUE);
					models[[idx]] <- list(strata=sta,feature=colnamesList[i],model=model,mean=avg,pval=p);
					idx= idx+1;
			
					if (!is.na(p))
					{
						switch(type, 
							SPLINE = 
							{ 
								if (p<pvalue)
								{
									if (length(tbbaseModel) > 9)
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata[,baseModel])$y;
									}
									else
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
									}
								}
							},
							RQNT = 
							{ 
								if (p<pvalue)
								{
									cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
								}
							},
							LM = 
							{ 
								if (p<pvalue)
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata);
								}
								else
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-avg;
								}
							},
							RLM = 
							{ 
								if (p<pvalue)
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata);
								}
								else
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-avg;
								}
							},
							GLS =
							{
								if (p<pvalue)
								{
									avg <-model$coef[1];
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata)+avg;
								}
							}
						)
					}
				}
			}
			if (created == 1) 
			{
				AdjustedFrame <- rbind(AdjustedFrame,cstrata);
			}
			else
			{
				created = 1;
				AdjustedFrame = cstrata;
			}
		}
	}
	attr(AdjustedFrame,"models") <- models;
	# for (i in 1:size)		
	# { 
		# var1 <- var(data[,colnamesList[i]],na.rm = TRUE);
		# var2 <- var(AdjustedFrame[,colnamesList[i]],na.rm = TRUE);
 		# cat(" Variable: \t",colnamesList[i],"\t Var Ini: \t",var1,"\t Var End:\t",var2,"\t F:\t",var1/var2,"\n");
	# }
	if (created ==0 ) 
	{
		AdjustedFrame=NULL;
	}

	return (AdjustedFrame);
}
