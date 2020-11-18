featureAdjustment <-
function(variableList,baseModel,strata=NA,data,referenceframe,type=c("LM","GLS","RLM","NZLM","SPLINE","MARS","LOESS"),pvalue=0.05,correlationGroup = "ID",...) 
{

if (!requireNamespace("nlme", quietly = TRUE)) {
   install.packages("nlme", dependencies = TRUE)
} 
if (!requireNamespace("MASS", quietly = TRUE)) {
	install.packages("MASS", dependencies = TRUE)
}
if (!requireNamespace("mda", quietly = TRUE)) {
	install.packages("mda", dependencies = TRUE)
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
	AdjustedFrame <- NULL;
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
			tbbaseModel <- table(cstrataref[,baseModel])
			for (i in 1:size)
			{ 
				tb <- table(cstrataref[,colnamesList[i]])
				if (length(tb) > 4)
				{
					ftm1 <- paste(colnamesList[i],paste(" ~ ",baseModel));
					ftmp <- formula(ftm1);

					switch(type,
					
    					LOESS =
						{
							rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)

							dgf = 0.5*nrow(cstrataref) - 2;
							model <- try(loess(ftmp,data=cstrataref,model=FALSE,...)
										  )
							modellm <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)

							if (inherits(model, "try-error"))
							{
#								cat("Error\n");
								model <- modellm;						
								rss2 <- var(model$residuals,na.rm = TRUE);
							}
							else
							{
#								cat(class(model));
								dgf = nrow(cstrataref)*model$pars$span-model$pars$degree;
								pred <- predict(model,cstrataref);
								predlm <- predict(modellm,cstrataref);
								pred[is.na(pred)] <- predlm[is.na(pred)];
								ress <- cstrataref[,colnamesList[i]] - pred;
								rss2 <- var(ress,na.rm = TRUE);
							}								
							f1 = rss1/rss2;
							p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
#							cat(colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},						
    					MARS =
						{
#							cat(baseModel,": Variable:\t ",colnamesList[i]);
							rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)
							dgf = nrow(cstrataref) - 1;
#								cat(baseModel,":",colnamesList[i],"\n")
							model <- try(mda::mars(cstrataref[,baseModel], 
											  y = cstrataref[,colnamesList[i]],
											  ...
											  )
										  )
#								cat(class(model),"\n")
							if (inherits(model, "try-error"))
							{
#								cat("Error\n");
								model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude);						
								rss2 <- var(model$residuals,na.rm = TRUE);
							}
							else
							{
#								cat(class(model));
								dgf = nrow(cstrataref) - length(model$coefficients);
								pred <- as.numeric(predict(model,cstrataref[,baseModel]));
								ress <- cstrataref[,colnamesList[i]] - pred;
								rss2 <- var(ress,na.rm = TRUE);
							}								
							f1 = rss1/rss2;
							p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
#							cat(colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},						
    					SPLINE =
						{
#							cat(baseModel,": Variable:\t ",colnamesList[i]);
							rss1 <- var(cstrataref[,colnamesList[i]],na.rm = TRUE)
							dgf = nrow(cstrataref) - 1;
							if (length(tbbaseModel) > 5)
							{	
#								cat(baseModel,":",colnamesList[i],"\n")
								model <- try(smooth.spline(cstrataref[,baseModel], 
												  y = cstrataref[,colnamesList[i]],
												  keep.data = FALSE,
												  ...
												  )
											  )
#								cat(class(model),"\n")
								if (inherits(model, "try-error"))
								{
#									cat("Error\n");
									model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude);						
									rss2 <- var(model$residuals,na.rm = TRUE);
								}
								else
								{
									dgf = nrow(cstrataref) - model$df;
									pred <- predict(model,cstrataref[,baseModel]);
									ress <- cstrataref[,colnamesList[i]] - pred$y;
									rss2 <- var(ress,na.rm = TRUE);
								}								
							}
							else
							{
								model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude);						
								dgf = nrow(cstrataref) - 1;
								rss2 <- var(model$residuals,na.rm = TRUE);
							}
							f1 = rss1/rss2;
							p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
#							cat(colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},						
						NZLM = 
						{ 
							model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)
							f <- summary(model)$fstatistic
							f1 = f[1];
							p <- pf(f[1],f[2],f[3],lower.tail=FALSE)										
						},
						LM = 
						{ 
							model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)
							f <- summary(model)$fstatistic
							f1 = f[1];
							p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						},
						RLM = 
						{ 
							if (length(tbbaseModel) > 5)
							{
								model <- try(MASS::rlm(ftmp,data=cstrataref,na.action=na.exclude, model = FALSE,method = "MM"))
								if (!inherits(model, "try-error"))
								{								
									sw <- sum(model$w);
									dgf = sw-length(model$coef)+1;
									m1 <- sum(model$w*cstrataref[,colnamesList[i]],na.rm = TRUE)/sw
									rss1 <- sum(model$w*(cstrataref[,colnamesList[i]]^2),na.rm = TRUE)/sw-m1*m1
									m2 <- sum(model$w*model$residuals,na.rm = TRUE)/sw
									rss2 <- sum(model$w*(model$residuals^2),na.rm = TRUE)/sw-m2*m2
									f1 = rss1/rss2;
									p <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE);
								}
								else
								{
									model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)
									f <- summary(model)$fstatistic
									f1 = f[1];
									p <- pf(f[1],f[2],f[3],lower.tail=FALSE)										
								}
							}
							else
							{
								model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)
								f <- summary(model)$fstatistic
								f1 = f[1];
								p <- pf(f[1],f[2],f[3],lower.tail=FALSE)										
							}	
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
							model <- lm(ftmp,data=cstrataref, model = FALSE,na.action=na.exclude)
							f <- summary(model)$fstatistic
							f1 = f[1];
							p <- pf(f[1],f[2],f[3],lower.tail=FALSE)			
	#						cat(" Variable:\t ",colnamesList[i],"\t F Stats:\t ",f1,"\t P-value:\t",p,"\n");
						}
					)
					environment(model$formula) <- NULL;
					environment(model$terms) <- NULL;
					avg <- mean(cstrataref[,colnamesList[i]],na.rm = TRUE);
					models[[idx]] <- list(strata=sta,feature=colnamesList[i],model=model,mean=avg,pval=p);
					idx= idx+1;
			
					if (!is.na(p))
					{
						switch(type, 
							LOESS = 
							{ 
								if (p < pvalue)
								{
									pred <- predict(model,cstrata);
									predlm <- predict(modellm,cstrata);
									pred[is.na(pred)] <- predlm[is.na(pred)];
									cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-pred;
								}
							},
							MARS = 
							{ 
								if (p < pvalue)
								{
									if (class(model) == "mars")
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-as.numeric(predict(model,cstrata[,baseModel]));
									}
									else
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
									}
								}
							},
							SPLINE = 
							{ 
								if (p < pvalue)
								{
									if (class(model) == "smooth.spline")
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata[,baseModel])$y;
									}
									else
									{
										cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
									}
								}
							},
							NZLM = 
							{ 
								if (p < pvalue)
								{
									cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
								}
							},
							LM = 
							{ 
								if (p < pvalue)
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
								if (p < pvalue)
								{
									cstrata[,colnamesList[i]] <- avg + cstrata[,colnamesList[i]]-predict(model,cstrata);
								}
							},
							GLS =
							{
								if (p < pvalue)
								{
									avg <-model$coef[1];
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]]-predict(model,cstrata) + avg;
								}
							}
						)
					}
				}
			}
			AdjustedFrame <- rbind(AdjustedFrame,cstrata);
			created = 1;
		}
	}
	# for (i in 1:size)		
	# { 
		# var1 <- var(data[,colnamesList[i]],na.rm = TRUE);
		# var2 <- var(AdjustedFrame[,colnamesList[i]],na.rm = TRUE);
 		# cat(" Variable: \t",colnamesList[i],"\t Var Ini: \t",var1,"\t Var End:\t",var2,"\t F:\t",var1/var2,"\n");
	# }
	if (created > 0) 
	{
		attr(AdjustedFrame,"models") <- models;
	}

	return (AdjustedFrame);
}
