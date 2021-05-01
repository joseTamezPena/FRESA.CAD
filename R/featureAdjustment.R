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
	isContinous <- TRUE
	datamodel <- NULL;
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
			if (sum(str_count(baseModel,"\\+")) == 0)
			{
				datamodel <- cstrataref[,baseModel]
				isContinous <- length(table(datamodel)) > 5;
			}
			for (i in 1:size)
			{ 
				dtacolumn <- cstrataref[,colnamesList[i]]
				tb <- table(dtacolumn)
				if (length(tb) > 4)
				{
					avgref <- mean(dtacolumn,na.rm = TRUE);
					ress1 <- dtacolumn - avgref;
					ftm1 <- paste(colnamesList[i],paste(" ~ ",baseModel));
					ftmp <- formula(ftm1);
					mfref <- model.frame(ftmp,cstrataref);
					mfstrata <- model.frame(ftmp,cstrata);
					modellm <- try(lm(ftmp,data=mfref, model = FALSE,na.action=na.exclude))
					model <- modellm;
					plm <- 1.0;
					p <- 1.0;
					if (!inherits(modellm, "try-error"))
					{	
						ress2 <- modellm$residuals
						plm <- .Call("improvedResidualsCpp",ress1,ress2,"Wilcox",0)$p.value
						model <- modellm;
						p <- plm;
						if (isContinous)
						{
							plm <- cor.test(datamodel,dtacolumn,method="spearman")$p.value
						}
					}
					modelRLM <- modellm;
#					f <- summary(modellm)$fstatistic
#					pft <- pf(f[1],f[2],f[3],lower.tail=FALSE);
#					if (is.na(pft)) pft <- 1.0;
#					if (is.na(plm)) plm <- 1.0;
#					plm <- min(plm,0.5*pft)
					switch(type,
    					LOESS =
						{
							if ((plm < pvalue) && isContinous )
							{
								modelRLM <- try(MASS::rlm(ftmp,data=mfref,na.action=na.exclude, model = FALSE,method = "MM"), silent=TRUE)
								if (inherits(modelRLM, "try-error"))
								{	
									modelRLM <- modellm
									if (any(is.na(model$residuals)))
									{
										modelRLM <- modellm;
									}
								}
								model <- try(loess(ftmp,data=mfref,model=FALSE,...), silent=TRUE)
								if (!inherits(model, "try-error"))
								{
									pred <- predict(model,mfref);
									predlm <- predict(modelRLM,mfref);
									pred[is.na(pred)] <- predlm[is.na(pred)];
									ress <- dtacolumn - pred;
									p <- .Call("improvedResidualsCpp",ress1,ress,"Wilcox",0)$p.value
#									sdd <- 3.0*median(abs(ress)) + 2.0*mean(abs(ress))
#									if (sdd == 0)
#									{
#										sdd <- 1.0
#									}									
#									wts <- exp(-(ress/sdd)^2)
#									dgf = length(ress)*min(model$pars$span,1.0)-model$pars$degree;
#									rss1 <- sum(wts*ress1^2)
#									rss2 <- sum(wts*ress^2)
#									pft <- 0;
#									if (rss2 > 0)
#									{
#										pft <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE)
#									}
#									if (is.na(p)) p <- 1.0;
#									if (is.na(pft)) pft <- 1.0;
#									p <- min(p,0.5*pft);
								}								
								else
								{
									model <- modelRLM;
								}								
							}
						},						
    					MARS =
						{
							if ((plm < pvalue) && isContinous )
							{
								model <- try(mda::mars(datamodel, 
												  y = dtacolumn,
												  ...
												  ), silent=TRUE
											  )
								if (!inherits(model, "try-error"))
								{
									pred <- as.numeric(predict(model,datamodel));
									ress <- dtacolumn - pred;
									p <- .Call("improvedResidualsCpp",ress1,ress,"Wilcox",0)$p.value

#									sdd <- 3.0*median(abs(ress)) + 2.0*mean(abs(ress))
#									if (sdd == 0)
#									{
#										sdd <- 1.0
#									}									
#									wts <- exp(-(ress/sdd)^2)
#									dgf = length(ress) - length(model$coefficients);
#									rss1 <- sum(wts*ress1^2)
#									rss2 <- sum(wts*ress^2)
#									pft <- 0;
#									if (rss2 > 0)
#									{
#										pft <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE)
#									}
#									if (is.na(p)) p <- 1.0;
#									if (is.na(pft)) pft <- 1.0;
#									p <- min(p,0.5*pft);
								}								
								else
								{
									model <- modellm;
								}								
							}
						},						
    					SPLINE =
						{
							if ((plm < pvalue) && isContinous )
							{
								model <- try(smooth.spline(datamodel, 
												  y = dtacolumn,
												  keep.data = FALSE,
												  ...
												  ), silent=TRUE
											  )
								if (!inherits(model, "try-error"))
								{
									pred <- predict(model,datamodel);
									ress <- dtacolumn - pred$y;
									p <- .Call("improvedResidualsCpp",ress1,ress,"Wilcox",0)$p.value

#									sdd <- 3.0*median(abs(ress)) + 2.0*mean(abs(ress))
#									if (sdd == 0)
#									{
#										sdd <- 1.0
#									}									
#									wts <- exp(-(ress/sdd)^2)
#									dgf = length(ress) - model$fit$nk;
#									rss1 <- sum(wts*ress1^2)
#									rss2 <- sum(wts*ress^2)
#									pft <- 0;
#									if (rss2 > 0)
#									{
#										pft <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE)
#									}
#									if (is.na(p)) p <- 1.0;
#									if (is.na(pft)) pft <- 1.0;
#									p <- min(p,0.5*pft);
								}
								else
								{
									model <- modellm;
								}								
							}
						},						
						RLM = 
						{ 
							if ((plm < pvalue) && isContinous)
							{
								model <- try(MASS::rlm(ftmp,data=mfref,na.action=na.exclude, model = FALSE,x.ret = FALSE,method = "MM"), silent=TRUE)
								if (!inherits(model, "try-error"))
								{								
									ress <- predict(model,mfref)-dtacolumn;
									if ( any(is.na(model$w)) | any(is.na(ress)) )
									{
										model <- modellm;
									}
									else
									{
										p <- .Call("improvedResidualsCpp",ress1,ress,"Wilcox",0)$p.value

#										model$w[is.na(model$w)] <- 1.0;
#										model$w[model$w == 0] <- 1.0e-5;
#										sw <- sum(model$w);
#										dgf = length(ress)-length(model$coef)+1;
#										m1 <- sum(model$w*dtacolumn,na.rm = TRUE)/sw
#										rss1 <- sum(model$w*(dtacolumn^2),na.rm = TRUE)/sw-m1*m1
#										rss2 <- sum(model$w*(ress^2),na.rm = TRUE)/sw
#										pft <- 0;
#										if (rss2 > 0)
#										{
#											pft <- pf(dgf*rss1/rss2-dgf,1,dgf,lower.tail=FALSE)
#										}
#										if (is.na(p)) p <- 1.0;
#										if (is.na(pft)) pft <- 1.0;
#										p <- min(p,0.5*pft);
									}
								}
								else
								{
									model <- modellm;
								}								
							}
						},
						GLS =
						{
							model <- eval(parse(text=paste("try(nlme::gls(formula(",ftm1,"),cstrataref,na.action=na.exclude,correlation = nlme::corAR1(0.9,form = ~ 1 | ",correlationGroup,")))")))
							ress <- model$residuals
							p <- .Call("improvedResidualsCpp",ress1,ress,"Wilcox",0)$p.value

							dgf = length(ress)-length(model$coef)+1;
							rss1 <- sum(ress1^2)
							rss2 <- sum(ress^2,na.rm = TRUE);
							p1 <- 1.0-pf(dgf*rss1/rss2-dgf,1,dgf);
							reg <- summary(model);						
							p <- min(c(p,p1,reg$tTable[-1,4]))
						}
					)
					if (is.na(p)) p <- 1.0;
					if (is.nan(p)) p <- 1.0;
					environment(model$formula) <- NULL;
					environment(model$terms) <- NULL;
					models[[idx]] <- list(strata=sta,feature=colnamesList[i],model=model,avgref=avgref,pval=p);
					idx= idx+1;
			
					if (!is.na(p))
					{
						orgdata <- cstrata[,colnamesList[i]];
						switch(type, 
							LOESS = 
							{ 
								if (p < pvalue)
								{
									pred <- predict(model,mfstrata);
									predlm <- predict(modelRLM,mfstrata);
									pred[is.na(pred)] <- predlm[is.na(pred)];
									cstrata[,colnamesList[i]] <-  avgref + cstrata[,colnamesList[i]] - pred;
								}
							},
							MARS = 
							{ 
								if (p < pvalue)
								{
									if (class(model) == "mars")
									{
										cstrata[,colnamesList[i]] <- avgref + cstrata[,colnamesList[i]] - as.numeric(predict(model,cstrata[,baseModel]));
									}
									else
									{
										cstrata[,colnamesList[i]] <- avgref + cstrata[,colnamesList[i]] - predict(model,mfstrata);
									}
								}
							},
							SPLINE = 
							{ 
								if (p < pvalue)
								{
									if (class(model) == "smooth.spline")
									{
										cstrata[,colnamesList[i]] <-  avgref + cstrata[,colnamesList[i]] - predict(model,cstrata[,baseModel])$y;
									}
									else
									{
										cstrata[,colnamesList[i]] <- avgref + cstrata[,colnamesList[i]] - predict(model,mfstrata);
									}
								}
							},
							NZLM = 
							{ 
								if (p < pvalue)
								{
									cstrata[,colnamesList[i]] <- model$coef[1] + cstrata[,colnamesList[i]] - predict(model,mfstrata);
								}
							},
							LM = 
							{ 
								if (p < pvalue)
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]] - predict(model,mfstrata);
								}
								else
								{
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]] - model$coef[1];
								}
							},
							RLM = 
							{ 
								if (p < pvalue)
								{
									pred <- predict(model,mfstrata) - model$coef[1];
									if (any(is.na(pred)))
									{
										pred <- predict(modellm,mfstrata) - modellm$coef[1];;
									}
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]] - pred;
								}
							},
							GLS =
							{
								if (p < pvalue)
								{
									avg <- model$coef[1];
									cstrata[,colnamesList[i]] <- cstrata[,colnamesList[i]] - predict(model,mfstrata) + avg;
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
