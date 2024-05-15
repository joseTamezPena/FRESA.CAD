baggedModelS <-
function(modelFormulas,data,type=c("LM","LOGIT","COX"),Outcome=NULL,timeOutcome=NULL)
{
	type <- match.arg(type)

	formulaLoops <- 1;
	observations <- nrow(data);
	avgLogPvalues <- NULL;
	forder <- NULL;
#	cat(length(modelFormulas)," :Bagging\n",modelFormulas[1],"\n");
	model <- NULL;
	wts <- 0.0;
	if (inherits(modelFormulas,"list"))
	{ 
		listformulas <- modelFormulas;
		modelFormulas <- character();
		for (i in 1:length(listformulas))
		{
			modelFormulas <- append(modelFormulas,paste(listformulas[[i]],collapse = ' + '))
		}
	}
	coefList <- NULL;
	formulaList <- NULL;
	wtsList <- NULL;
	oF <- NULL;
	nnmodel <- NULL;
	if (length(modelFormulas) > 0)
	{
		if (modelFormulas[1] != "=-=End=-=")
		{
			iscompletef <- (gregexpr(pattern ='~',modelFormulas[1])[1] > 0)
			if (iscompletef && is.null(Outcome))
			{
				varsmod <- all.vars(formula(modelFormulas[1]));
				if (substr(modelFormulas[1],1,5)!="Surv(")
				{
					Outcome <- varsmod[1];
				}
				else
				{
					Outcome <- varsmod[2];
					timeOutcome <- varsmod[1];
				}
			}
			outcomes <- c(Outcome,timeOutcome);



			theoutcome <- data[,Outcome];
			varOutcome <- var(theoutcome);
			binoutcome <- (length(table(theoutcome))==2) && (min(theoutcome)==0);
			predtype="linear";
			if (binoutcome) predtype="prob";
			if ( (type=="LM") && (binoutcome==TRUE) )
			{
				data[,Outcome] = 2*theoutcome-1.0;
				predtype="linear";
			}
			
			if (type!="COX")
			{
				baseForm <- paste(Outcome,"~ 1");
			}
			else
			{
				baseForm <- paste("Surv(",timeOutcome,",",Outcome,")~ 1");
			}
			EquTrainSet <- data;
			minTrainSamples <- nrow(data);
			maxTrainSamples = minTrainSamples;
			casesample  <- NULL;
			controlsample <- NULL;
			noequalSets <- FALSE;
			nrowcases <- minTrainSamples;
			nrowcontrols <- minTrainSamples;
			nrep <- 1;
			frma <- baseForm;
			Jaccard.SM <- 0;
			coefEvolution <- NULL;
			avgsize <- 0;
			formulaNetwork <- NULL;
			WformulaNetwork <- NULL;
			coefList <- list();
			formulaList <- list();
			wtsList <- numeric();
			coef_in <- 1;

#				cat(length(modelFormulas)," :Entering orderFeatures\n");

			oF <- orderFeatures(modelFormulas,univariate=NULL,baseForm);
			
			VarFrequencyTable <- oF$VarFrequencyTable;
			if (!is.null(VarFrequencyTable))
			{
#				print(VarFrequencyTable);
				forder <- oF$forder;
				theterms <- oF$theterms;
				features <- oF$features;
				modelFormulas <- oF$modelFormulas;
				
				if (length(names(VarFrequencyTable)) > 0)
				{
					loops <- length(modelFormulas);
					vnames <- rownames(VarFrequencyTable);
					formulaNetwork <- matrix(0,nrow=length(VarFrequencyTable),ncol = length(VarFrequencyTable))
					dimnames(formulaNetwork) <- list(names(VarFrequencyTable),names(VarFrequencyTable))
					m <- formulaNetwork
					WformulaNetwork <- formulaNetwork

					Jaccard.SM <- 0;
					tota <- 0;
					for (n in 1:loops)
					{
						feat <- theterms[[n]];
						lft <- length(feat);
						m[,] <- 0;
						if (lft>0)
						{
							m[feat,feat]=1;
							if (n<loops)
							{
								for (i in (n+1):loops)
								{
									feat2 <- theterms[[i]];
									if (length(feat2) > 0)
									{
										Jaccard.SM = Jaccard.SM+sum(duplicated(c(feat2,feat)))/length(unique(c(feat2,feat)));
										tota = tota + 1;
									}
								}
							}
						}
						formulaNetwork <- formulaNetwork+m
					}
					if (tota>1) Jaccard.SM = Jaccard.SM/tota;
					fnrom <- loops;
					if (oF$numberofBreaks>0) fnrom <- oF$numberofBreaks;

					formulaNetwork <- round(formulaNetwork/fnrom,digits = 3);

				#	cat("Size :",nrow(data)," Features :",length(VarFrequencyTable))
					
					nsize <- nrow(data)
					
					lastTopVariable = length(VarFrequencyTable);
			#		if (lastTopVariable >= 2*(nsize-2)) lastTopVariable <- 2*(nsize-2); #The largest model size
					
					frma <- baseForm;
					enterlist <- vector();
					toRemove <- vector();
					bmodelsize <- 1;
					removed <- 0;
					avgsize <- 0;
					coefEvolution <- NULL;
					zthr <- 0;
					fistfreq <- VarFrequencyTable[1];
					for ( i in 1:length(vnames))
					{
						if ((vnames[i] != " ") && (vnames[i] != ""))
						{
							frma <- paste(frma,"+",vnames[i]);	
							bmodelsize = bmodelsize + 1;
						}
					}
#					cat("\nNum. Models:",loops," To Test:",length(vnames)," TopFreq:",fistfreq," Thrf:",thrsfreq," Removed:",removed,"\n")
					model <- modelFitting(frma,data,type=type,fitFRESA=TRUE);
					nnmodel <- model;
					if (bmodelsize>1)
					{
						thevars <- all.vars(formula(frma));
						
						if (inherits(model, "try-error"))
						{
			#				cat("Fitting Error\n");
							warning(frma," Warning Bagging Fitting error\n")
						}
						else
						{
							msize <- length(model$coefficients)
							basecoef <- abs(model$coefficients)+1e-6;
							names(basecoef) <- names(model$coefficients);
							
							if ((type=="COX") && !inherits(model,"fitFRESA"))
							{
								avgLogPvalues <- numeric(length(model$coefficients));
								names(avgLogPvalues) <- names(model$coefficients);
							}
							else
							{
								avgLogPvalues <- numeric(length(model$coefficients)-1);
								names(avgLogPvalues) <- names(model$coefficients)[-1];
							}
							addedZvalues <- avgLogPvalues;
							addedCoef <- avgLogPvalues;
							baggingAnalysis <- list();
							baggingAnalysis$uMS_values <- avgLogPvalues;
							baggingAnalysis$rMS_values <- avgLogPvalues;
							baggingAnalysis$NeRI_values <- avgLogPvalues;
							baggingAnalysis$pt_values <- avgLogPvalues;
							baggingAnalysis$pWilcox_values <- avgLogPvalues;
							baggingAnalysis$pF_values <- avgLogPvalues;
							baggingAnalysis$pBin_values <- avgLogPvalues;
							baggingAnalysis$mMSE_values <- avgLogPvalues;
							baggingAnalysis$dMSE_values <- avgLogPvalues;
							baggingAnalysis$uAcc_values <- avgLogPvalues;
							baggingAnalysis$rAcc_values <- avgLogPvalues;
							baggingAnalysis$uAUC_values <- avgLogPvalues;
							baggingAnalysis$rAUC_values <- avgLogPvalues;
							baggingAnalysis$dAUC_values <- avgLogPvalues;
							baggingAnalysis$idi_values <- avgLogPvalues;
							baggingAnalysis$nri_values <- avgLogPvalues;
							baggingAnalysis$zidi_values <- avgLogPvalues;
							baggingAnalysis$znri_values <- avgLogPvalues;
							baggingAnalysis$mAUC_values <- avgLogPvalues;
							baggingAnalysis$mACC_values <- avgLogPvalues;
							baggingAnalysis$coefstd <- avgLogPvalues;
							baggingAnalysis$coefficients <- avgLogPvalues;
							baggingAnalysis$wts <- avgLogPvalues;
						#	print(basecoef);
							avgsize <- msize-1;
							mado <- NA;
							rnames <- 0;
							if ((msize > 1)&&(loops>1))
							{
								model$type=type;
								onames <- names(model$coefficients);
								mmult <- 1+1*(type=="COX");
								model$estimations <- numeric(mmult*msize); 
								wts <- 0;
								model$coefficients <- numeric(msize);
								names(model$coefficients) <- onames;

								modelmeans <- model$coefficients;
								coefEvolution <- c(0,model$coefficients);
								names(coefEvolution) <- c("Weight",names(model$coefficients));
							#	cat("\n");
								tot_cycles <- 0;
								b_casesample <- casesample;
								b_controlsample <- controlsample;
#								print(model$coefficients);

								avgsize = 0;
								preddata <- data[,outcomes];
								predformula <- baseForm
								prenames <- outcomes
								onlypred <- NULL;
#								print(predtype)
								bestfit <- NULL
								bestformula <- baseForm
								allprednames <- NULL
								for (n in 1:loops)
								{
									if ((n %% 10) == 0) cat(".");
									feat <- theterms[[n]]
									if (length(feat)>0)
									{
										pname <- sprintf("V%d",n)
										out <- modelFitting(formula(modelFormulas[n]),data,type,fitFRESA=TRUE);
										preddata <- cbind(preddata,predict.fitFRESA(out,data,predtype));
										prenames <- c(prenames,pname);
										predformula <- paste(predformula," + ",pname,sep="")
										onlypred <- c(onlypred,pname);
#										predformula <- paste(bestformula," + ",pname,sep="")
#										preddata <- as.data.frame(preddata)
#										colnames(preddata) <- prenames;
#										wtsfit <- modelFitting(formula(predformula),preddata,type,fitFRESA=TRUE);
#										print(wtsfit$coefficients)
#										if (wtsfit$coefficients[pname]>0)
#										{
#											bestformula <- predformula
#											bestfit <- wtsfit
#											onlypred <- c(onlypred,pname);
#										}
#										allprednames <- c(allprednames,pname);
									}
								}
								preddata <- as.data.frame(preddata)
								colnames(preddata) <- prenames;
								wtsfit <- modelFitting(formula(predformula),preddata,type,fitFRESA=TRUE);
#								print(summary(wtsfit)$coef)
								allwts <- wtsfit$coefficients[onlypred];
#								print(predformula)
#								print(colnames(preddata))
#								wtsfit <- modelFitting(formula(bestformula),preddata,type,fitFRESA=TRUE);
#								print(onlypred);
#								allwts <- numeric(length(allprednames));
#								names(allwts) <- allprednames
#								allwts[onlypred] <- wtsfit$coefficients[onlypred];
#								print(allwts);
								fidx=1;

								for (n in 1:loops)
								{
									if ((n %% 10) == 0) cat(".");
									feat <- theterms[[n]]
									avgsize = avgsize+length(feat);
									if (length(feat)>0)
									{
											Rwts <- allwts[fidx]
#											print(Rwts);
											fidx <- fidx + 1;
											out <- modelFitting(formula(modelFormulas[n]),data,type,fitFRESA=TRUE);
											if (!inherits(out, "try-error")) 
											{
												osize <- length(out$coefficients)
												if (osize > 1)
												{
													coefList[[coef_in]] <- out$coefficients;
													formulaList[[coef_in]] <- modelFormulas[n];
													coef_in <- coef_in + 1;
													tot_cycles = tot_cycles + 1;
													onames <- names(out$coefficients);
													znames <- onames;
													if ((type!="COX") || inherits(out,"fitFRESA")) znames <- onames[-1];
													if (predtype=="linear")
													{
														gvar <- getVar.Res(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
													}
													else
													{
														gvar <- getVar.Bin(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
													}
													if (predtype=="linear")
													{
														baggingAnalysis$uMS_values[znames] <- baggingAnalysis$uMS_values[znames] + Rwts*gvar$unitestMSE;
														baggingAnalysis$rMS_values[znames] <- baggingAnalysis$rMS_values[znames] + Rwts*gvar$redtestMSE;
														baggingAnalysis$NeRI_values[znames] <- baggingAnalysis$NeRI_values[znames] + Rwts*gvar$NeRIs;
														baggingAnalysis$pF_values[znames] <- baggingAnalysis$pF_values[znames] + Rwts*log(gvar$FP.value);
														baggingAnalysis$pt_values[znames] <- baggingAnalysis$pt_values[znames] + Rwts*log(gvar$tP.value);
														baggingAnalysis$pBin_values[znames] <- baggingAnalysis$pBin_values[znames] + Rwts*log(gvar$BinP.value);
														baggingAnalysis$pWilcox_values[znames] <- baggingAnalysis$pWilcox_values[znames] + 
																								Rwts*log(gvar$WilcoxP.value);
														baggingAnalysis$mMSE_values[znames] <- baggingAnalysis$mMSE_values[znames] + Rwts*gvar$FullTestMSE;
														baggingAnalysis$dMSE_values[znames] <- (baggingAnalysis$dMSE_values[znames] + 
																								Rwts*(gvar$redtestMSE - gvar$FullTestMSE) );
													}
													else
													{
														baggingAnalysis$uAcc_values[znames] <- baggingAnalysis$uAcc_values[znames] + Rwts*gvar$uniTestAccuracy;
														baggingAnalysis$rAcc_values[znames] <- baggingAnalysis$rAcc_values[znames] + Rwts*gvar$redtestAccuracy;
														baggingAnalysis$uAUC_values[znames] <- baggingAnalysis$uAUC_values[znames] + Rwts*gvar$uniTestAUC;
														baggingAnalysis$rAUC_values[znames] <- baggingAnalysis$rAUC_values[znames] + Rwts*gvar$redtestAUC;
														baggingAnalysis$idi_values[znames] <- baggingAnalysis$idi_values[znames] + Rwts*gvar$IDIs;
														baggingAnalysis$nri_values[znames] <- baggingAnalysis$nri_values[znames] + Rwts*gvar$NRIs;
														baggingAnalysis$zidi_values[znames] <- baggingAnalysis$zidi_values[znames] + Rwts*gvar$z.IDIs;
														baggingAnalysis$znri_values[znames] <- baggingAnalysis$znri_values[znames] + Rwts*gvar$z.NRIs;
														baggingAnalysis$mAUC_values[znames] <- baggingAnalysis$mAUC_values[znames] + Rwts*gvar$fullTestAUC;
														baggingAnalysis$dAUC_values[znames] <- (baggingAnalysis$dAUC_values[znames] + 
																								Rwts*(gvar$fullTestAUC-gvar$redtestAUC) );
														baggingAnalysis$mACC_values[znames] <- baggingAnalysis$mACC_values[znames] + Rwts*gvar$fullTestAccuracy;
													}
													rnames <- append(rnames,tot_cycles)
													outmeans <- out$coefficients;
													wts = wts + Rwts;
													model$coefficients[onames] <- model$coefficients[onames] + Rwts*out$coefficients[onames];
													baggingAnalysis$coefficients[znames] <- baggingAnalysis$coefficients[znames] + Rwts*out$coefficients[znames];
													baggingAnalysis$coefstd[znames] <- baggingAnalysis$coefstd[znames] + Rwts*(out$coefficients[znames]^2);
													baggingAnalysis$wts[znames] <- baggingAnalysis$wts[znames] + rep(Rwts,length(znames));
													coefEvolution <- rbind(coefEvolution,c(Rwts,model$coefficients/wts));
#														print(model$coefficients/wts);
													addedZvalues[znames] <- addedZvalues[znames] + rep(Rwts,length(znames));
													addedCoef[znames] <- addedCoef[znames] + rep(1,length(znames));
													if (type=="COX")
													{
														fullmodelmeans <- abs(modelmeans[onames]);
														names(fullmodelmeans) <- onames 
														for (ei in 1:osize) 
														{
															outmeans[ei] <- out$estimations[osize+ei];
														}
														for (ei in onames) 
														{
															if (fullmodelmeans[ei]>0)
															{
																modelmeans[ei] <- 0.5*(modelmeans[ei] + outmeans[ei]); 
															}
															else
															{
																modelmeans[ei] <- outmeans[ei]; 
															}
														}
													}
#													print(model$coefficients)
												}
											}
											else
											{
												cat("+");
							#					print(out$coef);
											}
									}
								}
								avgsize = avgsize/loops;
								if( wts > 0)
								{
									conames <- names(model$coefficients);
									model$coefficients <- model$coefficients/wts;

									names(model$coefficients) <- conames;
#									print(model$coefficients);
									baggingAnalysis$coefficients <- baggingAnalysis$coefficients/baggingAnalysis$wts;
									gain <- model$coefficients[names(baggingAnalysis$coefficients)]/baggingAnalysis$coefficients;
									
									baggingAnalysis$coefstd <- gain*sqrt(abs(baggingAnalysis$coefstd/baggingAnalysis$wts-(baggingAnalysis$coefficients)^2));
									baggingAnalysis$coefficients <- gain*baggingAnalysis$coefficients;
									
									coefEvolution <- as.data.frame(coefEvolution);
									rownames(coefEvolution) <- rnames
									if (type == "COX")
									{
										model$estimations <- c(model$coefficients,modelmeans);
									}
									else
									{
										model$estimations <- model$coefficients;
									}
									avgLogPvalues <- avgLogPvalues/addedCoef;
									baggingAnalysis$formula.list <- modelFormulas;
									baggingAnalysis$uMS_values <- baggingAnalysis$uMS_values/addedZvalues;
									baggingAnalysis$rMS_values <- baggingAnalysis$rMS_values/addedZvalues;
									baggingAnalysis$NeRI_values <- baggingAnalysis$NeRI_values/addedZvalues;
									baggingAnalysis$pt_values <- exp(baggingAnalysis$pt_values/addedZvalues);
									baggingAnalysis$pWilcox_values <- exp(baggingAnalysis$pWilcox_values/addedZvalues);
									baggingAnalysis$pF_values <- exp(baggingAnalysis$pF_values/addedZvalues);
									baggingAnalysis$pBin_values <- exp(baggingAnalysis$pBin_values/addedZvalues);
									baggingAnalysis$uAcc_values <- baggingAnalysis$uAcc_values/addedZvalues;
									baggingAnalysis$rAcc_values <- baggingAnalysis$rAcc_values/addedZvalues;
									baggingAnalysis$uAUC_values <- baggingAnalysis$uAUC_values/addedZvalues;
									baggingAnalysis$rAUC_values <- baggingAnalysis$rAUC_values/addedZvalues;
									baggingAnalysis$idi_values <- baggingAnalysis$idi_values/addedZvalues;
									baggingAnalysis$nri_values <- baggingAnalysis$nri_values/addedZvalues;
									baggingAnalysis$zidi_values <- baggingAnalysis$zidi_values/addedZvalues;
									baggingAnalysis$znri_values <- baggingAnalysis$znri_values/addedZvalues;
									baggingAnalysis$mAUC_values <- baggingAnalysis$mAUC_values/addedZvalues;
									baggingAnalysis$dAUC_values <- baggingAnalysis$dAUC_values/addedZvalues;
									baggingAnalysis$mACC_values <- baggingAnalysis$mACC_values/addedZvalues;
									baggingAnalysis$mMSE_values <- baggingAnalysis$mMSE_values/addedZvalues;
									baggingAnalysis$dMSE_values <- baggingAnalysis$dMSE_values/addedZvalues;
									baggingAnalysis$wts <- baggingAnalysis$wts/addedCoef;
									baggingAnalysis$RelativeFrequency <- VarFrequencyTable/fnrom;
									baggingAnalysis$Jaccard.SM <- Jaccard.SM;
									baggingAnalysis$coeff_n_samples <- addedCoef;
									baggingAnalysis$observations <- observations;
									baggingAnalysis$avgLogPvalues <- avgLogPvalues;

#									model$linear.predictors <- predict(model);
									model$linear.predictors <- predict.fitFRESA(out,data,"linear");
									model$formula.List <- formulaList;
#									model$coefficients.List <- coefList;
									model$wts.List <- wtsList;
									model$fraction <- 1.0/nrow(data);

									model$baggingAnalysis <- baggingAnalysis;
									WformulaNetwork <- sqrt(WformulaNetwork);
									diag(WformulaNetwork) <- 0;
									WformulaNetwork <- WformulaNetwork/max(WformulaNetwork);
								}
								else
								{
									
									model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE)
						#			print(model$coefficients)
									model$coefficients[is.nan(model$coefficients)] <- 0.0;
									model$coefficients[is.na(model$coefficients)] <- 0.0;
									model$estimations[is.nan(model$estimations)] <- 0.0;
									model$estimations[is.na(model$estimations)] <- 0.0;
								}
							}
							else
							{
								model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE)
								model$coefficients[is.nan(model$coefficients)] <- 0.0;
								model$coefficients[is.na(model$coefficients)] <- 0.0;
								model$estimations[is.nan(model$estimations)] <- 0.0;
								model$estimations[is.na(model$estimations)] <- 0.0;
							}
						}
					}
				}
				else
				{
					model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE);
				}
			}
			else
			{
				model <- modelFitting(formula(frma),data,type=type,fitFRESA=TRUE);
			}
#			print(model$coefficients);
			environment(model$formula) <- globalenv();
			environment(model$terms) <- globalenv();		
		}
		else
		{
			warning("No Formulas\n");
			model <- NULL;
			frma <- NULL;
			VarFrequencyTable <- NULL;
			avgsize <- 0;
			formulaNetwork <- NULL;
			Jaccard.SM <- 0;
			coefEvolution <- NULL;
		}
	}
	else
	{
		warning("No Formulas\n");
		model <- NULL;
		frma <- NULL;
		VarFrequencyTable <- NULL;
		avgsize <- 0;
		formulaNetwork <- NULL;
		Jaccard.SM <- 0;
		coefEvolution <- NULL;
	}
	

  	result <- list(bagged.model=model,
				   formula=frma,
				   frequencyTable=VarFrequencyTable,
				   averageSize=avgsize,
				   formulaNetwork=formulaNetwork,
				   WformulaNetwork=WformulaNetwork,
				   Jaccard.SM = Jaccard.SM,
				   coefEvolution=coefEvolution,
				   featureLocation=forder,
				   predfit = wtsfit
				   );
  
	return (result);
}
