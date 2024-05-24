baggedModel <-
function(modelFormulas,data,type=c("LM","LOGIT","COX"),Outcome=NULL,timeOutcome=NULL,frequencyThreshold=0.025,univariate=NULL,useFreq=TRUE,n_bootstrap=1,equifreqCorrection=0)
{
	type <- match.arg(type)

	formulaLoops <- 1;
	if (is.numeric(useFreq))
	{
		formulaLoops <- useFreq;
		useFreq <- TRUE;
	}
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

#		cat(length(modelFormulas)," :in 1 Bagging\n");

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
			if (type != "LM")
			{
				casesample = subset(data,get(Outcome)  == 1);
				controlsample = subset(data,get(Outcome) == 0);
				nrowcases <- nrow(casesample);
				nrowcontrols <- nrow(controlsample);
				
				minTrainSamples <- min(c(nrowcases,nrowcontrols));
				maxTrainSamples <- max(c(nrowcases,nrowcontrols));
#				noequalSets <- (minTrainSamples < (0.75*maxTrainSamples));
				nrep <- 1;
				if (noequalSets) nrep <- 1 + (as.integer(3*maxTrainSamples/minTrainSamples + 0.5));
#				cat(nrowcases,":",nrowcontrols,":",noequalSets,"Boot :",nrep,"\n")
		#		cat(minTrainSamples,":",maxTrainSamples,":",noequalSets,"\n")
			}
#			cat(length(modelFormulas)," :in 2 Bagging\n");
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

			oF <- orderFeatures(modelFormulas,univariate,baseForm,useFreq,n_bootstrap);
			
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
				#						print(feat2)
				#						print(feat)
				#						cat("Dup:",sum(duplicated(c(feat2,feat)))," U:",length(unique(c(feat2,feat)))," JI:",Jaccard.SM,"\n")
									}
								}
							}
						}
						formulaNetwork <- formulaNetwork+m
					}
					if (tota>1) Jaccard.SM = Jaccard.SM/tota;
					fnrom <- loops;
					if (oF$numberofBreaks>0) fnrom <- oF$numberofBreaks;
					if (equifreqCorrection>0) fnrom <- fnrom*equifreqCorrection;

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
					if (!is.null(univariate))
					{
						nfeat <- nrow(univariate);
						zthr <- abs(qnorm(0.1/(formulaLoops*nfeat)));
					}
					thrsfreq <- as.integer(frequencyThreshold*VarFrequencyTable[1]+0.5);
					if (formulaLoops==1) 
					{
						fthrsfreq <- as.integer(frequencyThreshold*loops+0.5);
					}
					else
					{
						fthrsfreq <- as.integer(frequencyThreshold*formulaLoops+0.5);
					}
					fistfreq <- VarFrequencyTable[1];
					for ( i in 1:length(vnames))
					{
						if ((vnames[i] != " ") && (vnames[i] != ""))
						{
							enterlist <- append(enterlist,vnames[i]);
							passZ <- TRUE;
							if (!is.null(univariate) && (VarFrequencyTable[i]<=(fthrsfreq+1)) && (useFreq>0) && (frequencyThreshold>0))
							{
								passZ <- (univariate[vnames[i],"ZUni"]>zthr); # only features that have strong unit pvalue for very low frequencies
		#						if (passZ==FALSE)
		#						{
		#							cat(fthrsfreq," Removing:",vnames[i],":",VarFrequencyTable[i]," Z=",univariate[vnames[i],"ZUni"],"Zthr=",zthr,"\n");
		#						}
							}
							if (passZ&&(i<=lastTopVariable)&&(VarFrequencyTable[i] > thrsfreq))  # Only features with a given frequency
							{
								if ((fistfreq == loops)&&(VarFrequencyTable[i] > (loops/3)))
								{
									fistfreq <- VarFrequencyTable[i];
									thrsfreq <- as.integer(frequencyThreshold*fistfreq+0.5);
								}
								frma <- paste(frma,"+",vnames[i]);	
								bmodelsize = bmodelsize + 1;
							}
							else
							{
								toRemove <- append(toRemove,paste(" ",vnames[i]," ",sep=""));
								removed = removed+1;
							}
						}
					}
#					cat("\nNum. Models:",loops," To Test:",length(vnames)," TopFreq:",fistfreq," Thrf:",thrsfreq," Removed:",removed,"\n")
					model <- modelFitting(frma,data,type=type,fitFRESA=TRUE);
					nnmodel <- model;
					if (bmodelsize>1)
					{
						thevars <- all.vars(formula(frma));
#						data <- data[,thevars];
						if (noequalSets)
						{
							casesample <- casesample[,thevars]
							controlsample <- controlsample[,thevars]
							trainCaseses <- casesample;
							trainControls <- controlsample;
							EquTrainSet <- rbind(trainCaseses,trainControls)
#							cat("Cases: ",nrow(trainCaseses),"Controls: ",nrow(trainControls),"\n");
						}
						else
						{
							EquTrainSet <- data;
						}
						
						
					#	print(toRemove);
					#	print(model$coefficients);
					#	cat(frma,"\n");
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

								for (m in 1:n_bootstrap)
								{
									if (n_bootstrap>1)
									{
										if (type!="LM")
										{
											b_casesample <- casesample[sample(1:nrowcases,nrowcases,replace=TRUE),]
											b_controlsample <- controlsample[sample(1:nrowcontrols,nrowcontrols,replace=TRUE),]						
											EquTrainSet <- rbind(b_casesample,b_controlsample)
										}
										else
										{
											EquTrainSet <- data[sample(1:nrow(data),nrow(data),replace=TRUE),];
										}
										theoutcome <- EquTrainSet[,Outcome];
										varOutcome <- var(theoutcome);
									}
									avgsize = 0;
									for (n in 1:loops)
									{
										if ((n %% 10) == 0) cat(".");
										feat <- theterms[[n]]
										avgsize = avgsize+length(feat);
										if (m==1)
										{
								#			cat(modelFormulas[n],"\n");

											if (length(toRemove)>0)
											{
												modelFormulas[n] <- paste(modelFormulas[n],"  ",sep="");
												for (rml in 1:length(toRemove))
												{
													modelFormulas[n] <- sub(toRemove[rml]," 1 ",modelFormulas[n],fixed=TRUE);
												}
							#					cat("After Rem:",modelFormulas[n],"\n");
												feat <- attr(terms(formula(modelFormulas[n])),"term.labels");
											}
										}
										if (length(feat)>0)
										{
											for (replicates in 1:nrep)
											{
												if (noequalSets)
												{
													if (minTrainSamples < nrowcases) trainCaseses <- b_casesample[sample(1:nrowcases,minTrainSamples,replace=FALSE),]
													if (minTrainSamples < nrowcontrols) trainControls <- b_controlsample[sample(1:nrowcontrols,minTrainSamples,replace=FALSE),]
													EquTrainSet <- rbind(trainCaseses,trainControls)
													theoutcome <- EquTrainSet[,Outcome];
													varOutcome <- var(theoutcome);
												}
												out <- modelFitting(formula(modelFormulas[n]),EquTrainSet,type,fitFRESA=TRUE);
												coef_Panalysis <- NULL;
												if (!inherits(out, "try-error")) 
												{
													osize <- length(out$coefficients)
													if (osize > 1)
													{
														coefList[[coef_in]] <- out$coefficients;
														formulaList[[coef_in]] <- modelFormulas[n];
														coef_in <- coef_in +1;
														tot_cycles = tot_cycles+1;
#														curprediction <- predict.fitFRESA(out,EquTrainSet,predtype)
#														mpred <- cbind(mpred,predict.fitFRESA(out,data,predtype));
#														residual <- as.vector(abs(curprediction-theoutcome));
														onames <- names(out$coefficients);
														znames <- onames;
														if ((type!="COX") || inherits(out,"fitFRESA")) znames <- onames[-1];
														if (predtype=="linear")
														{
															gvar <- getVar.Res(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
															coef_Panalysis <- -log10(gvar$FP.value);
															fitScore <- (varOutcome-gvar$FullTestMSE)/varOutcome;
															if (fitScore < 0) fitScore <- 0;
															coefscore <- (gvar$redtestMSE-gvar$FullTestMSE)/varOutcome;
														}
														else
														{
															gvar <- getVar.Bin(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
			#												cat("Equ: ",mean(EquTrainSet[,Outcome])," Data: ",mean(data[,Outcome]),"\n");
															coef_Panalysis <- -log10(1.0-pnorm(gvar$z.IDIs));
															fitScore <- 2.0*(gvar$fullTestAUC - 0.5);
															if (fitScore < 0) fitScore <- 0;
															fitScore <- fitScore*fitScore;
															coefscore <- 2.0*(gvar$fullTestAUC-gvar$redtestAUC);
															coefscore[coefscore < 0] <- 0;
															coefscore <- coefscore^2;
														}
#														print(coef_Panalysis);
														coefscore[coefscore <= 0] <- 1.0e-4;
														infnum <- is.infinite(coef_Panalysis)
														if (sum(infnum)>0)
														{
#															print(coef_Panalysis);
															coef_Panalysis[coef_Panalysis == Inf] <- 100.0;
														}
														coef_Panalysis[coef_Panalysis > 100.0] <- 100.0;
														WformulaNetwork[znames,znames] <- WformulaNetwork[znames,znames] + coefscore%*%t(coefscore)
														Rwts <- sum(coef_Panalysis) - 2*length(coef_Panalysis);
														if (Rwts <= 0) Rwts <- 1.0e-4;
														Rwts <- Rwts*(fitScore + mean(coefscore));
														if (Rwts <= 0) Rwts <- 1.0e-4;
														wtsList <- c(wtsList,Rwts);
														if (predtype=="linear")
														{
#															gvar <- getVar.Res(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
#															coef_Panalysis <- -log(gvar$FP.value); 
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
#															gvar <- getVar.Bin(out,data=EquTrainSet,Outcome=Outcome,type=type,testData=data)
			#												cat("Equ: ",mean(EquTrainSet[,Outcome])," Data: ",mean(data[,Outcome]),"\n");
#															coef_Panalysis <- -log(1.0-pnorm(gvar$z.IDIs));
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
														avgLogPvalues[znames] <- avgLogPvalues[znames] + coef_Panalysis;
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
								#						print(model$coefficients)
													}
												}
												else
												{
													cat("+");
								#					print(out$coef);
												}
											}
										}
									}
									avgsize = avgsize/loops;
				#					cat("*");
								}
								if( wts > 0)
								{
				#					print(baggingAnalysis$coefficients^2);
				#					print(baggingAnalysis$coefstd);
				#					print(addedZvalues);

#									print(wts);
#									print(model$coefficients);
									conames <- names(model$coefficients);
									model$coefficients <- model$coefficients/wts;
									model$formula.List <- formulaList;
									model$coefficients.List <- coefList;
									model$wts.List <- wtsList;
									model$fraction <- 1.0/nrow(data);

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
									baggingAnalysis$n_bootstrap <- n_bootstrap;
									baggingAnalysis$coeff_n_samples <- addedCoef;
									baggingAnalysis$observations <- observations;
									baggingAnalysis$avgLogPvalues <- avgLogPvalues;

									model$baggingAnalysis <- baggingAnalysis;
									model$linear.predictors <- predict(model);
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
			model$model <- NULL;
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
	
	if (!is.null(oF))
	{
		modelFormulas <- oF$modelFormulas;
		if (length(names(VarFrequencyTable)) > 0)
		{
			loops <- length(modelFormulas);
			if (loops>1)
			{
				mpred <- NULL;
				out <- list();
				for (n in 1:loops)
				{
					out[[n]] <- modelFitting(formula(modelFormulas[n]),data,type,fitFRESA=TRUE);
					pred <- predict.fitFRESA(out[[n]],data,predtype);
					mpred <- cbind(mpred,as.numeric(pred));
					out[[n]] <- out[[n]]$coef
					out[[n]] <- out[[n]][names(out[[n]]) %in% names(nnmodel$coefficients)]
				}
				mpred <- as.data.frame(mpred)
				mpred$Outcome <- data[,Outcome] 
				nnmodel$coefficients <- nnmodel$coefficients*0;
				
				if (!is.null(timeOutcome))
				{
					mpred$timeOutcome <- data[,timeOutcome] 
				}
				twts <- 0;
				if (type!="COX")
				{
#					print(modelFormulas)
					if (loops<10)
					{
						wnnmodel <- modelFitting(formula("Outcome~."),mpred,type,fitFRESA=FALSE);
						for (n in 1:loops)
						{
							nnmodel$coefficients[names(out[[n]])] <- nnmodel$coefficients[names(out[[n]])] + wnnmodel$coef[n+1]*out[[n]]
						}
						twts <- sum(wnnmodel$coef[2:(loops+1)])
						wnnmodel <- NULL
					}
					else
					{
						for (lrep in 1:(floor(3*loops/10+1.0)))
						{
							whichsample <- c(1:3,sample(c(4:loops),6,replace=FALSE));
							snames <- c(colnames(mpred)[whichsample],"Outcome")
							wnnmodel <- modelFitting(formula("Outcome~."),mpred[,snames],type,fitFRESA=FALSE);
#							print(wnnmodel$coef)
							for (n in 1:length(whichsample))
							{
								wmodel <- out[[whichsample[n]]];
								nnmodel$coefficients[names(wmodel)] <- nnmodel$coefficients[names(wmodel)] + wnnmodel$coef[n+1]*wmodel
								twts <- twts + wnnmodel$coef[n+1];
							}
							wnnmodel <- NULL
						}
					}
					nnmodel$coefficients <- nnmodel$coefficients/twts;
					nnmodel$estimations <- nnmodel$coefficients
				}
				else
				{
					if (loops<10)
					{
						wnnmodel <- modelFitting(formula("Surv(timeOutcome,Outcome)~."),mpred,type,fitFRESA=FALSE);
						for (n in 1:loops)
						{
							nnmodel$coefficients[names(out[[n]])] <- nnmodel$coefficients[names(out[[n]])] + wnnmodel$coef[n]*out[[n]]
						}
						twts <- sum(wnnmodel$coef)
					}
					else
					{
						for (lrep in 1:(floor(3*loops/10+1.0)))
						{
							whichsample <- c(1:3,sample(c(4:loops),6,replace=FALSE));
							snames <- c(colnames(mpred)[whichsample],"Outcome","timeOutcome")
							wnnmodel <- modelFitting(formula("Surv(timeOutcome,Outcome)~."),mpred[,snames],type,fitFRESA=FALSE);
#							print(wnnmodel$coef)
							for (n in 1:length(whichsample))
							{
								wmodel <- out[[whichsample[n]]];
								nnmodel$coefficients[names(wmodel)] <- nnmodel$coefficients[names(wmodel)] + wnnmodel$coef[n]*wmodel
								twts <- twts + wnnmodel$coef[n];
							}
							wnnmodel <- NULL
						}
					}
					nnmodel$coefficients <- nnmodel$coefficients/twts;
					nnmodel$estimations <- model$estimations
					nnmodel$estimations[1:length(nnmodel$coefficients)] <- nnmodel$coefficients
					wnnmodel <- NULL
				}
			}
			environment(nnmodel$formula) <- globalenv();
			environment(nnmodel$terms) <- globalenv();		
		}
	}

  	result <- list(bagged.model=model,
				   formula=frma,
				   frequencyTable=VarFrequencyTable,
				   averageSize=avgsize,
				   formulaNetwork=formulaNetwork,
				   WformulaNetwork=WformulaNetwork,
				   Jaccard.SM = Jaccard.SM,
				   coefEvolution=coefEvolution,
				   avgLogPvalues=avgLogPvalues,
				   featureLocation=forder,
				   nnmodel=nnmodel
				   );
  
	return (result);
}
