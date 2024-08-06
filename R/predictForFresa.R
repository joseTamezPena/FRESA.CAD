#' @method predict fitFRESA
predict.fitFRESA <-
function (object,...) 
{

	testData <- NULL
	impute=FALSE;
	parameters <- list(...);
	if (length(parameters)==3)
	{
		if (!is.null(parameters$impute))
		{
			impute=parameters$impute;
		}
	}
	if (length(parameters)==2)
	{
		testData <- parameters[[1]];
		predictType <- parameters[[2]];
	}
	else
	{
		if (length(parameters)==1)
		{
			testData <- parameters[[1]];
			predictType <- "linear";
			if (!is.null(object$type))
			{
				if (object$type== "LOGIT")
				{
					predictType <- "prob";
				}
			}			
#			if (!is.null(object$coefficients.List))
#			{
#				predictType <- "bagg";
#			}
		}
		else
		{
			if (is.null(parameters$testData))
			{
				testData <- object$model;
			}
			else
			{
				testData <- parameters$testData
			}
			if (is.null(parameters$predictType))
			{
				predictType <- "linear";
				if (!is.null(object$type))
				{
					if (object$type== "LOGIT")
					{
						predictType <- "prob";
					}
				}			
#				if (!is.null(object$coefficients.List))
#				{
#					predictType <- "bagg";
#				}
			}
			else
			{
				predictType <- parameters$predictType
			}
		}
	}
	
#	print(as.character(frm))
#	cat(as.character(frm),":",nrow(testData),":",predictType,":",class(object),"\n");
	if (is.null(testData))
	{
		stop("No test data");
	}
	out <-numeric(nrow(testData))

	
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			pobj <- object;
			switch(predictType[1], 
				linear = 
					{		
						out <- predict(pobj,testData,type = 'lp');
					},
				prob = 
					{
						out <- 1.0/(1.0+exp(-predict(pobj,testData,type = 'lp')));
					},
					{
						out <- predict(pobj,testData,type = 'lp');
					}
			)
		},
		fi =
		{
			frm <- formula(terms(object));
			cf <- object$estimations;
			if (object$type=="COX")
			{			
				mm <- model.matrix(frm, testData)
			}
			else
			{
				mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
				mm[,1] <- rep(1,nrow(mm));
			}
			if ((any(is.na(mm))) && impute)
			{
				mm[,-1]=nearestNeighborImpute(tobeimputed=mm[,-1],referenceSet = object$model[,-1]);
			}
			if ((predictType[1] == "bagg") && (!is.null(object$coefficients.List)))
			{
				pred <- NULL;
				if (object$type == "LOGIT") 
				{
					pred <- .Call("predictForFresa",cf,mm,"prob",object$type);
				}
				if (length(object$coefficients.List) > 0)
				{
#					cat("Bagging");
					premin <- 2.0*object$fraction;
					premax <- 1.0 - premin;
					mm <- model.matrix(formula(object$formula.List[[1]]), testData);
					maout <- as.vector(mm  %*% object$coefficients.List[[1]]);
					npout <- maout;
					avgout <- maout;
#					print(object$coefficients.List);
					totwts <- sum(object$wts.List);
					wts <- object$wts.List[1]/totwts;
					fraction <- numeric(nrow(testData));
					if (object$type == "LOGIT") 
					{
						fraction <- fraction + wts*(maout > 0.0);
						tgext <- (maout > 36.0);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- 36.0;
						}
						tgext <- (maout < -36.0);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- -36.0;
						}
						maout <- 1.0/(1.0+exp(-maout));
						tgext <- (maout > premax);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- premax;
						}
						tgext <- (maout < premin);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- premin;
						}
						avgout <- maout*wts;
						npout <- log(1.0-maout)*wts;
						maout <- log(maout)*wts; 
					}
					else
					{
						maout <- maout*wts;
					}
					if (length(object$coefficients.List) > 1)
					{
						for (nm in 2:length(object$coefficients.List))
						{
							mm <- model.matrix(formula(object$formula.List[[nm]]), testData)
							nout <- as.vector(mm  %*% object$coefficients.List[[nm]]);
							wts <- object$wts.List[nm]/totwts;
							fraction <- fraction + wts*(nout > 0.0);
							if (object$type == "LOGIT") 
							{
								tgext <- (nout > 36.0);
								if (sum(tgext) > 0 ) 
								{
									nout[tgext] <- 36.0;
								}
								tgext <- (nout < -36.0);
								if (sum(tgext) > 0 ) 
								{
									nout[tgext] <- -36.0;
								}
								nout <- 1.0/(1.0+exp(-nout));
								tgext <- (nout > premax);
								if (sum(tgext) > 0 ) 
								{
									nout[tgext] <- premax;
								}
								tgext <- (nout < premin);
								if (sum(tgext) > 0 ) 
								{
									nout[tgext] <- premin;
								}
								avgout <- avgout+nout*wts;
								maout <- maout+log(nout)*wts;
								npout <- npout+log(1.0-nout)*wts;
							}
							else
							{
								maout <- maout + wts*nout;
							}
						}
					}
#					out <- maout;
					if (object$type == "LOGIT") 
					{
						tgext <- (maout > 36.0);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- 36.0;
						}
						tgext <- (maout < -36.0);
						if (sum(tgext) > 0 ) 
						{
							maout[tgext] <- -36.0;
						}
						tgext <- (npout > 36.0);
						if (sum(tgext) > 0 ) 
						{
							npout[tgext] <- 36.0;
						}
						tgext <- (npout < -36.0);
						if (sum(tgext) > 0 ) 
						{
							npout[tgext] <- -36.0;
						}
						nout <- exp(maout)/(exp(maout)+exp(npout));
						out <- avgout;
						if (length(object$coefficients.List) >= 10)
						{
#							selcheck <- ((avgout > 0.5) & (fraction < 0.5) ) | ((avgout < 0.5) & (fraction > 0.5)); 
#							out[selcheck] <- 0.5*(avgout[selcheck]+fraction[selcheck]);
							out <- 0.5*(avgout+fraction);
						}
						attr(out,"WeightedAverage") <- avgout;
						attr(out,"wLOP") <- nout;
						attr(out,"LOGIT") <- as.vector(pred$prediction);
						attr(out,"Fraction") <- fraction;
#						print(out);
					}
					else
					{
						out <- maout;
					}
				}
				else
				{
					mm <- model.matrix(frm, testData)
					out <- as.vector(mm  %*% cf);
					if (object$type=="LOGIT") 
					{
						tgext <- (out > 36.0);
						if (sum(tgext) > 0 ) 
						{
							out[tgext] <- 36.0;
						}
						tgext <- (out < -36.0);
						if (sum(tgext) > 0 ) 
						{
							out[tgext] <- -36.0;
						}
						out <- 1.0/(1.0+exp(-out));
					}
				}
			}
			else
			{
				pred <-.Call("predictForFresa",cf,mm,predictType[1],object$type);
				out <- as.vector(pred$prediction);
			}
#			names(out) <- rownames(testData);
			if (length(out)==nrow(testData)) names(out) <- rownames(testData);
		},
		tr =
		{
			frm <- formula(terms(object));
			if (!is.null(frm))
			{
				cf <- object$estimations;
				if (object$type=="COX")
				{			
					mm <- model.matrix(frm, testData)
				}
				else
				{
					mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
					mm[,1] <- rep(1,nrow(mm));
				}
				if (!is.null(cf))
				{
					if (length(cf)>0)
					{
						warning("Warning: Fitting error. Object: ",class(object),"NA coff set to zero\n");
		#				print(cf);
						pred <-.Call("predictForFresa",cf,mm,predictType[1],object$type);
						out <- as.vector(pred$prediction);
					}
					else
					{
						out <- rep(NA,nrow(testData));
					}
				}
				else
				{
					warning("Warning: Fitting error. Object: ",class(object),"All predictions set to NA \n");
					out <- rep(NA,nrow(testData));
				}
			}
			else
			{
				out <- rep(NA,nrow(testData));
			}
#			names(out) <- rownames(testData);
			if (length(out)==nrow(testData)) names(out) <- rownames(testData);
		},
		sv =
		{
			out <- predict(object,testData);
		},
		or =
		{
			totClases <- length(object$theScores);
			totModels <- totClases-1;
			theScores <- as.numeric(object$theScores);

			outOrder <- matrix(nrow=nrow(testData),ncol=totModels);
			stepOrder <- matrix(nrow=nrow(testData),ncol=totModels);
			outClass <- matrix(nrow=nrow(testData),ncol=totClases);
			out <- matrix(nrow=nrow(testData),ncol=7+totClases);


			colc <- 1;
			for (s in theScores)
			{
				outClass[,colc] <- predict(object$theClassBaggs[[colc]]$bagged.model,testData,"prob");
				if (colc <= totModels)
				{
					outOrder[,colc] <- predict(object$theBaggedModels[[colc]]$bagged.model,testData,"prob");
					stepOrder[,colc] <- predict(object$redBaggedModels[[colc]]$bagged.model,testData,"prob");
				}
				colc <- colc + 1;
			}
#			print(outClass[1,])
#			print(outOrder[1,])
#			print(stepOrder[1,])
			uscores <- 1:totClases;
			for (n in 1:nrow(testData))
			{
				pClass <- outClass[n,];

				pOrder <- rep(1.0,totClases);
				for (pc in 1:totClases)
				{
					 for (p in 1:totModels)
					 {
						 if (p >= pc)
						 {
							 pOrder[pc] <- pOrder[pc]*(1.0-outOrder[n,p]);
						 }
						 else
						 {
							 pOrder[pc] <- pOrder[pc]*(outOrder[n,p]);
						 }
					 }
				}
				pOrder <- pOrder/sum(pOrder);

				pCond <- rep(1.0,totClases);
				pCond[1] <- (1.0-stepOrder[n,1]) > 0.5;
				pCond[totClases] <- stepOrder[n,totModels] > 0.5;
				for (p in 2:totModels)
				{
					pCond[p] <- 0.5*(((1.0-stepOrder[n,p]) > 0.5) + (stepOrder[n,p-1] > 0.5));
				}

				maxOrder <- max(pOrder);
				whoOrder <- pOrder > 0.999*maxOrder;
				wts <- pOrder[whoOrder] - 0.999*maxOrder;
				wts <- wts*wts;
				scoreOrder <- sum(uscores[whoOrder]*wts)/sum(wts);

				maxClass <- max(pClass)
				whoClass <- pClass > 0.999*maxClass;
				wts <- pClass[whoClass] - 0.999*maxClass;
				wts <- wts*wts;
				scoreClass <- sum(uscores[whoClass]*wts)/sum(wts);

				pComb <- (pOrder+1.0e-5)*(pCond+1.0e-5)*(pClass+1.0e-5);
				pComb <- pComb/sum(pComb);
				maxComb <- max(pComb);
				whoComb <- pComb > 0.999*maxComb;
				wts <- pComb[whoComb] - 0.999*maxComb;
				wts <- wts*wts;
				scoreComb <- sum(uscores[whoComb]*wts)/sum(wts);
	
				iscoreTotal <- as.integer(scoreComb + 0.5);
				if (iscoreTotal == 1)
				{
					sadj <- stepOrder[n,1];
				}
				else
				{
					if (iscoreTotal == totClases)
					{
						sadj <- stepOrder[n,totModels] - 1.0;
					}
					else
					{
						sadj <- stepOrder[n,iscoreTotal-1] - 1.0 + stepOrder[n,iscoreTotal];
					}
				}
				iscoreTotal <- as.integer(iscoreTotal + sadj + 0.5);

				out[n,1] <- theScores[iscoreTotal];
				out[n,2] <- sadj;
				out[n,3] <- scoreOrder;
				out[n,4] <- maxOrder;
				out[n,5] <- scoreClass;
				out[n,6] <- maxClass;
				out[n,7] <- scoreComb;
				out[n,8:(7+totClases)] <- pOrder;
			}
			rownames(out) <- rownames(testData);
			colnames(out) <- c("Class","Adjs","O.Class","p.Ord","C.Class","p.Class","T.Class",theScores);
		},
		BS =
		{
			if (is.null(object$oridinalModels))
			{
				if (is.null(object$bagging) || is.null(object$bagging$bagged.model))
				{
					out <- predict(object$forward.model$final.model,...);
					attr(out,"model") <- "forward.update";
				}
				else
				{
					# if (object$equivalent)
					# {
						# if (length(object$formula.list>1))
						# {
							# if (object$bagging$bagged.model$type == "LOGIT")
							# {
								# pred <- ensemblePredict(object$formula.list,
								# object$BSWiMS.model$bootCV$data,
								# testData,"prob",
								# object$bagging$bagged.model$type)
							# }
							# else
							# {
								# pred <- ensemblePredict(object$formula.list,
								# object$BSWiMS.model$bootCV$data,
								# testData,predictType,
								# object$bagging$bagged.model$type)
							# }						
							# out <- as.numeric(pred$ensemblePredict)
	 #						names(out) <- rownames(testData);
							# attr(out,"MeanEnsemblePredict") <- pred$wPredict;
							# attr(out,"model") <- "Ensemble";
						# }
						# else
						# {
							# out <- predict(object$bagging$bagged.model,...)
						# }
						# if (length(out)==nrow(testData)) names(out) <- rownames(testData);
					# }
					# else
					{
						out <- predict(object$bagging$bagged.model,...);
#						out <- predict(object$bagging,...);
						attr(out,"model") <- "bagged";
						if (object$bagPredictType[1]=="wNN")
						{
							if (!is.null(object$bagging$nnmodel))
							{
								out <- predict(object$bagging,...);
								attr(out,"model") <- "wNN";
							}
						}
						if (object$bagPredictType[1]=="Ens")
						{
							if (!is.null(object$bagging$nnmodel))
							{
								out <- (out + predict(object$bagging,...))/2.0;
								attr(out,"model") <- "Ens";
							}
						}
					}
				}
			}
			else
			{
				out <- predict(object$oridinalModels,...)[,1];
				attr(out,"model") <- "ordinal";
			}
		},
		{
			cf <- coef(object)
			if (is.null(cf))
			{
				warning("Warning: Null object in predict. Object: ",cobj,"\n");
				out <- rep(NA,nrow(testData));
			}
			else
			{
				frm <- formula(terms(object));
				out <- rep(NA,nrow(testData));
				if (!is.null(frm))
				{
					cf <- object$coefficients;
					if (length(cf)>0)
					{
						s <- is.na(cf);
						if (any(s)) 
						{
							cf[s] <- 0;
						}
						if (object$type=="COX")
						{			
							mm <- model.matrix(frm, testData)
						}
						else
						{
							mm <- as.matrix(model.frame(frm, testData,na.action=na.pass))
							mm[,1] <- rep(1,nrow(mm));
						}
						if (!is.null(mm))
						{
#							if ((class(mm) == "matrix") && (ncol(mm)>0) &&  (class(cf) == "vector"))
							{
								switch(predictType[1], 
									linear = 
										{
										   out <- as.vector(mm  %*% cf);
										}, 
									prob = 
										{
										  out <- as.vector(mm  %*% cf);
										  out <- 1.0/(1.0+exp(-out));
										},
									pbagg =
										{
											if (!is.null(object$coefficients.List))
											{
												if (length(object$coefficients.List)>0)
												{
													mm <- model.matrix(formula(object$formula.List[[1]]), testData)
													maout <- as.vector(mm  %*% object$coefficients.List[[1]]);
													maout <- 1.0/(1.0+exp(-maout));
#													miout <- maout;
													if (length(object$coefficients.List)>1)
													{
														for (nm in 2:length(object$coefficients.List))
														{
															mm <- model.matrix(formula(object$formula.List[[nm]]), testData)
															nout <- as.vector(mm  %*% object$coefficients.List[[nm]]);
															nout <- 1.0/(1.0+exp(-nout));
															maout <- maout + nout;
#															maout <- pmax(maout,nout);
#															miout <- pmin(miout,nout);
#															print(c(maout[1],miout[1],nout[1]));
														}
													}
													out <- maout/length(object$coefficients.List);
#													whodif <- ((1.0-miout) > out);
#													out[whodif] <- miout[whodif];
#													print(c(maout[1],miout[1],out[1]));
												}
												else
												{
													  out <- as.vector(mm  %*% cf);
													  out <- 1.0/(1.0+exp(-out));
												}
											}
											else
											{
												  out <- as.vector(mm  %*% cf);
												  out <- 1.0/(1.0+exp(-out));
											}
											
										},
										{
										  out <- as.vector(mm  %*% cf);
										}
								)
							}
						}
					}
				}
				else
				{
					warning(paste(as.character(match.call()),"No formula \n"));
				}
			}
			if (length(out)==nrow(testData)) names(out) <- rownames(testData);
		}
	)
	s <- is.na(out);
	if (any(s)) 
	{
		warning(paste(as.character(match.call()),"Warning NA predict.fitFRESA \n"));
		switch(predictType[1], 
			linear = 
				{
				   out[s] <- 0;
				}, 
			prob = 
				{
				  out[s] <- 0.5;
				}
		)
	}
	if (length(out)!=nrow(testData))
	{
		warning(object$type,": Different number of rows:",length(out),"(",nrow(testData),"). Setting to NA \n");
		tout <- out;
		out=rep(NA,nrow(testData));
		names(out) <- rownames(testData);
		out[names(tout)] <- tout;
	}
    return (out)
}
