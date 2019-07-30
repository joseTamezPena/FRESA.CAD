uniRankVar <-
function(variableList,formula,Outcome,data,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail","Tail","RawRaw"),type=c("LOGIT","LM","COX"),rankingTest=c("zIDI","zNRI","IDI","NRI","NeRI","Ztest","AUC","CStat","Kendall"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,description=".",uniType=c("Binary","Regression"),FullAnalysis=TRUE,acovariates=NULL,timeOutcome=NULL) 
{

	if (is.null(timeOutcome)) timeOutcome="";
	if (is.na(timeOutcome)) timeOutcome="";
#	cat(categorizationType,"Catego \n")
#	cat(uniType," Uni \n")

	if (is.null(raw.dataFrame))  raw.dataFrame <- data;
	type <- match.arg(type);
	uniType <- match.arg(uniType);
	categorizationType <- match.arg(categorizationType);
#	cat(categorizationType,"\n")
	rankingTest <- match.arg(rankingTest);
	colnamesList <- as.vector(variableList[,1]);
	varinserted <- length(colnamesList);

	Name <- colnamesList;
	parent <- colnamesList;
	descripList <- colnamesList;
	
	if (FullAnalysis || (categorizationType != "Raw"))
	{
		varinserted <- 0;

		if (description == ".")
		{
			descripList <- colnamesList;
		}
		else
		{
			descripList <- as.vector(variableList[,description]);
		}
		
		if (uniType=="Binary")
		{
			caserawsample <- subset(raw.dataFrame,get(Outcome)  == 1);
			controlrawsample <- subset(raw.dataFrame,get(Outcome) == 0);

			caseZsample <- subset(data,get(Outcome)  == 1);
			controlZsample <- subset(data,get(Outcome) == 0);

			sizecaseZsample <- nrow(caseZsample);
			sizecontrolZsample <- nrow(controlZsample);
		}
		else
		{
			caserawsample <- NULL;
			controlrawsample <- NULL;

			caseZsample <- NULL;
			controlZsample <- NULL;

			sizecaseZsample <- NULL;
			sizecontrolZsample <- NULL;
		}

		size = length(colnamesList);
		maxsize = size*(length(cateGroups)+1);
		if (categorizationType=="RawRaw") maxsize = size*(as.integer(size/2)+1)
		Name <- character(maxsize);

		parent <- character(maxsize);
		descrip <- character(maxsize);


		IDI <- numeric(maxsize);
		NRI <- numeric(maxsize);
		zIDI <- numeric(maxsize);
		zNRI <- numeric(maxsize);
		ROCAUC <- numeric(maxsize);
		ZGLM <- numeric(maxsize);
		Beta <- numeric(maxsize);

		NeRI <- numeric(maxsize);
		BinRes.p <- numeric(maxsize);
		WilcoxRes.p <- numeric(maxsize);
		TstudentRes.p <- numeric(maxsize);
		FRes.p <- numeric(maxsize);


		cohortMean <- numeric(maxsize);
		cohortStd <- numeric(maxsize);
		cohortKSD <- numeric(maxsize);
		cohortKSP <- numeric(maxsize);
		cohortZKSP <- numeric(maxsize);
		cohortZKSD <- numeric(maxsize);

		
		caseMean <- numeric(maxsize);
		caseStd <- numeric(maxsize);
		caseKSD <- numeric(maxsize);
		caseKSP <- numeric(maxsize);
		caseZKSP <- numeric(maxsize);
		caseZKSD <- numeric(maxsize);
		caseN_Z_Low_Tail <- numeric(maxsize);
		caseN_Z_Hi_Tail <- numeric(maxsize);

		controlMean <- numeric(maxsize);
		controlStd <- numeric(maxsize);
		controlKSD <- numeric(maxsize);
		controlKSP <- numeric(maxsize);
		controlZKSP <- numeric(maxsize);
		controlZKSD <- numeric(maxsize);
		controlN_Z_Low_Tail <- numeric(maxsize);
		controlN_Z_Hi_Tail <- numeric(maxsize);
		kendall.r <- numeric(maxsize);
		spearman.r <- numeric(maxsize);
		pearson.r <- numeric(maxsize);
		kendall.p <- numeric(maxsize);
		cStatCorr <- numeric(maxsize);
		
		t.Rawvalue <- numeric(maxsize);
		t.Zvalue <- numeric(maxsize);
		wilcox.Zvalue <- numeric(maxsize);
		Sensitivity <- numeric(maxsize);
		Specificity <- numeric(maxsize);

		
		frm1 <- formula;
		ftmp <- formula(frm1);

		varlist <- attr(terms(ftmp),"variables")
		dependent <- as.character(varlist[[2]])
		if (length(dependent)==3)
		{
			type = "COX"
			timeOutcome = dependent[2];
			baseData <- as.data.frame(data[,c(timeOutcome,Outcome)]);
			colnames(baseData) <- c(timeOutcome,Outcome);
		}
		else
		{
			baseData <- as.data.frame(data[[Outcome]]);
			colnames(baseData) <- c(Outcome);
		}
		if (length(varlist)>2)
		{
			for (i in 3:length(varlist))
			{
				bnames <- colnames(baseData)
				baseData <- cbind(baseData,data[[as.character(varlist[[i]])]])
				colnames(baseData) <- c(bnames,as.character(varlist[[i]]));
			}
		}


		if (type=="COX")
		{
			termslist <- attr(terms(ftmp),"term.labels");
			if (length(termslist)==0)
			{
				ftmp <- formula(paste(frm1,"+ dummy"));
#				cat("Dummy Cox: ",frm1,": ",paste(frm1,"+ dummy"),"\n");
				dummy <-  rnorm(nrow(baseData));
				coxframe <- cbind(baseData,dummy);
				colnames(coxframe) <- c(colnames(baseData),"dummy")
				
				bmodel <- modelFitting(ftmp,coxframe,type,fitFRESA=FALSE)
				basepredict <- predict.fitFRESA(bmodel,coxframe, 'prob');
				baseResiduals <- residualForFRESA(bmodel,coxframe,Outcome);
			}
			else
			{
				bmodel <- modelFitting(ftmp,data,type,fitFRESA=FALSE)
				baseResiduals <- residualForFRESA(bmodel,data,Outcome)+rnorm(nrow(data),0,1e-10);
				basepredict <- predict.fitFRESA(bmodel,data, 'prob');
			}
		}
		else
		{
			bmodel <- modelFitting(ftmp,data,type,fitFRESA=FALSE)
			baseResiduals <- residualForFRESA(bmodel,data,Outcome)+rnorm(nrow(data),0,1e-10);
			basepredict <- predict.fitFRESA(bmodel,data, 'prob');
		}
#		cat("Start Ranking\n");
		
	#	if (!FullAnalysis)
	#	{
	#		rankingTest="Ztest";
	#	}
		
		casesOutcome <- (data[,Outcome]==1)
		controlOutcome <- (data[,Outcome]==0)
		for (j in 1:size)
		{
			if ((j %% 100)==0) cat (j,":",colnamesList[j]," ");
			if ((j %% 500)==0) cat ("\n");
			frm1 = formula;
			categories = 1;
			catlist <- character(1);
			caseCount1 = 0;
			caseCount2 = 0;
			controlCount1 = 0;
			controlCount2 = 0;
			stddf=0;
			kendcor <- NA;
			pearcor <- NA;
			speacor <- NA;
			cstat <- NA;
			kstZdf <- NA;
			kstdf <- NA;
			meCa <- NA;
			stdCa <- NA;
			kstCa <- NA;
			kstZCa <- NA;
			meCo <- NA;
			stdCo <- NA;
			kstCo <- NA;
			kstZCo <- NA;
			rtt <- NA;
			ztt <- NA;
			medf <- NA;
			datacolumn <- data[,colnamesList[j]];
			dataoutcome <- data[,Outcome];

			if ((uniType=="Binary")&& FullAnalysis)
			{
				wtt <- -qnorm(wilcox.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude)$p.value,0,1);
			}
			else
			{
				wtt <- NA;
			}
	#		cat (colnamesList[j],"Wilcox: ",wtt,"\n")
			if (FullAnalysis)
			{
				stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
				if (stddf>0)
				{
					kendcor <- try(cor.test(datacolumn,dataoutcome,method="kendall",na.action=na.omit));
					pearcor <- try(cor.test(datacolumn,dataoutcome,method="pearson",na.action=na.omit));
					speacor <- try(cor.test(datacolumn,dataoutcome,method="spearman",na.action=na.omit));
					if (inherits(kendcor, "try-error")) kendcor$estimate = 0;
					if (is.na(kendcor$estimate)) kendcor$estimate = 0;
					if (kendcor$estimate > 0)
					{
						cstat <- rcorr.cens(datacolumn,dataoutcome, outx=FALSE)
					}
					else
					{
						cstat <- rcorr.cens(datacolumn,-dataoutcome, outx=FALSE)
					}
				}
			}
	#		cat (colnamesList[j],"cstat: ",cstat,"\n")
			datalength <- length(table(datacolumn));
			if (datalength>2)
			{

				if (FullAnalysis)
				{
					medf <- mean(datacolumn,na.rm = TRUE);
					stddf <- sd(datacolumn,na.rm = TRUE);
					kstZdf <- ks.test(datacolumn,"pnorm",medf,stddf);
					medf <- mean(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
					stddf <- sd(raw.dataFrame[,colnamesList[j]],na.rm = TRUE);
					kstdf <- ks.test(raw.dataFrame[,colnamesList[j]],"pnorm",medf,stddf);
				}
				if (uniType=="Binary")
				{
					if (FullAnalysis)
					{
						meCa <- mean(caseZsample[,colnamesList[j]],na.rm = TRUE);
						stdCa <- sd(caseZsample[,colnamesList[j]],na.rm = TRUE);
						kstZCa <- ks.test(caseZsample[,colnamesList[j]],"pnorm",meCa,stdCa);
						
						meCo <- mean(controlZsample[,colnamesList[j]],na.rm = TRUE);
						stdCo <- sd(controlZsample[,colnamesList[j]],na.rm = TRUE);
						kstZCo <- ks.test(controlZsample[,colnamesList[j]],"pnorm",meCo,stdCo);
						
						meCa <- mean(caserawsample[,colnamesList[j]],na.rm = TRUE);
						stdCa <- sd(caserawsample[,colnamesList[j]],na.rm = TRUE);
						kstCa <- ks.test(caserawsample[,colnamesList[j]],"pnorm",meCa,stdCa);
						
						meCo <- mean(controlrawsample[,colnamesList[j]],na.rm = TRUE);
						stdCo <- sd(controlrawsample[,colnamesList[j]],na.rm = TRUE);
						kstCo <- ks.test(controlrawsample[,colnamesList[j]],"pnorm",meCo,stdCo);
						
						rtt <- try(t.test(controlrawsample[,colnamesList[j]],caserawsample[,colnamesList[j]],na.action=na.exclude));
						ztt <- try(t.test(controlZsample[,colnamesList[j]],caseZsample[,colnamesList[j]],na.action=na.exclude));
					}

					if (!is.na(cateGroups[1]))
					{
						if (categorizationType != "Raw")
						{
							zthr = sprintf("[,'%s'] < %5.3f )",colnamesList[j],qnorm(cateGroups[1]));
							
#							cat(zthr,"\n")

							caseCount1 <- eval(parse(text = paste("sum(caseZsample",zthr)));
							controlCount1 <- eval(parse(text = paste("sum(controlZsample",zthr)));

							categories=length(cateGroups);
							if (!is.na(cateGroups[categories]))
							{
								zthr = sprintf("[,'%s'] > %5.3f )",colnamesList[j],qnorm(cateGroups[categories]));
	#						cat(zthr,"\n")
								caseCount2 <- eval(parse(text = paste("sum(caseZsample",zthr)));
								controlCount2 <- eval(parse(text = paste("sum(controlZsample",zthr)));

							}
							else
							{
								zthr = sprintf("[,'%s'] > %5.3f )",colnamesList[j],1-qnorm(cateGroups[1]));
	#						cat(zthr,"\n")
								caseCount2 <- eval(parse(text = paste("sum(caseZsample",zthr)));
								controlCount2 <- eval(parse(text = paste("sum(controlZsample",zthr)));
							}
						}
					}
				}


				if (((min(datacolumn)<0) && (datalength>4)) || (categorizationType=="RawRaw"))
				{
#					if (j<10) print(colnamesList[j]);
					catlist <- character(1);
					switch(categorizationType,
						Raw =
						{
							categories=1;
							catlist[1] <- colnamesList[j];				
						},
						RawRaw =
						{
							catlist[1] <- colnamesList[j];
							for (jj in j:size)
							{
								catlist <- append(catlist,sprintf("I(%s*%s)",colnamesList[j],colnamesList[jj]));
							}							
							categories = length(catlist);
#							if (j<5) print(catlist);
						},
						Categorical =
						{
							categories=length(cateGroups);

							zthr = sprintf("%5.3f",qnorm(cateGroups[1]));


							for (n in 1:categories)
							{
								if (n==1)
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,")");
									catlist[1] <- catvar;
								}
								else
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
									zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I((",colnamesList[j]);
									catvar = paste(catvar," >= ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,") & (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr2);
									catvar = paste(catvar,"))");
									catlist <- append(catlist,catvar);
								}
							}
							zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,")");
							catlist <- append(catlist,catvar);
							categories = categories + 1;
						},
						ZCategorical =
						{
							categories=length(cateGroups);
							zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

							for (n in 1:categories)
							{

								if (n==1)
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,"))");
									catlist[1] <- catvar;
								}
								else
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
									zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* ((");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," >= ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,") & (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr2);
									catvar = paste(catvar,")))");
									catlist <- append(catlist,catvar);
								}
							}
							zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,"))");
							catlist <- append(catlist,catvar);
							categories = categories+1;
						},
						RawZCategorical =
						{
							categories=length(cateGroups);

							zthr = sprintf("%5.3f",qnorm(cateGroups[1]));

							catlist[1] <- colnamesList[j];				
							for (n in 1:categories)
							{

								if (n==1)
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,"))");
									catlist <- append(catlist,catvar);
								}
								else
								{
									zthr = sprintf("%5.3f",qnorm(cateGroups[n-1]));
									zthr2 = sprintf("%5.3f",qnorm(cateGroups[n]));
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* ((");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," >= ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,") & (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr2);
									catvar = paste(catvar,")))");
									catlist <- append(catlist,catvar);
								}
							}
							zthr = sprintf("%5.3f",qnorm(cateGroups[categories]));
							catvar = paste("I(",colnamesList[j]);
							catvar = paste(catvar,"* (");
							catvar = paste(catvar,colnamesList[j]);
							catvar = paste(catvar," >= ");
							catvar = paste(catvar,zthr);
							catvar = paste(catvar,"))");
							catlist <- append(catlist,catvar);
							categories = categories+2;
						},			
						RawZTail =
						{
							categories = 1;
							catlist[1] <- colnamesList[j];
#							print(catlist);
#							cat(colnamesList[j],"\n");					

							if (!is.null(sizecaseZsample))
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
								f1= caseCount1/sizecaseZsample;
								f2= controlCount1/sizecontrolZsample;
								if ((f1>f2)&&(f1>0.1))
								{
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,"))");
									catlist <- append(catlist,catvar);
									categories = categories+1;
								}				


								zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
								f1= caseCount2/sizecaseZsample;
								f2= controlCount2/sizecontrolZsample;
								if ((f1>f2)&&(f1>0.1))
								{
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar,"* (");
									catvar = paste(catvar,colnamesList[j]);
									catvar = paste(catvar," > ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,"))");
									catlist <- append(catlist,catvar);
									categories = categories+1;
								}
							}
						},			
						RawTail =
						{
							categories = 1;
							catlist[1] <- colnamesList[j];				
							if (!is.null(sizecaseZsample))
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
								f1 = caseCount1/sizecaseZsample;
								f2 = controlCount1/sizecontrolZsample; 
								if ((f1>f2)&&(f1>0.1))   # will add only if fraction is greater
								{
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar," < ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,")");
									catlist <- append(catlist,catvar);
									categories = categories+1;
								}				


								zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
								f1= caseCount2/sizecaseZsample;
								f2= controlCount2/sizecontrolZsample;
								if ((f1>f2)&&(f1>0.1)) # will add only if fraction is greater
								{
									catvar = paste("I(",colnamesList[j]);
									catvar = paste(catvar," > ");
									catvar = paste(catvar,zthr);
									catvar = paste(catvar,")");
									catlist <- append(catlist,catvar);
									categories = categories+1;
								}
							}
						},
						Tail =
						{
							categories = 0;
							if (!is.null(sizecaseZsample))
							{
								zthr = sprintf("%5.3f",qnorm(cateGroups[1]));
								if (caseCount1 > caseCount2)
								{
									f1= caseCount1/sizecaseZsample;
									f2= controlCount1/sizecontrolZsample;
									if (f1>0.00)
									{
										catvar = paste("I(",colnamesList[j]);
										catvar = paste(catvar,"* (");
										catvar = paste(catvar,colnamesList[j]);
										catvar = paste(catvar," < ");
										catvar = paste(catvar,zthr);
										catvar = paste(catvar,"))");
										catlist[1] <- catvar;
										categories = categories+1;
									}
								}
								else
								{
									zthr = sprintf("%5.3f",qnorm(1.0-cateGroups[1]));
									f1= caseCount2/sizecaseZsample;
									f2= controlCount2/sizecontrolZsample;
									if (f1>0.00)
									{
										catvar = paste("I(",colnamesList[j]);
										catvar = paste(catvar,"* (");
										catvar = paste(catvar,colnamesList[j]);
										catvar = paste(catvar," > ");
										catvar = paste(catvar,zthr);
										catvar = paste(catvar,"))");
										catlist[1] <- catvar;
										categories = categories+1;
									}
								}							
							}
						},			
						{
							categories = 1;
							catlist[1] <- colnamesList[j];
						}
					)
				}
				else
				{
					categories=1;
					catlist[1] <- colnamesList[j];
				}
			}
			else
			{
				categories=1;
				catlist[1] <- colnamesList[j];
				if (FullAnalysis)
				{
					medf = table(datacolumn)[1];			
					stddf = table(datacolumn)[2];			
					if (uniType=="Binary")
					{
						meCa <- table(caseZsample[,colnamesList[j]])[1];
						stdCa <- table(caseZsample[,colnamesList[j]])[2];
						
						meCo <- table(controlZsample[,colnamesList[j]])[1];
						stdCo <- table(controlZsample[,colnamesList[j]])[2];
					}
				}
			}
			if (categories>0)
			{
#				cat(catlist[1],"\n");
				if (categorizationType != "RawRaw")
				{
					mdata <- cbind(baseData,datacolumn)
					colnames(mdata) <- c(colnames(baseData),colnamesList[j]);
				}
				else
				{
					mdata <- data;
				}
				
				for (n in 1:categories)
				{
					termName <- str_replace_all(catlist[n]," ","");
					termName <- str_replace_all(termName,"<"," < ");
					termName <- str_replace_all(termName,">"," > ");
					termName <- str_replace_all(termName,"&"," & ");
					termName <- str_replace_all(termName,"=","= ");
					termName <- str_replace_all(termName,fixed("> ="),">=");
					termName <- str_replace_all(termName,fixed("*")," * ");
					frmg <- paste( formula,paste("+",termName));
#					if (j<5) cat(catlist[n],":",frmg,"\n")
					ftmg <- formula(frmg);
					if (type=="COX") 
					{
						zcol=4;
					}
					else
					{
						zcol=3;
					}
					
					lmodel <- modelFitting(ftmg,mdata,type,fitFRESA=FALSE)
					sen=0;
					spe=0;
					if (!inherits(lmodel, "try-error"))
					{
						modcoef <- summary(lmodel)$coefficients;
						sizecoef <- length(lmodel$coef);
						if ((uniType=="Binary") && FullAnalysis) 
						{
							sen = sum( 1.0*(lmodel$linear.predictors>0) & casesOutcome )/sizecaseZsample;
							spe = sum( 1.0*(lmodel$linear.predictors<0) & controlOutcome )/sizecontrolZsample;
						}
					}
					else
					{
						modcoef <- NULL;
						sizecoef <- NULL;
					}
						varinserted <- varinserted+1;
						Name[varinserted] <- termName;
						parent[varinserted] <- colnamesList[j];
						descrip[varinserted] <- descripList[j];
						if ((uniType=="Binary")&&(FullAnalysis))
						{
							caseMean[varinserted] <- meCa; 
							caseStd[varinserted] <- stdCa;
							if (!is.na(kstCa[[1]])) 
							{
								caseKSD[varinserted] <- kstCa$statistic;
								caseKSP[varinserted] <- kstCa$p.value;
								caseZKSP[varinserted] <- kstZCa$p.value;
								caseZKSD[varinserted] <- kstZCa$statistic;
							}
							else
							{
								caseKSD[varinserted] <- NA;
								caseKSP[varinserted] <- NA;
								caseZKSP[varinserted] <- NA;
								caseZKSD[varinserted] <- NA;
							}
							controlMean[varinserted] <- meCo; 
							controlStd[varinserted] <- stdCo;
							if (!is.na(kstCo[[1]]))
							{
								controlKSD[varinserted] <- kstCo$statistic;
								controlKSP[varinserted] <- kstCo$p.value;
								controlZKSP[varinserted] <- kstZCo$p.value;
								controlZKSD[varinserted] <- kstZCo$statistic;
							}
							else
							{
								controlKSD[varinserted] <- NA;
								controlKSP[varinserted] <- NA;
								controlZKSP[varinserted] <- NA;
								controlZKSD[varinserted] <- NA;
							}
							if (!is.na(rtt[[1]]))
							{
								if ( !inherits(rtt, "try-error"))
								{
									t.Rawvalue[varinserted] <- rtt$statistic;
								}
								else
								{
									t.Rawvalue[varinserted] <- NA;
								}
								if ( !inherits(ztt, "try-error"))
								{
									t.Zvalue[varinserted] <- ztt$statistic;
								}
								else
								{
									t.Zvalue[varinserted] <- NA;
								}
							}
							else
							{
								t.Rawvalue[varinserted] <- NA;
								t.Zvalue[varinserted] <- NA;
							}
							if (!is.na(wtt[[1]]))
							{
								wilcox.Zvalue[varinserted] <- wtt;
							}
							else
							{
								wilcox.Zvalue[varinserted] <- NA;
							}
						}
					
					if (FullAnalysis)
					{

						cohortMean[varinserted] <- medf; 
						cohortStd[varinserted] <- stddf;
						if (!is.na(kstdf[[1]]))
						{
							cohortKSD[varinserted] <- kstdf$statistic;
							cohortKSP[varinserted] <- kstdf$p.value;
							cohortZKSP[varinserted] <- kstZdf$p.value;
							cohortZKSD[varinserted] <- kstZdf$statistic;
						}
						else
						{
							cohortKSD[varinserted] <- NA;
							cohortKSP[varinserted] <- NA;
							cohortZKSP[varinserted] <- NA;
							cohortZKSD[varinserted] <- NA;
						}
						
						if (!is.na(kendcor[[1]])) 
						{
							kendall.r[varinserted] <- kendcor$estimate;
							kendall.p[varinserted] <- kendcor$p.value;
							pearson.r[varinserted] <- pearcor$estimate;
							spearman.r[varinserted] <- speacor$estimate;
							cStatCorr[varinserted] <- cstat[1];
						}
						else
						{
							kendall.r[varinserted] <- NA;
							kendall.p[varinserted] <- NA;
							pearson.r[varinserted] <- NA;
							spearman.r[varinserted] <- NA;
							cStatCorr[varinserted] <- NA;
						}
					}
					
					if (is.null(sizecoef) || is.na(lmodel$coef[sizecoef])) 
					{
						if (FullAnalysis)
						{
							test=NA;
							if (uniType=="Binary")
							{
								IDI[varinserted] <- test;
								NRI[varinserted] <- test;
								zIDI[varinserted] <- test;
								zNRI[varinserted] <- test;
								ROCAUC[varinserted] <- test;
								caseN_Z_Low_Tail[varinserted] <- test;
								caseN_Z_Hi_Tail[varinserted] <- test;
								controlN_Z_Low_Tail[varinserted] <- test;
								controlN_Z_Hi_Tail[varinserted] <- test;
								Sensitivity[varinserted] <- sen;
								Specificity[varinserted] <- spe;
							}
							ZGLM[varinserted] <- test;
							NeRI[varinserted] <- test;
							BinRes.p[varinserted] <- test;
							WilcoxRes.p[varinserted] <- test;
							TstudentRes.p[varinserted] <- test;
							FRes.p[varinserted] <- test;
						}
					}
					else
					{
						if (FullAnalysis)
						{
							if (uniType=="Binary")
							{
								Sensitivity[varinserted] <- sen;
								Specificity[varinserted] <- spe;
								spredict <- predict.fitFRESA(lmodel,data, 'prob');
								iprob <- .Call("improveProbCpp",basepredict,spredict,dataoutcome);
								IDI[varinserted] <- iprob$idi;
								NRI[varinserted] <- iprob$nri;
								zIDI[varinserted] <- iprob$z.idi;
								zNRI[varinserted] <- iprob$z.nri;
								if (length(dataoutcome)==length(spredict))
								{
									ROCAUC[varinserted] <- pROC::roc( dataoutcome, spredict,plot=FALSE,auc=TRUE,quiet = TRUE)$auc[1];
								}
								else 
								{
									ROCAUC[varinserted] <- NA;
								}
								caseN_Z_Low_Tail[varinserted] <- caseCount1;
								caseN_Z_Hi_Tail[varinserted] <- caseCount2;
								controlN_Z_Low_Tail[varinserted] <- controlCount1;
								controlN_Z_Hi_Tail[varinserted] <- controlCount2;
							}

							varResiduals <- residualForFRESA(lmodel,data,Outcome);
							rprob <- .Call("improvedResidualsCpp",baseResiduals,varResiduals," ",0);
							NeRI[varinserted] <- rprob$NeRI;
							BinRes.p[varinserted] <- rprob$BinP.value;
							WilcoxRes.p[varinserted] <- rprob$WilcoxP.value;
							TstudentRes.p[varinserted] <- rprob$tP.value;
							FRes.p[varinserted] <- rprob$FP.value;
						}
						else
						{
							if (uniType=="Binary")
							{
#								spredict <- predict.fitFRESA(lmodel,mdata, 'prob');
								spredict <- predict.fitFRESA(lmodel,mdata, 'prob');
								iprob <- .Call("improveProbCpp",basepredict,spredict,dataoutcome);
								zIDI[varinserted] <- iprob$z.idi;
							}
						}
						Beta[varinserted] <- modcoef[sizecoef,1];
						ZGLM[varinserted] <- abs(modcoef[sizecoef,zcol]);
					}
				}
			}
		}
	}
#	cat(" Finished unitable \n");
	
	if (FullAnalysis)
	{
		if (uniType=="Binary")
		{
			orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,caseMean,
			caseStd,caseKSD,caseKSP,caseZKSD,caseZKSP,controlMean,controlStd,controlKSD,controlKSP,controlZKSD,
			controlZKSP,Beta,t.Rawvalue,t.Zvalue,wilcox.Zvalue,ZGLM,zNRI,zIDI,ROCAUC,cStatCorr,NRI,IDI,NeRI,kendall.r,
			kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p,caseN_Z_Low_Tail,caseN_Z_Hi_Tail,controlN_Z_Low_Tail,controlN_Z_Hi_Tail,Sensitivity,Specificity);
			orderframe <- orderframe[1:varinserted,];
			switch(rankingTest,
				zIDI=
				{
					orderframe <- with(orderframe,orderframe[order(-zIDI),]);
				},
				zNRI=
				{
					orderframe <- with(orderframe,orderframe[order(-zNRI),]);
				},
				IDI=
				{
					orderframe <- with(orderframe,orderframe[order(-IDI),]);
				},
				NRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NRI),]);
				},
				NeRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NeRI),]);
				},
				Ztest=
				{
					orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
				},
				AUC=
				{
					orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
				},
				Kendall=
				{
					orderframe <- with(orderframe,orderframe[order(kendall.p),]);
				},
				{
					orderframe <- with(orderframe,orderframe[order(-ROCAUC),]);
				}
			)

		}
		else
		{
			orderframe <- data.frame(Name,parent,descrip,cohortMean,cohortStd,cohortKSD,cohortKSP,cohortZKSD,cohortZKSP,Beta,ZGLM,
			NeRI,cStatCorr,spearman.r,pearson.r,kendall.r,kendall.p,BinRes.p,TstudentRes.p,WilcoxRes.p,FRes.p);
			orderframe <- orderframe[1:varinserted,];
			switch(rankingTest,
				NeRI=
				{
					orderframe <- with(orderframe,orderframe[order(-NeRI),]);
				},
				Ztest=
				{
					orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
				},
				CStat=
				{
					orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
				},
				Kendall=
				{
					orderframe <- with(orderframe,orderframe[order(kendall.p),]);
				},
				{
					orderframe <- with(orderframe,orderframe[order(-cStatCorr),]);
				}
			)
		}

	}
	else
	{
		if (is.null(acovariates)) 
		{
			acovariates <- "1";
		}
		if (uniType=="Binary")
		{
			if (categorizationType == "Raw")
			{
#				cat(timeOutcome,"\n",rankingTest,"\n",acovariates,"\n");
				univariateModels <- ForwardSelection.Model.Bin(nrow(variableList),1.0,0.0,1,acovariates,Outcome,variableList,data,1,type=type,timeOutcome=timeOutcome,selectionType="zIDI");
				orderframe <- data.frame(variableList[,1],variableList[,1],univariateModels$base.Zvalues);
				colnames(orderframe) <- c("Name","RName","ZUni");
				orderframe <- with(orderframe,orderframe[order(-ZUni),]);
			}
			else
			{
#				cat(timeOutcome,"\n",rankingTest,"\n",acovariates,"\n");
				orderframe <- data.frame(Name,parent,descrip,ZGLM,zIDI,zIDI);
				colnames(orderframe) <-  c("Name","parent","descrip","ZGLM","zIDI","ZUni");
				orderframe <- orderframe[1:varinserted,];
				switch(rankingTest,
					zIDI=
					{
						orderframe <- with(orderframe,orderframe[order(-zIDI),]);
					},
					Ztest=
					{
						orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
					}
				)
			}
		}
		else
		{
			if (categorizationType == "Raw")
			{
#				cat(timeOutcome,"\n",rankingTest,"\n",acovariates,"\n");
				univariateModels <- ForwardSelection.Model.Res(nrow(variableList),1.0,0.0,1,acovariates,Outcome,variableList,data,1,type=type,testType="Ftest",timeOutcome=timeOutcome);

				orderframe <- data.frame(variableList[,1],variableList[,1],univariateModels$base.Zvalues);
				colnames(orderframe) <- c("Name","RName","ZUni");
				orderframe <- with(orderframe,orderframe[order(-ZUni),]);
			}
			else
			{
				orderframe <- data.frame(Name,parent,descrip,ZGLM,ZGLM);
				colnames(orderframe) <-  c("Name","parent","descrip","ZGLM","ZUni");
				orderframe <- orderframe[1:varinserted,];
				orderframe <- with(orderframe,orderframe[order(-ZGLM),]);
			}
		}
#		print(parent);
	}
	row.names(orderframe) <- orderframe$Name;


	
	result <- list(orderframe=orderframe,
					variableList=variableList,
					formula=formula,
					Outcome=Outcome,
					data=data,
					categorizationType=categorizationType,
					type=type,
					rankingTest=rankingTest,
					cateGroups=cateGroups,
					raw.dataFrame=raw.dataFrame,
					description=description,
					uniType=uniType,
					acovariates=acovariates,
					timeOutcome=timeOutcome
					)
	

	return (result);
	
}
