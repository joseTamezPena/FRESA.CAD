FRESAScale <- function(data,refFrame=NULL,method=c("Norm","Order","OrderLogit","RankInv","LRankInv"),refMean=NULL,refDisp=NULL,strata=NA)
{
	if (is.null(refFrame))
	{
		refFrame <- as.data.frame(data);
	}
	
	if (inherits(refFrame,"character"))
	{
		refFrame <- refFrame[refFrame %in% rownames(data)];
		if (length(refFrame)>10)
		{
			refFrame <- as.data.frame(data[refFrame,]);
		}
		else
		{
			warning("Ref IDs not in data. All data will be used");
			refFrame <- as.data.frame(data);
		}
	}

	if (!is.null(refMean))
	{
		usedFeatures <- names(refMean);
	}
	else
	{
#		print(class(refFrame));
		usedFeatures <- colnames(refFrame);
		usedFeatures <- usedFeatures[sapply(refFrame,is.numeric)];
		if (!is.na(strata)) usedFeatures <- usedFeatures[usedFeatures != strata]
#		print(usedFeatures)
		outs <- lapply(refFrame[,usedFeatures],table);
		outl <- numeric(length(outs));
		for (i in 1:length(outs)) outl[i] <- length(outs[[i]]);
#		print(outl)

		usedFeatures <- usedFeatures[outl > 5];
#		print(usedFeatures);
	}
#	print(method);
    if (length(usedFeatures) > 0)
	{
		method <- match.arg(method);
#		srownames <- rownames(data)
		thesta <- c(1)
		if (!is.na(strata))
		{
			thesta <- unique(refFrame[,strata])
		}
		for (sta in thesta)
		{
			if (!is.na(strata))
			{
				stracondition = paste (strata,paste('==',sta));
				strastatement = paste ("subset(refFrame,",paste(stracondition,")"));
	#			cat ("Strata:",stracondition,"\n");
				cstrataref <- eval(parse(text=strastatement));
				strastatement = paste ("subset(data,",paste(stracondition,")"));
				cstrata <- eval(parse(text=strastatement));
	#			cat ("Rows:",nrow(cstrataref),"Rows 2",nrow(cstrata)," \n");
			}
			else
			{
				cstrataref = refFrame;
				cstrata = data;
			}
			datRefUses <-  as.data.frame(cstrataref[,usedFeatures]);
			colnames(datRefUses) <- usedFeatures;
			rownames(datRefUses) <- rownames(cstrataref);
			switch(method,
					Norm =
					{
						if (is.null(refMean))
						{
							refMean <- apply(datRefUses,2,mean, na.rm = TRUE);
							refDisp <- apply(datRefUses,2,sd, na.rm = TRUE);
							refDisp[refDisp == 0] <- 1.0;
						}

						meanmat <- matrix(rep(refMean,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						sdmat <- matrix(rep(refDisp,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						cstrata[,usedFeatures] <- (cstrata[,usedFeatures]-meanmat)/sdmat;
					},
					Order =
					{
						if (is.null(refMean))
						{
							refmin <- apply(datRefUses,2,min, na.rm = TRUE);
							refmax <- apply(datRefUses,2,max, na.rm = TRUE);
							refRange <- 0.5*(refmax-refmin);
							meanRange <- 0.5*(refmax+refmin);
							refRange[refRange == 0] <- 1.0;
							refMean <- apply(datRefUses,2,median, na.rm = TRUE);
							refDisp <- apply(datRefUses,2,IQR, na.rm = TRUE);
							refMean[refDisp == 0] <- meanRange[refDisp == 0];
							refDisp[refDisp == 0] <- refRange[refDisp == 0];
							refDisp <- refDisp/abs(2*qnorm(0.25));
						}

						meanmat <- matrix(rep(refMean,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						sdmat <- matrix(rep(refDisp,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						cstrata[,usedFeatures] <- (cstrata[,usedFeatures]-meanmat)/sdmat;
					},
					OrderLogit =
					{
						if (is.null(refMean))
						{
		#					cat("Here");
							refmin <- apply(datRefUses,2,min, na.rm = TRUE);
							refmax <- apply(datRefUses,2,max, na.rm = TRUE);
							refRange <- 0.5*(refmax-refmin);
							meanRange <- 0.5*(refmax+refmin);
							refRange[refRange == 0] <- 1.0;
							refMean <- 0.5*(apply(datRefUses,2,median, na.rm = TRUE) + apply(datRefUses,2,mean,trim = 0.25, na.rm = TRUE));
							refDisp <- apply(datRefUses,2,IQR, na.rm = TRUE);
							refMean[refDisp == 0] <- meanRange[refDisp == 0];
							refDisp[refDisp == 0] <- refRange[refDisp == 0];
							refDisp <- refDisp/abs(2*qnorm(0.25));
						}
						iqrsdratio = 0.5;

						meanmat <- matrix(rep(refMean,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						sdmat <- matrix(rep(refDisp,nrow(cstrata)),nrow=nrow(cstrata),ncol=length(usedFeatures),byrow = TRUE);
						cstrata[,usedFeatures] <- 4.0*(1.0/(1.0+exp(-iqrsdratio*(cstrata[,usedFeatures]-meanmat)/sdmat)) - 0.5)/iqrsdratio;
					},
					RankInv =
					{
#						data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame,strata); 			
						cstrata <- rankInverseNormalDataFrame(usedFeatures,cstrata,cstrataref); 			
					},
					LRankInv =
					{
						iqrsdratio = 0.5;
#						data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame,strata);
#						data[,usedFeatures] <- 4.0*(1.0/(1.0+exp(-iqrsdratio*data[,usedFeatures])) - 0.5)/iqrsdratio;
						cstrata <- rankInverseNormalDataFrame(usedFeatures,cstrata,cstrataref);
						cstrata[,usedFeatures] <- 4.0*(1.0/(1.0+exp(-iqrsdratio*cstrata[,usedFeatures])) - 0.5)/iqrsdratio;
					}
				)
			data[rownames(cstrata),] <- cstrata
		}
#		data <- data[srownames,]
		if (!is.null(refMean))
		{	
			names(refMean) <- usedFeatures;
			names(refDisp) <- usedFeatures;
		}
	}
	if (!is.na(strata))
	{
		refMean <- NULL
		refDisp <- NULL
	}
	scaledData=as.data.frame(data);
	attr(scaledData,"usedFeatures") <- usedFeatures;
	result <- list(scaledData=scaledData,refMean=refMean,refDisp=refDisp,strata=strata,method=method,refFrame=refFrame);
	return (result);
}
