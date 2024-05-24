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
		datRefUses <-  as.data.frame(refFrame[,usedFeatures]);
		colnames(datRefUses) <- usedFeatures;
		rownames(datRefUses) <- rownames(refFrame);
		srownames <- rownames(data)

		switch(method,
				Norm =
				{
					if (is.null(refMean))
					{
						refMean <- apply(datRefUses,2,mean, na.rm = TRUE);
						refDisp <- apply(datRefUses,2,sd, na.rm = TRUE);
						refDisp[refDisp == 0] <- 1.0;
					}

					meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					data[,usedFeatures] <- (data[,usedFeatures]-meanmat)/sdmat;
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

					meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					data[,usedFeatures] <- (data[,usedFeatures]-meanmat)/sdmat;
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

					meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
					data[,usedFeatures] <- 4.0*(1.0/(1.0+exp(-iqrsdratio*(data[,usedFeatures]-meanmat)/sdmat)) - 0.5)/iqrsdratio;
				},
				RankInv =
				{
					data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame,strata); 			
				},
				LRankInv =
				{
					iqrsdratio = 0.5;
					data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame,strata);
					data[,usedFeatures] <- 4.0*(1.0/(1.0+exp(-iqrsdratio*data[,usedFeatures])) - 0.5)/iqrsdratio;
				}
			)
		data <- data[srownames,]
		if (!is.null(refMean))
		{	
			names(refMean) <- usedFeatures;
			names(refDisp) <- usedFeatures;
		}
	}
	scaledData=as.data.frame(data);
	attr(scaledData,"usedFeatures") <- usedFeatures;
	result <- list(scaledData=scaledData,refMean=refMean,refDisp=refDisp,strata=strata,method=method,refFrame=refFrame);
	return (result);
}
