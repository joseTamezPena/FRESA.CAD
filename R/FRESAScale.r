FRESAScale <- function(data,refFrame=NULL,method=c("Norm","Order","OrderLogit","RankInv"),refMean=NULL,refDisp=NULL,strata=NA)
{
	if (is.null(refFrame))
	{
		refFrame <- as.data.frame(data);
	}

	if (!is.null(refMean))
	{
		usedFeatures <- names(refMean);
	}
	else
	{
#		print(class(refFrame));
		usedFeatures <- colnames(refFrame);
#		print(usedFeatures)
		outs <- lapply(refFrame,table);
		outl <- numeric(length(outs));
		for (i in 1:length(outs)) outl[i] <- length(outs[[i]]);
#		print(outl)
		usedFeatures <- usedFeatures[outl > 5];
#		print(usedFeatures);
	}
#	print(method);
	method <- match.arg(method);
	datRefUses <-  as.data.frame(refFrame[,usedFeatures]);
	colnames(datRefUses) <- usedFeatures;
	rownames(datRefUses) <- rownames(refFrame);

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
					refMean <- apply(datRefUses,2,median, na.rm = TRUE);
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
			}
		)
	if (!is.null(refMean))
	{	
		names(refMean) <- usedFeatures;
		names(refDisp) <- usedFeatures;
	}

	result <- list(scaledData=data,refMean=refMean,refDisp=refDisp,strata=strata,method=method);
	return (result);
}
