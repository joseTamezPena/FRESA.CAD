FRESAScale <- function(data,refFrame=NULL,method=c("Norm","Order","OrderLogit","RankInv"),refMean=NULL,refDisp=NULL,strata=NA)
{
	if (is.null(refFrame))
	{
		refFrame <- data;
	}

	usedFeatures <- colnames(refFrame);

	outs <- sapply(refFrame,table);
	outl <- numeric(length(outs));
	for (i in 1:length(outs)) outl[i] <- length(outs[[i]]);
	usedFeatures <- usedFeatures[outl > 10];
	method <- match.arg(method);
	switch(method,
			Norm =
			{
				if (is.null(refMean))
				{
					refMean <- apply(refFrame[,usedFeatures],2,mean, na.rm = TRUE);
					refDisp <- apply(refFrame[,usedFeatures],2,sd, na.rm = TRUE);
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
					refSD <- apply(refFrame[,usedFeatures],2,sd, na.rm = TRUE);
					refSD[refSD == 0] <- 1.0;
					refMean <- apply(refFrame[,usedFeatures],2,median, na.rm = TRUE);
					refDisp <- apply(refFrame[,usedFeatures],2,IQR, na.rm = TRUE);
					refDisp[refDisp == 0] <- refSD[refDisp == 0];
				}

				meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
				sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
				data[,usedFeatures] <- (data[,usedFeatures]-meanmat)/sdmat;
			},
			OrderLogit =
			{
				if (is.null(refMean))
				{
					refSD <- apply(refFrame[,usedFeatures],2,sd, na.rm = TRUE);
					refSD[refSD == 0] <- 1.0;
					refMean <- apply(refFrame[,usedFeatures],2,median, na.rm = TRUE);
					refDisp <- apply(refFrame[,usedFeatures],2,IQR, na.rm = TRUE);
					refDisp[refDisp == 0] <- refSD[refDisp == 0];
				}
				iqrsdratio = 2.0;

				meanmat <- matrix(rep(refMean,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
				sdmat <- matrix(rep(refDisp,nrow(data)),nrow=nrow(data),ncol=length(usedFeatures),byrow = TRUE);
				data[,usedFeatures] <- 1.0/(1.0+exp(-iqrsdratio*(data[,usedFeatures]-meanmat)/sdmat));
			},
			RankInv =
			{
				data <- rankInverseNormalDataFrame(usedFeatures,data,refFrame,strata); 			
			}
		)
	result <- list(scaledData=data,refMean=refMean,refDisp=refDisp,strata=strata);
	return (result);
}
