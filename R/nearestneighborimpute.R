nearestNeighborImpute <- function(tobeimputed,referenceSet=NULL,catgoricCol=NULL,distol=1.05,useorder=TRUE)
{	
	if (is.null(referenceSet))
	{
		trainset <- tobeimputed;
	}
	else
	{	
		rowsnotin <- !(rownames(referenceSet) %in% rownames(tobeimputed))
		trainset <- rbind(referenceSet[rowsnotin,colnames(tobeimputed)],tobeimputed);
	}
	imputeddata <- tobeimputed;
	medianvalues <-  as.numeric(apply(trainset,2,median, na.rm = TRUE));
	IQRvalues <-  as.numeric(apply(trainset,2,IQR, na.rm = TRUE));
	sdvalues <-  as.numeric(apply(trainset,2,sd, na.rm = TRUE));
#	trainset <- trainset[complete.cases(trainset),]
	IQRvalues[IQRvalues==0] <- sdvalues[IQRvalues==0];
	IQRvalues[IQRvalues==0] <- 1;
	catvalues <- NULL
	if (!is.null(catgoricCol))
	{
		catvalues <- imputeddata[1,catgoricCol]
	}
	else
	{
		useorder <- FALSE
	}
#	print(IQRvalues)
	
	for (i in 1:nrow(imputeddata))
	{
#		print(rownames(imputeddata)[i]);
		nacol <- is.na(imputeddata[i,]);
		if ((i %% 10)==0) cat(".");
		if ((i %% 500)==0) cat(i,"\n");
		if (any(nacol))
		{	
			if (useorder)
			{
				if (sum(abs(imputeddata[i,catgoricCol] - catvalues)) == 0)
				{
					if (i>1)
					{
						imputeddata[i,nacol] <- imputeddata[i-1,nacol];
					}
				}
				catvalues <- imputeddata[i,catgoricCol]
				nacol <- is.na(imputeddata[i,]);
			}
			if (any(nacol))
			{
				if (sum(1*(!nacol)) == 0)
				{
					imputeddata[i,] <- medianvalues;
				}
				else
				{
					redtrain <- trainset[,!nacol];
					datatrain <- as.data.frame(trainset[,nacol]);
					theCompleteCases <- complete.cases(datatrain)
					datatrain <- as.data.frame(datatrain[theCompleteCases,])
					redtrain <- redtrain[theCompleteCases,]
					redimputed <- as.numeric(imputeddata[i,!nacol]);
					distance <- abs(sweep(redtrain,2,redimputed,"-"))
					distance <- sweep(distance,2,IQRvalues[!nacol],"/");
					distance[distance > 1.0] <- 1.0
					if (!is.null(catgoricCol))
					{
						colnames(distance) <- colnames(redtrain)
						distance[,catgoricCol] <- 100.0*(distance[,catgoricCol] > 0)
					}
					
					if (sum(!nacol) > 1)
					{
						distance <- apply(distance,1,mean, na.rm = TRUE);
					}
					distance <- as.numeric(distance);
					mindistance <- distol*min(distance);
					wsmaller <- (distance<=mindistance);
					utrainset <- datatrain[wsmaller,];

					if (sum(wsmaller)==1)
					{
						imputeddata[i,nacol] <- as.numeric(utrainset);
					}
					else
					{
						if (sum(nacol) > 1)
						{
							mv <- as.numeric(apply(utrainset,2,median, na.rm = TRUE));
						}
						else
						{
							mv <- as.numeric(median(utrainset,na.rm = TRUE));
						}
						imputeddata[i,nacol] <- mv;
					}
				}
			}
		}
	}
	return (imputeddata)
}
