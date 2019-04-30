nearestNeighborImpute <- function(tobeimputed,referenceSet=NULL,distol=1.05)
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
	trainset <- trainset[complete.cases(trainset),]
	imputeddata <- tobeimputed;
	medianvalues <-  as.numeric(apply(trainset,2,median, na.rm = TRUE));
	IQRvalues <-  as.numeric(apply(trainset,2,IQR, na.rm = TRUE));
	IQRvalues[IQRvalues==0] <- 1;
	for (i in 1:nrow(imputeddata))
	{
		nacol <- is.na(imputeddata[i,]);
		if ((i %% 10)==0) cat(".");
		if ((i %% 500)==0) cat(i,"\n");
		if (any(nacol))
		{	
			if (sum(1*(!nacol)) == 0)
			{
				imputeddata[i,] <- medianvalues;
			}
			else
			{
				redtrain <- trainset[,!nacol];
				redimputed <- as.numeric(imputeddata[i,!nacol]);
				distance <- abs(sweep(redtrain,2,redimputed,"-"))
				distance <- sweep(distance,2,IQRvalues[!nacol],"/");
				if (sum(!nacol) > 1)
				{
					distance <- apply(distance,1,mean, na.rm = TRUE);
				}
				distance <- as.numeric(distance);
				mindistance <- distol*min(distance);
				wsmaller <- (distance<=mindistance);
				utrainset <- trainset[wsmaller,nacol];
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
	return (imputeddata)
}
