#' labels data sets based on nearest centroid
#'
#' Based on the mean and covariance matrixes it will lables the dataset points.
#' @param dataset Dataset from which to label
#' @param clustermean the mean of clusters
#' @param clustercov the covariance of the cluster
#' @param p.threshold Threshold for outliers detection
#' @return Returns the labels of each dataset point:
#' @examples
#' data	<- load("my_data.RData")
#' labels <-	 <- nearestCentroid(dataset = data)
#' @importFrom 
#' @export

nearestCentroid <- function(dataset,clustermean=NULL,clustercov=NULL, p.threshold=1.0e-6)
{
	if (class(clustermean) == "matrix")
	{
		centers <- clustermean;
		clustermean <- list();
		for (j in 1:nrow(centers))
		{
			clustermean[[j]] <- as.numeric(centers[j,]);
		}	
	}
	distance <- numeric(nrow(dataset));
	ClusterLabels <- numeric(nrow(dataset));
	if (is.null(clustermean) && (class(clustercov) == "list"))
	{
		if (length(clustercov)>0)
		{
			if (!is.null(clustercov[[1]]$centroid))
			{
				clustermean <-	list();
				robCov <- clustercov;
				clustercov <- list();
				for (j in 1:length(robCov))
				{
					clustermean[[j]] <- robCov[[j]]$centroid
					clustercov[[j]] <- robCov[[j]]$cov
				}
			}
		}
	}
	if (!is.null(clustermean))
	{
		k <- length(clustermean);
		## assign clusters labels to all points
		if (k > 0)
		{
			if (p.threshold>0)
			{
				chithreshold_out <- qchisq(1.0-p.threshold,length(clustermean[[1]]));
			}
			else
			{
				chithreshold_out <- 1.0e100;
			}
			if (is.null(clustercov))
			{
				clustercov <- list();
				ones <- rep(1,length(clustermean[[1]]));
				for (i in 1:k)
				{
					clustercov[[i]] <- diag(ones);
				}
				chithreshold_out <- 1.0e100;
			}
			distancemaha <- matrix(0,k,nrow(dataset));
			for (i in 1:k)
			{
				distancemaha[i,] <- mahalanobis(dataset,clustermean[[i]],clustercov[[i]]);
			}
			for (i in 1:nrow(dataset))
			{
				ClusterLabels[i] <- which.min(distancemaha[,i])[1];
				distance[i] <- distancemaha[ClusterLabels[i],i];
				if (distancemaha[ClusterLabels[i],i] > chithreshold_out)
				{
					ClusterLabels[i] <- 0;
				}
			}
		}
	}
	attr(ClusterLabels,"distance") <- distance;
	return (ClusterLabels)
}