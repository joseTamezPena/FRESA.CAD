#' returns the jaccard index matrix of two labeled sets
#'
#' jaccard = intersection/union
#' @param clustersA cluster labels of set 
#' @param clustersB second cluster labels of the same set
#' @return Returns the jaccard matrix of the labeled set
#' @examples
#' @importFrom 
#' @export

clusterISODATA <- function(dataset,clusteringMethod=GMVECluster,trainFraction=0.99,randomTests=10,jaccardThreshold=0.45,isoDataThreshold=0.75,plot=TRUE,...)
{

	classLabels <- list();
	clustercenters <- list();
	ClusterLabels <- numeric();
	isoDataCov <- list();
	sumDataCov <- list();
	clusterPvals <- list();


	for (i in 1:randomTests) 
	{ 
		cat("[");
		res <- try(clusteringMethod(dataset[sample(nrow(dataset),trainFraction*nrow(dataset)),],...));
		cat("]");
		classLabels[[i]] <- nearestCentroid(dataset,clustercov = res$robCov);
		clustercenters[[i]] <- res$robCov;
		clusterPvals[[i]] <- res$pvals;	
		if (plot) plot(dataset[,c(1:2)],col = classLabels[[i]],main=as.character(i));
	}
	
	n = 0;
	perTestJaccard <- list();
	pointjaccard <- numeric(nrow(dataset));
	for (i in 1:randomTests)
	{
		attestjaccard <- numeric(nrow(dataset));
		for (j in 1:randomTests)
		{ 
			if (i != j)
			{
				jaccard <- jaccardMatrix(classLabels[[i]],classLabels[[j]]);
				attestjaccard <- attestjaccard + jaccard$elementJaccard;
				n <- n + 1;
			}
		}
		pointjaccard <- pointjaccard + attestjaccard;
		perTestJaccard[[i]] <-  attestjaccard/(randomTests-1);
	}
	pointjaccard <- pointjaccard/n;
	if (plot) plot(dataset[,c(1:2)],col = 1*(pointjaccard >= jaccardThreshold) );
	inpoints <- (pointjaccard > jaccardThreshold);
	goodpoints <- sum(inpoints);
	isoclusters <- 0;
	if (goodpoints > 0)
	{
		p <- ncol(dataset)
		zeroCov <- list(centroid=rep(0,p),cov=diag(rep(0,p)));
		goodjaccard <- pointjaccard[inpoints];
		orderJaccard <- order(-goodjaccard);
		goodjaccard <- goodjaccard[orderJaccard];
#		print(c(mean(pointjaccard),goodpoints));
		for (i in 1:randomTests)
		{	
			classLabels[[i]] <- classLabels[[i]][inpoints];
			classLabels[[i]] <- classLabels[[i]][orderJaccard];
			perTestJaccard[[i]] <- perTestJaccard[[i]][inpoints];
			perTestJaccard[[i]] <- perTestJaccard[[i]][orderJaccard];
		}

		inclass <- 0;
		pcov <- zeroCov;
		for (j in 1:randomTests)
		{
			label <- classLabels[[j]][1];
			if (label > 0)
			{
				wts <- (perTestJaccard[[j]][1]-jaccardThreshold)*(perTestJaccard[[j]][1] > jaccardThreshold);
				wts <- wts*wts*clusterPvals[[j]][label];
				pcov$centroid <-  pcov$centroid + wts*clustercenters[[j]][[label]]$centroid;
				pcov$cov <- pcov$cov + wts*clustercenters[[j]][[label]]$cov;
				inclass <- inclass + wts;
			}
		}
		if (inclass > 0)
		{
			isoclusters <- 1;
			totwts <- numeric();
			sumDataCov[[isoclusters]] <- pcov;
			totwts <-  append(totwts,inclass);
			pcov$centroid <-  pcov$centroid/inclass;
			pcov$cov <- pcov$cov/inclass;
			isoDataCov[[isoclusters]] <- pcov;
			if (goodpoints > 1)
			{
				disThreshold <- qchisq(isoDataThreshold,p);
				for (i in 2:goodpoints)
				{
					inclass <- 0;
					pcov <- zeroCov;
					for (j in 1:randomTests)
					{
						label <- classLabels[[j]][i];
						if (label > 0)
						{
							wts <- (perTestJaccard[[j]][i]-jaccardThreshold)*(perTestJaccard[[j]][i] > jaccardThreshold);
							wts <- wts*wts*clusterPvals[[j]][label];
							pcov$centroid <-  pcov$centroid + wts*clustercenters[[j]][[label]]$centroid;
							pcov$cov <- pcov$cov + wts*clustercenters[[j]][[label]]$cov;
							inclass <- inclass + wts;
						}
					}
					if (inclass > 0)
					{
						rcov <- pcov;
						pcov$centroid <-  pcov$centroid/inclass;
						pcov$cov <- pcov$cov/inclass;
						clusteroverlap <- FALSE;
						j <- 0;
						while (!clusteroverlap && (j < isoclusters))
						{
							j <- j + 1;
							distancecluster1 <- mahalanobis(pcov$centroid,isoDataCov[[j]]$centroid,isoDataCov[[j]]$cov);
							distancecluster2 <- mahalanobis(isoDataCov[[j]]$centroid,pcov$centroid,pcov$cov);
							clusteroverlap <- (distancecluster1 < disThreshold) || (distancecluster2 < disThreshold);
						}
						if (!clusteroverlap)
						{
							isoclusters <- isoclusters + 1;
							isoDataCov[[isoclusters]] <- pcov;
							sumDataCov[[isoclusters]] <- rcov;
							totwts <-  append(totwts,inclass);
							cat("(",i,",",isoclusters,"):")
						}
						else
						{
							sumDataCov[[j]]$centroid <- sumDataCov[[j]]$centroid + inclass*pcov$centroid;
							sumDataCov[[j]]$cov <- sumDataCov[[j]]$cov + inclass*pcov$cov;
							totwts[j] <-  totwts[j] + inclass;
						}
					}
				}
			}
		}
	}
	if (isoclusters > 0)
	{
		for (i in 1:isoclusters)
		{
			isoDataCov[[i]]$centroid <- sumDataCov[[i]]$centroid/totwts[i];
			isoDataCov[[i]]$cov <- sumDataCov[[i]]$cov/totwts[i];
		}
		ClusterLabels <- nearestCentroid(dataset,clustercov=isoDataCov);
	}



	result <- list(
		cluster = ClusterLabels,
		classification = ClusterLabels,
		robustCovariance = isoDataCov,
		pointjaccard = pointjaccard
		)
	
	return (result)
}