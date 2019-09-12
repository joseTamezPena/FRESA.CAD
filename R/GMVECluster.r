#' Partitions the data based on the GMVE
#'
#' Automatically finds the best Gaussian partition of the data.
#' @param dataset Dataset from which to generate partitions.
#' @param p.threshold Threshold for clustering
#' @param p.samplingthreshold p value for cluster candidates estimation
#' @param samples Number of random samples
#' @param sampling.rate the uniform sampling rate
#' @param jitter if true, add noise to the internal copy data set
#' @param tryouts the number of random samples per candidate point
#' @return Returns a list containing:
#' \describe{
#'   \item{\code{cluster}}{The cluster label of each data set}
#'   \item{\code{clusterMean}}{A list with the cluster centroid}
#'   \item{\code{clusterCov}}{A list with the cluster covariance}
#'   \item{\code{robCov}}{The list of the robust covariances}
#'   \item{\code{pvals}}{The vector of cluster p-values}
#'   \item{\code{k}}{The number of clusters}
#' }
#' @examples
#' data  <- load("my_data.RData")
#' res   <- GMVECluster(dataset = data)
#' labels <- res$cluster
#' @importFrom 
#' @export

GMVECluster <- function(dataset, p.threshold=0.975,samples=10000,p.samplingthreshold=0.50,sampling.rate = 3,jitter=TRUE,tryouts=25,verbose=FALSE)
{

  if (!requireNamespace("robustbase", quietly = TRUE)) {
	  install.packages("robustbase", dependencies = TRUE)
	  }
	  
	a = as.numeric(Sys.time());
	set.seed(a);

	intdata <- dataset
	features <- colnames(dataset);
	ndata <- nrow(intdata);
	ClusterLabels <- numeric(ndata);
	p <- ncol(intdata);
	if (ndata > samples)
	{
		intdata <- dataset[sample(ndata,samples),]
	}
	ndata <- nrow(intdata);
	chithreshold <- qchisq(p.threshold,p);
	chithreshold2 <- qchisq(p.threshold/2,p);
	chithreshold3 <- qchisq(p.threshold/4,p);
	chithreshold_out <- qchisq(0.999,p);
	samplingthreshold <- qchisq(p.samplingthreshold,p);
#	print(samplingthreshold);
	maxp <- 1.0;
	minD <- 1.0;
	k <- 0;
	bestCov <- list();
	bestmean <- list();
	robCov <- list();
	pvals <- numeric();
	## the list of possible volume fractions
	alphalist <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50);
	hlist <- as.integer(ndata*alphalist);
	p1 <- p + 1;
	minpvalThr <- 1.0e-2;
	globalcov <- cov(intdata);
	vvar <- diag(globalcov);
	minvvar <- vvar/(ndata*ndata);
	gmincov <- diag(minvvar);
	minminVar <- det(gmincov);
	detmincov <- cov(gmincov);
	jitteredData <- NULL;
	if (jitter)
	{
		minunit <- 2.0*apply(intdata,2,IQR)/ndata;
		gmincov <- diag(pmin(minvvar,minunit*minunit));
		for (i in 1:p)
		{
			if (class(intdata[,i]) == "integer")
			{
				intdata[,i] <- as.numeric(intdata[,i]+max(1.0,minunit[i])*rnorm(ndata));
			} 
			else 
			{
				intdata[,i] <- intdata[,i]+max(sqrt(minvvar[i]),minunit[i])*rnorm(ndata);
			}
		}
		jitteredData <- intdata;
	}
		
	JClusters <- 1;
	cycles <- 0;
	maxMahadis <- numeric(ndata);
	h0 <- max(15,p+2); #at least 15 samples in a cluster
	andata <- ndata;
	## Loop unit all clusters are found
	while ((andata >= h0) && (cycles < 5))
	{
		k <- k + 1;
		maxp <- 0.0;
		minD <- 1.0;
		colmean <- list();
		covmat <- list();
		mdistlist <- list();
		detcovmat <- numeric();
		JClusters <- 0;
		auxdata <- intdata[maxMahadis == 0,];
		andata <- nrow(auxdata);
		auxdata <- auxdata[sample(andata),];
	## Loop for cluster candidates
		cat(":");
		for (i in sampling.rate*(1:as.integer(andata/sampling.rate)))
		{
			datao <- as.numeric(auxdata[i,]);
			smdist <- mahalanobis(intdata,datao,globalcov);
			qdata <- intdata[(smdist > 0) & (smdist < samplingthreshold),];
#			print(nrow(qdata))
			if (length(qdata) > 0)
			{
				if (nrow(qdata) >= h0)
				{
					for (j in 1:tryouts)
					{
						sdata <- rbind(datao,qdata[sample(nrow(qdata),p),]);
						jmean <- apply(sdata,2,mean);
						jcov <- cov(sdata)+gmincov;
						jcovDet <- try(det(jcov));
						if ( !inherits(jcovDet, "try-error") && !is.nan(jcovDet) && !is.na(jcovDet) )
						{
							mdist <- try(mahalanobis(intdata,jmean,jcov));
							if ((!inherits(mdist, "try-error")) && (jcovDet > minminVar))
							{
								JClusters <- JClusters + 1;
								mdistlist[[JClusters]] <- mdist[order(mdist)];
								colmean[[JClusters]] <- jmean;
								covmat[[JClusters]] <- jcov;
								detcovmat[JClusters] <- jcovDet^(1.0/(2.0*p));
							}
						}
					}
				}
			}
		}
		cat("-");
		atalpha <- 0;
		ptsinside <- 0;
		inside <- numeric(0);
		ondata <- ndata;
		inside.centroid <- 0;
		maxpD <- -1;
		refinecount <- 0;
		if (JClusters > 0)
		{
	## Loop for different cluster volumes fractions and select the one with the closer Gaussian distribution and minimum volume
			bestCMean <- numeric(p);
			bestCCov <- diag(bestCMean);
			Ellipsoidvol <- numeric(JClusters);
			cat(":");
			for ( alpha in alphalist )
			{
				h <- as.integer(alpha*ndata+0.5);
				if (h >= p1) 
				{
					for (i in 1:JClusters)
					{
						Ellipsoidvol[i] <- mdistlist[[i]][h]*detcovmat[i];
					}
					minEllipVol <- which.min(Ellipsoidvol)[1];
					mincentroid <- colmean[[minEllipVol]];
					minCov <- covmat[[minEllipVol]];
					mdist <- mdistlist[[minEllipVol]];
					correction <- mdist[h]/qchisq(alpha,p);
					minCov <- minCov*correction;
					mdist <- mahalanobis(intdata,mincentroid,minCov);
					inside <- (mdist < chithreshold);
					if (!is.na(sum(inside)))
					{
						ptsinside <- sum(inside)
						if (ptsinside >= h0)
						{
							newCentroid <- apply(intdata[inside,],2,mean);
							newCovariance <- cov(intdata[inside,])+gmincov;
							distanceupdate <- mahalanobis(intdata[inside,],newCentroid,newCovariance);
							distanceupdate <- distanceupdate[order(distanceupdate)];
							mvsd <- as.integer(ptsinside/p.threshold+0.5);
							dsample <- (1:mvsd)/mvsd;
							if (dsample[ptsinside] > p.threshold)
							{
								dsample[ptsinside] <- p.threshold;
							}
							disTheoretical <- qchisq(dsample[1:ptsinside],p);
							kst <- ks.test(disTheoretical,distanceupdate + rnorm(length(distanceupdate),0,1e-10));
							if (kst$p.value > maxp)
							{
								 bestmean[[k]] <- newCentroid;
								 bestCov[[k]] <- newCovariance;
								 bestCMean <- newCentroid;
								 bestCCov <- newCovariance;
								 maxp <- kst$p.value;
								 minD <- kst$statistic;
								 pvals[k] <- maxp;
								 atalpha <- alpha;
							}
						}
					}
				}
			}
			cat("-");
## Now use the similar centroids to refine the candidate clusters
			refinecount <- 1;
			cat(":");
			if (maxp > minpvalThr)
			{
				for ( alpha in alphalist )
				{
					h <- as.integer(alpha*ndata+0.5);
					if (h >= p1) 
					{
						Ellipsoidvol <- numeric(JClusters);
						for (i in 1:JClusters)
						{
							Ellipsoidvol[i] <- mdistlist[[i]][h]*detcovmat[i];
						}
						minEllipVol <- which.min(Ellipsoidvol)[1];
						mincentroid <- colmean[[minEllipVol]];
						minCov <- covmat[[minEllipVol]];
						mdist <- mdistlist[[minEllipVol]];
						correction <- mdist[h]/qchisq(alpha,p);
						minCov <- minCov*correction;
						mdist <- mahalanobis(intdata,mincentroid,minCov);
						inside <- (mdist < chithreshold);
						if (!is.na(sum(inside)))
						{
							ptsinside <- sum(inside)
							if (ptsinside >= h0)
							{
								newCentroid <- apply(intdata[inside,],2,mean);
								newCovariance <- cov(intdata[inside,]);
								distanceupdate <- mahalanobis(intdata[inside,],newCentroid,newCovariance);
								distancecluster1 <- mahalanobis(bestCMean,newCentroid,newCovariance);
								distancecluster2 <- mahalanobis(newCentroid,bestCMean,bestCCov);
								distanceupdate <- distanceupdate[order(distanceupdate)]
								mvsd <- as.integer(ptsinside/p.threshold+0.5);
								dsample <- (1:mvsd)/mvsd;
								if (dsample[ptsinside] > p.threshold)
								{
									dsample[ptsinside] <- p.threshold;
								}
								disTheoretical <- qchisq(dsample[1:ptsinside],p);
								kst <- ks.test(disTheoretical,distanceupdate + rnorm(length(distanceupdate),0,1e-10));
								if ( ((kst$p.value >= minpvalThr) || (kst$statistic <= minD)) && (distancecluster1 < chithreshold2) && (distancecluster2 < chithreshold2) )
								{
									refinecount <- refinecount + kst$p.value;
									bestmean[[k]] <- bestmean[[k]] + kst$p.value*newCentroid;
									bestCov[[k]] <-bestCov[[k]] + kst$p.value*newCovariance;
									pvals[k] <- pvals[k]+ kst$p.value*maxp;
									atalpha <- atalpha + kst$p.value*alpha;
								}
							}
						}
					}
				}
			}
			cat("-");
			bestmean[[k]] <- bestmean[[k]]/refinecount;
			bestCov[[k]] <- (bestCov[[k]]/refinecount)/p.threshold;
			pvals[k] <- pvals[k]/refinecount;
			atalpha <- atalpha/refinecount;
			inside.centroid <- 0;
			## Check for cluster overlap
			if ((k > 1) && (maxp >= minpvalThr))
			{
				for (i in 1:(k-1))
				{
					inside.centroid <- inside.centroid + 0.25*(mahalanobis(bestmean[[i]],bestmean[[k]],bestCov[[k]]) < chithreshold3) + 0.75*(mahalanobis(bestmean[[k]],bestmean[[i]],bestCov[[i]]) < chithreshold3);
				}
			}
			inside <- numeric(0);
			## Include cluster only if p value is significant and no overlap with already discovered clusters
			if ((maxp < minpvalThr) || (inside.centroid > 0.5))
			{
				cycles <- cycles + 1;
				bestmean[[k]] <- NULL;
				bestCov[[k]] <- NULL;
				robCov[[k]] <- NULL;
				k <- k - 1;
				p.threshold <- 0.9*p.threshold;
				minpvalThr <- minpvalThr/2.0;
			}
			else
			{
				mdist <- mahalanobis(intdata,bestmean[[k]],bestCov[[k]]);
				inside <- (mdist < chithreshold);
				cludata <- intdata[inside,];
				if (nrow(cludata) > p1)
				{
					bestCov[[k]] <- cov(cludata);
					bestmean[[k]] <- apply(cludata,2,mean);
					lcov <- try(robustbase::covMcd(cludata));
					if (!inherits(lcov, "try-error"))
					{
						robCov[[k]] <- list(centroid=lcov$center,cov=lcov$cov);
					}
					else
					{
						robCov[[k]] <- list(centroid=bestmean[[k]],cov=bestCov[[k]]);
					}
					intdata <- intdata[!inside,];
					ndata <- nrow(intdata);
					maxMahadis <- maxMahadis[!inside] + 1*(mahalanobis(intdata,bestmean[[k]] ,bestCov[[k]]) < chithreshold_out);
				}
				else
				{
					cycles <- cycles + 1;
					bestmean[[k]] <- NULL;
					bestCov[[k]] <- NULL;
					robCov[[k]] <- NULL;
					k <- k - 1;
					p.threshold <- 0.9*p.threshold;
					minpvalThr <- minpvalThr/2.0;
				}
			}
		}
		else
		{
			cycles <- cycles + 1;
			k <- k - 1;
			p.threshold <- 0.9*p.threshold;
			minpvalThr <- minpvalThr/2.0;
		}
		if (verbose) 
		{
			cat(cycles,":",k,":",inside.centroid,":(",ondata,"->",andata,"):",atalpha,":",maxp,":",minD,":",JClusters,":",refinecount,":",sum(inside),"\n");
		}
		else 
		{
			cat("~");
		}
		chithreshold <- qchisq(p.threshold,p);
		chithreshold2 <- qchisq(p.threshold/2,p);
	}
	k <- length(bestmean);
	## assign clusters labels to all points
	if (k > 0)
	{
		ClusterLabels <- nearestCentroid(dataset,bestmean,bestCov,0);
	}
	cat("(",k,")\n");
	
	result <- list(
		cluster = ClusterLabels,
		classification = ClusterLabels,
		centers = bestmean,
		covariances = bestCov,
		robCov = robCov,
		pvals = pvals,
		k = k,
		features = features,
		jitteredData = jitteredData
	  )
	 class(result) <- "GMVE"
	 return (result);
}

predict.GMVE <- function(object,...)
{
	parameters <- list(...);
	testData <- parameters[[1]];
	thr <- 0;
	if (length(parameters) > 1)
	{
		thr <- parameters[[2]];
	}
	pLS <- list(classification=nearestCentroid(testData[,object$features],object$centers,object$covariances,thr));
	return (pLS);
}
