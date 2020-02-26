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

GMVECluster <- function(dataset, p.threshold=0.975,samples=10000,p.samplingthreshold = 0.50,sampling.rate = 3,jitter=TRUE,tryouts=25,verbose=FALSE)
{

  if (!requireNamespace("robustbase", quietly = TRUE)) {
	  install.packages("robustbase", dependencies = TRUE)
	  }
	  
	intdata <- dataset
	features <- colnames(dataset);
	ndata <- nrow(intdata);
	ClusterLabels <- numeric(ndata);
	p <- ncol(intdata);
	proot <- 1.0/(2.0*p);
	if (ndata > samples)
	{
		intdata <- dataset[sample(ndata,samples),]
	}
	ndata <- nrow(intdata);
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
	alphalist <- c(0.050,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,0.950);
	alphalist2 <- c(0.50,0.60,0.70,0.80,0.90);
	hlist <- as.integer(alphalist*ndata+0.5);
	p1 <- p + 1;
	minpvalThr <- 0.001;
	globalcov <- cov(intdata);
	vvar <- diag(globalcov);
	globalcov <- diag(vvar);
	minvvar <- vvar/(ndata*ndata);
	gmincov <- diag(minvvar);
	minminVar <- det(gmincov);
	JClusters <- 1;
	cycles <- 0;
	maxMahadis <- numeric(ndata);
	h0 <- max(15,p+2); #at least 15 samples in a cluster
	thedataCpy <- intdata;
	andata <- ndata;
	## Loop unit all clusters are found
	cat("{");
	while ((andata >= h0) && (cycles < 5))
	{
		
		chithreshold <- qchisq(p.threshold,p);
		chithreshold3 <- qchisq(0.75*p.threshold,p);
		chithreshold_out <- qchisq(0.5*(0.999+p.threshold),p);

		k <- k + 1;
		maxp <- minpvalThr;
		minD <- 1.0;
		colmean <- list();
		covmat <- list();
		mdistlist <- list();
		detcovmat <- numeric();
		JClusters <- 0;
		bestmean[[k]] <- apply(intdata,2,mean);
		bestCov[[k]] <- diag(diag(cov(intdata)));
		robCov[[k]] <- list(centroid=bestmean[[k]],cov=bestCov[[k]]);

		if (jitter)
		{
			minunit <- 2.0*apply(intdata,2,IQR)/ndata;
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
			if (k == 1)
			{
				thedataCpy <- intdata;
			}
		}

		if (sum(1*(maxMahadis == 0)) < h0)
		{
			auxdata <- intdata;
			andata <- 0;
		}
		else
		{
			auxdata <- intdata[maxMahadis == 0,];
			andata <- nrow(auxdata);
		}
		if (length(andata) == 0 ) 
		{
			andata <- 0;
		}

		if (andata >= h0)
		{
	#		print(andata);
			auxdata <- auxdata[sample(andata),];
		## Loop for cluster candidates
			cat(":");
			for (i in sampling.rate*(1:as.integer(andata/sampling.rate)))
			{
				datao <- as.numeric(auxdata[i,]);
				names(datao) <- colnames(auxdata);
#				jmean <- matrix(rep(datao,p1),ncol=p,nrow=p1,byrow=TRUE);
				smdist <- mahalanobis(intdata,datao,globalcov);
				qdata <- intdata[(smdist > 0) & (smdist < samplingthreshold),];
	#			print(nrow(qdata))
				if (!is.null(qdata))
				{
					if (nrow(qdata) >= h0)
					{
						for (j in 1:tryouts)
						{
							sdata <-  qdata[sample(nrow(qdata),p1),];
							jmean <- apply(sdata,2,mean);
							jcov <- cov(sdata);
							jcovDet <- try(det(jcov));
							if ( !inherits(jcovDet, "try-error") && !is.nan(jcovDet) && !is.na(jcovDet) )
							{
								if (jcovDet > minminVar)
								{
									mdist <- try(mahalanobis(intdata,jmean,jcov));
									if (!inherits(mdist, "try-error"))
									{
										JClusters <- JClusters + 1;
										mdistlist[[JClusters]] <- mdist[order(mdist)];
										colmean[[JClusters]] <- jmean;
										covmat[[JClusters]] <- jcov;
										detcovmat[JClusters] <- jcovDet^proot;
									}
								}
								else
								{
									cat("|",jcovDet,"|")
								}
							}
						}						
					}
				}
			}
		}
		else
		{
			andata <- 0;
		}
		cat("[");
		atalpha <- 0;
		ptsinside <- 0;
		inside <- numeric(0);
		ondata <- ndata;
		inside.centroid <- 0;
		maxpD <- -1;
		optsinside <- 0;
		if (JClusters > 0)
		{
	## Loop for different cluster volumes fractions and select the one with the closer Gaussian distribution and minimum volume
			Ellipsoidvol <- numeric(JClusters);
			cat("\\");
			bcorrection <- 1;
			for (h in hlist)
			{
				if ((h >= p1) && (h < andata))
				{
					alpha <- h/andata;
					for (i in 1:JClusters)
					{
						Ellipsoidvol[i] <- mdistlist[[i]][h]*detcovmat[i];
					}
					minEllipVol <- which.min(Ellipsoidvol)[1];
					mincentroid <- colmean[[minEllipVol]];
					minCovA <- covmat[[minEllipVol]];
					mdistA <- mdistlist[[minEllipVol]];
					for ( palpha in alphalist2 )
					{
						if (palpha >= alpha)
						{
							correction <- mdistA[h]/qchisq(palpha,p);
							minCov <- minCovA*correction;
#							cat("(",minEllipVol,":",h,":",palpha,":",correction,":",det(minCov),")");
							mdist <- try(mahalanobis(intdata,mincentroid,minCov));
							if (!inherits(mdist, "try-error"))
							{
								inside <- (mdist < chithreshold);
								if (!is.na(sum(inside)))
								{
									ptsinside <- sum(inside)
									if (ptsinside >= h0)
									{
										auxdata <- intdata[inside,];
										newCentroid <- apply(auxdata,2,mean);
										newCovariance <- cov(auxdata);
										distanceupdate <- try(mahalanobis(auxdata,newCentroid,newCovariance));
										if (!inherits(distanceupdate, "try-error"))
										{
											dstamples <- ptsinside;
											if (dstamples > 250)
											{	
												dstamples <- 250; # no more than 100 samples for p-value estimation
												distanceupdate <- distanceupdate[sample(nrow(auxdata),dstamples)];
											}
											dsample <- p.threshold*( 1:dstamples - 0.5)/dstamples;
											distanceupdate <- distanceupdate[order(distanceupdate)];
											correction <- distanceupdate[dstamples]/chithreshold;
											disTheoretical <- qchisq(dsample,p);
											kst <- ks.test(disTheoretical,distanceupdate/correction + rnorm(dstamples,0,1e-10));
											if ((kst$statistic < 0.75*minD) || (kst$p.value > maxp))
#											if (kst$p.value > maxp)
											{
				#								plot(disTheoretical,distanceupdate);
												 bcorrection <- correction;
												 bestmean[[k]] <- newCentroid;
												 bestCov[[k]] <- newCovariance*correction;
												 pvals[k] <- kst$p.value;
												 if (minD > kst$statistic) 
												 {
													minD <- kst$statistic;
												 }
												 if (kst$p.value > maxp)
												 {
													maxp <- kst$p.value;
												 }
												 atalpha <- palpha;
												 optsinside <- ptsinside;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			cat("/");
#			if (verbose) 
#			{
#				cat(sprintf("%3d : %5.3f : %5.3f : %5.3f : %5d->",cycles,atalpha,maxp,bcorrection,optsinside));
#			}
			inside.centroid <- 0;
## Check for cluster overlap
			if ((k > 1) && (maxp > minpvalThr))
			{
				for (i in 1:(k-1))
				{
					inside.centroid <- inside.centroid + 1.0*(mahalanobis(bestmean[[i]],bestmean[[k]],bestCov[[k]]) < chithreshold3) + 1.0*(mahalanobis(bestmean[[k]],bestmean[[i]],bestCov[[i]]) < chithreshold3);
				}
			}
			inside <- numeric(0);
			## Include cluster only if p value is significant and no overlap with already discovered clusters
			if ((maxp <= minpvalThr) || (inside.centroid > 0))
			{
				cat("-");
				cycles <- cycles + 1;
				bestmean[[k]] <- NULL;
				bestCov[[k]] <- NULL;
				robCov[[k]] <- NULL;
				k <- k - 1;
			}
			else
			{
				robCov[[k]] <- list(centroid=bestmean[[k]],cov=bestCov[[k]]);
				mdist <- mahalanobis(thedataCpy,bestmean[[k]],bestCov[[k]]);
				inside <- (mdist < chithreshold);
				cludata <- thedataCpy[inside,];
				if (nrow(cludata) >= h0)
				{
					bestmean[[k]] <- apply(cludata,2,mean);
					bestCov[[k]] <- cov(cludata);					
					mdist <- mahalanobis(cludata,bestmean[[k]],bestCov[[k]]);
					correction <- max(mdist)/chithreshold;
					bestCov[[k]] <- bestCov[[k]]*correction;
					robCov[[k]] <- list(centroid=bestmean[[k]],cov=bestCov[[k]]);
#					cat(sprintf("|(%4.2f) S:%4d|",correction,nrow(cludata)));
					mdist <- mahalanobis(intdata,bestmean[[k]],bestCov[[k]]);
					inside <- (mdist < chithreshold3);
					ndata <- sum(!inside);
					if (ndata >= h0)
					{
						intdata <- intdata[!inside,];
						maxMahadis <- 1*(mdist[!inside] < chithreshold_out);
					}
					else
					{	
						ndata <- 0;
						maxMahadis <- rep(1,nrow(intdata));
					}
				}
				else
				{
					cycles <- cycles + 1;
					bestmean[[k]] <- NULL;
					bestCov[[k]] <- NULL;
					robCov[[k]] <- NULL;
					k <- k - 1;
				}
			}
		}
		else
		{
			cycles <- cycles + 1;
			bestmean[[k]] <- NULL;
			bestCov[[k]] <- NULL;
			robCov[[k]] <- NULL;
			k <- k - 1;
			ndata <- nrow(intdata);
			maxMahadis <- numeric(ndata);
		}
		p.threshold <- 0.99*p.threshold;
		minpvalThr <- minpvalThr/10.0;
		andata <- ndata;
		if (verbose) 
		{
			cat(sprintf("%3d:%3d:( %4d -> %4d ) a: %5.3f p: %5.3f D: %5.3f C: %5.3f K: %4d In: %4d \n",k,inside.centroid,ondata,andata,atalpha,maxp,minD,correction,JClusters,sum(inside)));
		}
		else 
		{
			cat("]");
		}
	}
	k <- length(bestmean);
	## assign clusters labels to all points
	if (k > 0)
	{
		ClusterLabels <- nearestCentroid(dataset,bestmean,bestCov,0);
	}
	cat("(",k,")}");
	
	result <- list(
		cluster = ClusterLabels,
		classification = ClusterLabels,
		centers = bestmean,
		covariances = bestCov,
		robCov = robCov,
		pvals = pvals,
		k = k,
		features = features
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
