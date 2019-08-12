#' returns the jaccard index matrix of two labeled sets
#'
#' jaccard = intersection/union
#' @param clustersA cluster labels of set 
#' @param clustersB second cluster labels of the same set
#' @return Returns the jaccard matrix of the labeled set
#' @examples
#' @importFrom 
#' @export

jaccardMatrix <- function(clustersA=NULL,clustersB=NULL)
{
	jaccardM <- NULL;
	jaccardA <- NULL;
	meanJaccard <- NULL;
	if (length(clustersA) == length(clustersB))
	{
		minclassA <- min(clustersA)
		minclassB <- min(clustersB)
		maxlabelA <- max(clustersA);
		maxlabelB <- max(clustersB);
		jaccardM <- matrix(0,nrow=(maxlabelA-minclassA+1),ncol=(maxlabelB-minclassB)+1);
		ii <- 1;
		for (i in minclassA:maxlabelA)
		{
			setA <- (clustersA == i);
			jj <- 1;
			for (j in minclassB:maxlabelB)
			{
				setB <- (clustersB == j);
				unionsum <- sum(setA | setB);
				if (unionsum > 0) jaccardM[ii,jj] <- sum(setA & setB)/unionsum;
				jj <- jj + 1;
			}
			ii <- ii + 1;
		}
		meanJaccard <- mean(c(apply(jaccardM,2,max),apply(jaccardM,1,max)));
		rownames(jaccardM) <- minclassA:maxlabelA
		colnames(jaccardM) <- minclassB:maxlabelB
		jaccardA <- clustersA;
		classA <- as.character(clustersA);
		classB <- as.character(clustersB);
		for (i in 1:length(clustersA)) 
		{
			jaccardA[i] <- jaccardM[classA[i],classB[i]];
		}
	}
	result <- list(jaccardMat=jaccardM,elementJaccard=jaccardA,balancedMeanJaccard=meanJaccard);
	
	return (result)
}