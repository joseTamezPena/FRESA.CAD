barPlotCiError <-  function(ciTable, metricname, thesets, themethod, main, angle=0, offsets=c(0.1,0.1),scoreDirection=">",ho=NULL, ...)
{
  op <- par(no.readonly=TRUE);

error.bar <- function(x, y, upper, lower=upper, length=0.03, ...){
  if (length(x) != length(y) | length(y) != length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y + upper, x, y - lower, angle = 90, code = 3, length = length, ...)
}

  colnames(ciTable) <- c(metricname,"lower","upper")
  rownames(ciTable) <- rep(thesets,length(themethod))
  if (ciTable[1,3] < ciTable[1,2])
  {
    upper <- ciTable[,2]
    ciTable[,2] <- ciTable[,3]
    ciTable[,3] <- upper
  }
  ciTable <- ciTable[order(rep(1:length(thesets),length(themethod))),]
  
#  pander::pander(ciTable,caption = main,round = 3)
  barmatrix <- matrix(ciTable[,1],length(themethod),length(thesets))
  upmatrix <- matrix(ciTable[,3],length(themethod),length(thesets))
  infmatrix <- matrix(ciTable[,2],length(themethod),length(thesets))
  colnames(barmatrix) <- thesets
  rownames(barmatrix) <- themethod
  colnames(upmatrix) <- thesets
  rownames(upmatrix) <- themethod
  colnames(infmatrix) <- thesets
  rownames(infmatrix) <- themethod
#  pander::pander(barmatrix,caption = main,round = 3)

	SupMethod <- barmatrix;
	InfMethod <- barmatrix;
	zeroMethod <- barmatrix;
	
	for (i in 1:nrow(barmatrix))
	{
		for (j in 1:ncol(barmatrix))
		{
			SupMethod[i,j] <- sum( ( barmatrix[i,j] > upmatrix  ) & ( infmatrix[i,j] > barmatrix ) );
			InfMethod[i,j] <- sum( ( barmatrix[i,j] < infmatrix  ) & ( upmatrix[i,j] < barmatrix ) );
			zeroMethod[i,j] <- 0.0 - ( infmatrix[i,j] > 0 )*(1+sum( barmatrix[i,j] > upmatrix )) - ( upmatrix[i,j] < 0 )*(1+sum( barmatrix[i,j] < infmatrix ));
		}
	}

	if (scoreDirection == ">") 
	{
		interMethodScore <- SupMethod - InfMethod;
		if (!is.null(ho))
		{
			interMethodScore <- SupMethod*(infmatrix > ho ) - InfMethod - (infmatrix < ho );
		}
	}
	if (scoreDirection == "<") 
	{
		interMethodScore <- InfMethod-SupMethod;
		if (!is.null(ho))
		{
			interMethodScore <- InfMethod*( upmatrix < ho ) - SupMethod - ( upmatrix > ho );
		}
	}
	if (scoreDirection == "==") 
	{
		interMethodScore <- zeroMethod;
	}

	cmm <- colMeans(interMethodScore)
	rmm <- rowMeans(interMethodScore)
	if (scoreDirection == ">")
	{
		cmm <- cmm + 1.0e-10*colMeans(barmatrix)
		rmm <- rmm + 1.0e-10*rowMeans(barmatrix)
	}
	else
	{
		cmm <- cmm - 1.0e-10*colMeans(barmatrix)
		rmm <- rmm - 1.0e-10*rowMeans(barmatrix)
	}

	interMethodScore <- interMethodScore[order(-rmm),order(-cmm)];
	SupMethod <- SupMethod[order(-rmm),order(-cmm)];
	InfMethod <- InfMethod[order(-rmm),order(-cmm)];
	barmatrix <- barmatrix[order(-rmm),order(-cmm)]
	upmatrix <- upmatrix[order(-rmm),order(-cmm)]
	infmatrix <- infmatrix[order(-rmm),order(-cmm)]
	themethod <- themethod[order(-rmm)]
	thesets <- thesets[order(-cmm)];
  
  ymin <- min(0,1.25*min(ciTable));
  ymax <- max(0,1.10*max(ciTable));
  
#  par(mfrow = c(1,1));
  par(pty="m")

  if (length(thesets) > 1)
  {
	nf <- layout(matrix(c(1,3,1,3,1,4,1,5,1,5,1,6,2,0), 7, 2, byrow = TRUE),widths = c(3,1),heights = c(2,2,1.5,2,2,1.5,3))
  }
  else
  {
    nf <- layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights = c(5,1));
  }
  mar = par("mar")
  par(mar = c(0.0, mar[2],mar[3], mar[4]))
  barpo <- barplot(barmatrix,cex.names=0.7,las=2,ylim = c(ymin,ymax),main=main,ylab=metricname,beside= (length(thesets) > 1),legend = themethod, xaxt="n",...)
  
  xpos <- t(barpo[as.integer(length(themethod)/2) + 1,])
    
  error.bar(barpo,as.vector(barmatrix),as.vector((upmatrix - barmatrix)),as.vector((barmatrix-infmatrix)))

	par(mar = c(2.0, mar[2], 0.65, mar[4]))
	col = heat.colors(1+max(interMethodScore)-min(interMethodScore));
	barplot(interMethodScore,las=2,cex.axis=0.75,beside = TRUE,axisnames = FALSE,ylab="Score",col=col[as.vector(interMethodScore)-min(interMethodScore)+1]);
  #text(cex=0.7, x=xpos, y=-0.1, thesets, xpd=TRUE, srt=45)
  text(cex=0.8, x = xpos - offsets[1], y = min(interMethodScore)-offsets[2], thesets, xpd = TRUE, srt = angle)
  if (length(thesets) > 1)
  {
	par(mar = c(0.0, mar[2],mar[3], mar[4]))
	barp <- barplot(colMeans(barmatrix),cex.axis=0.95,cex.names=0.95,las=2,ylim = c(ymin,ymax),ylab=metricname, xaxt="n",cex.lab=0.85)
    xpos <- barp;
	par(mar = c(2.0, mar[2], 0.65, mar[4]))
	colmax <- apply(interMethodScore,2,max,na.rm = TRUE);
	barplot(colmax,cex.axis=0.95,cex.names=0.75,las=2,ylab="Score",xlab="Sets",col=col[colmax-min(interMethodScore)+1],cex.lab=0.85);

	par(mar = c(0.0, mar[2],mar[3], mar[4]))
	barp <- barplot(rowMeans(barmatrix),cex.axis=0.95,cex.names=0.95,las=2,ylim = c(ymin,ymax),ylab=metricname, xaxt="n",cex.lab=0.85)
    xpos <- barp;
	par(mar = c(2.0, mar[2], 0.65, mar[4]));
	rowmax <- apply(interMethodScore,1,max,na.rm = TRUE);
	barplot(rowmax,cex.axis=0.95,cex.names=0.75,las=2,ylab="Score",xlab="Methods",col=col[rowmax-min(interMethodScore)+1],cex.lab=0.85);

 }
  par(mar=mar);
#  par(mfrow = c(1,1),mar=mar);
  par(op)

  
  return(list(barplot=barpo,ciTable = list(mean=barmatrix,low95=infmatrix,top95=upmatrix), supMethod = SupMethod, infMethod = InfMethod, interMethodScore=interMethodScore ))
}
