heatMaps <- function(variableList=NULL,varRank=NULL,Outcome,data,title="Heat Map",hCluster=FALSE,prediction=NULL,Scale=FALSE,theFiveColors=c("blue","cyan","black","yellow","red"),outcomeColors = c("blue","lightgreen","yellow","orangered","red"),transpose=FALSE,...) 
{
  
  opg <- par(no.readonly=TRUE)
  
  color.bar <- function(lut, min, max=-min, nticks=3, ticks=seq(min, max, len=nticks),transpose=FALSE, Title='',...) 
  {
    op <- par(no.readonly=TRUE)
    par(new=TRUE,cex=0.35,mai=c(0,0,0,0),...)
    plot(c(0,5), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='',main='')
    
    dy <- (grconvertY(1.00,"npc")-grconvertY(0.00,"npc"))/(length(lut));
    axis(2, ticks, las=1,bty='n',outer=TRUE,labels=ifelse(ticks >= 1, sprintf("%0.1f", ticks),
sprintf("%0.1f", ticks)))
    for (i in 1:length(lut)) 
    {
      y = dy*(i-1)+min;
      rect(0,y,grconvertX(0.75,"npc"),(y+dy), col=lut[i], border=NA)
    }
    title(Title, line = 0.25,outer=TRUE)
    par(op)
  }
  
  Rowv = TRUE;
  Colv = TRUE;
  dendrogram = "both";
  if (hCluster == "col")
  {
		hCluster=FALSE;
		dendrogram="col";
		Rowv=NULL;
  }
  if (hCluster == "row")
  {
		hCluster=FALSE;
		dendrogram="row";
		Colv=NULL;
  }
  
  par(new=FALSE,pty='m')
  if (!requireNamespace("gplots", quietly = TRUE)) {
    install.packages("gplots", dependencies = TRUE)
  } 
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer", dependencies = TRUE)
  } 
  
  # creates a own color palette from red to blue
  my_palette <- colorRampPalette(theFiveColors)(n = 224)
  #	my_palette <- colorRampPalette(c("blue","light blue","black","pink","red"))(n = 224)
  
  # defines the color breaks manually for a "skewed" color transition
  col_breaks = c(
    seq(-2.0,-1.101,length=50),  		# for blue
    seq(-1.10,-0.201,length=50),          # for cyan
    seq(-0.20,0.20,length=25),          # for black
    seq(0.201,1.10,length=50),            # for yellows
    seq(1.101,2.0,length=50))             # for red
  
  
  if (is.null(variableList))
  {
	variableList <- colnames(data)
	names(variableList) <- colnames(data)
  }
  if (class(variableList)=="data.frame")
  {
    vnames <- as.vector(variableList[,1]);
  }
  else
  {
    vnames <- names(variableList[!(names(variableList) %in% Outcome)]);
  }
  if (length(vnames)>2)
  {
	  frm <- paste(Outcome," ~ ",1);
	  added = 0;
	  for (i in 1:length(vnames))
	  {
		if ((vnames[i] != "")&&(vnames[i] != "1"))
		{
		  frm <- paste(frm,"+",vnames[i])
		  added = added + 1;
		}
	  }
	  modelFrame <- model.frame(formula(frm),data);
	  rownames(modelFrame) <- rownames(data)
	  
	  if (is.null(varRank)) 
	  {
		topvarID <- seq(1, added, 1)
		hits <- seq(1, added, 1)
	  }
	  else
	  {
		topvarID <- as.numeric(rownames(varRank));
		hits <- as.vector(as.matrix(varRank))
	  }
	  orderData <- cbind(data[,Outcome]);
	  orderData <- as.data.frame(orderData);
	  rownames(orderData) <- rownames(data);

	  colnames(orderData) <- c(Outcome);
	  cn = 1; 
	  if (!is.null(prediction)) 
	  {
		if (prediction[1]!=Outcome)
		{
		  cn = cn + 1;
		  orderData <- cbind(orderData,as.vector(prediction));
		  colnames(orderData)[cn] <- "prediction";
		}
	  }
	  for ( i in 1:length(topvarID))
	  {
		if (hits[i]>0)
		{
		  orderData <- cbind(orderData,modelFrame[,topvarID[i]+1])
		  cn = cn + 1;
		  colnames(orderData)[cn] <- colnames(modelFrame)[topvarID[i]+1];
		}
	  }
	  
	  dataMat <- as.matrix(orderData);
	  minout <- min(data[,Outcome]);
	  maxout <- max(data[,Outcome]);
	  
	  
	  if (!is.null(prediction))
	  {
		if (prediction[1]!=Outcome)
		{
		  if (transpose)
		  {
			orderData <- eval(parse(text=paste("with(orderData,orderData[order(",Outcome,",orderData[,2]),])")))
			
		  }
		  else
		  {
			orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,",-orderData[,2]),])")))
		  }
		}
	  }
	  else
	  {
		sortsum <- rowSums (abs(orderData), na.rm = TRUE)
		orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,",-sortsum),])")))
		
		cm <- as.data.frame(cor(t(orderData[,-1]),method="pearson"));
		cm[is.na(cm)] <- -2;
#		save(cm,file="cormat.RDATA");


		rownames (cm) <- c(1:nrow(cm))
		colnames (cm) <- c(1:nrow(cm))
		indx <- 1
		ix <- 1
		rm <- cm;
		thenames <- colnames(cm)
		for (i in 1:(nrow(cm)-2))
		{
		  rm <- cm[ix,-indx]
		  ix <- as.integer(names(which.max(rm)))
		  indx <- append(indx,ix)
		}
		indx <- append(indx,as.integer(thenames[-indx]))
	#    cat(indx,"\n")
		if (transpose)
		{
		  orderData <- orderData[rev(indx),];
		  orderData <- eval(parse(text=paste("with(orderData,orderData[order(",Outcome,"),])")))
		}
		else
		{
		  orderData <- orderData[indx,];
		  orderData <- eval(parse(text=paste("with(orderData,orderData[order(-",Outcome,"),])")))
		}
	  }

	  
	  orderData <- as.matrix(orderData);
	  
	  nr <- nrow(orderData);
	  index <- floor(nr*(orderData[,1]-minout)/(1.0001*maxout-minout))+1;
	  rowcolors <- colorRampPalette(outcomeColors)(n = nr)[index];
	  
	  orderData <- orderData[,-1];
	  if (class(Scale) == "logical")
	  {
		if (Scale) orderData <- scale(orderData);
	  }
	  nclass <- length(table(data[,Outcome]))
	  nticks <- min(5,nclass);

#	  print(nclass)
	  
	  colorange = FALSE;
	  if (class(Scale) == "logical")
	  {
		colorange = !Scale;
	  }
	  if (class(Scale) == "numeric")
	  {
		colorange = TRUE;
	  }
	  if (colorange)
	  {
			if (class(Scale) == "numeric")
			{
				kmin <- Scale[1];
				kmax <- Scale[2];
			}
			else
			{
				kmin <- min(orderData);
				kmax <- max(orderData);
			}
			brange <- (kmax-kmin)/5;
			  col_breaks = c(
				seq(kmin,kmin+0.99*brange,length=45),  		# for blue
				seq(kmin+brange,kmin+1.99*brange,length=45),          # for cyan
				seq(kmin+2*brange,kmin+2.99*brange,length=45),          # for black
				seq(kmin+3*brange,kmin+3.99*brange,length=45),            # for yellows
				seq(kmin+4*brange,kmax,length=45))             # for red
			
		}

	  
	  if (transpose)
	  {
		if (is.null(prediction))
		{
		  if (hCluster==TRUE)
		  {
			heatMap <- gplots::heatmap.2(t(orderData), 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram="both", 		# raw and column dendrogram
										 ColSideColors=rowcolors,
										 #                                     lmat=rbind(c(0,5,0,4,0), c(0,3,1,2,0)),
										 #                                     lhei=c(1.5,4.5),
										 #                                     lwid=c(0.10,3,3,2,0.10),
										 ...) 				
		  }
		  else
		  {
			#			rownames(orderData) <- format(orderData[,1],digits=4);
			
			heatMap <- gplots::heatmap.2(t(orderData), 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 ColSideColors=rowcolors,
										 dendrogram=dendrogram, 		
										 Colv=Colv,
										 Rowv=Rowv,
										 ...) 				# turn off row clustering				
		  }
		}
		else
		{
		  if (hCluster==FALSE)
		  {
			heatMap <- gplots::heatmap.2(t(orderData), 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram=dendrogram, 		
										 Colv=Colv,
										 Rowv=Rowv,
										 ColSideColors=rowcolors,
										 ...) 				
		  }
		  else
		  {
			#			rownames(orderData) <- format(orderData[,1],digits=4);
			#        index <- ((orderData[,1]+2)/4)
			#        index <- index*(index>0);
			#        index <- index*(index<1)+0.999*(index>=1);
			#        index <- floor(nrow(orderData)*index)+1;			
			#        rowcolors <- colorRampPalette(theFiveColors)(n = nrow(orderData))[index];
			#        orderData <- orderData[,-1];
			heatMap <- gplots::heatmap.2(t(orderData), 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram="col", 		# 
										 ColSideColors=rowcolors,
										 #                                     lmat=rbind(c(0,5,0,4,0), c(0,3,1,2,0)),
										 #                                     lhei=c(1.5,4.5),
										 #                                     lwid=c(0.10,2,0.25,6.55,0.10),
										 Rowv="NA",				# turn off col clustering
										 ...)
		  }
		}
		rowcolors <- colorRampPalette(outcomeColors)(n = min(100,nclass));
		color.bar(rowcolors,min=minout,max=maxout,nticks=nticks,Title=Outcome,omd=c(0.96,0.99,0.83,0.95))
	  }
	  else
	  {
		if (is.null(prediction))
		{
		  if (hCluster==TRUE)
		  {
			heatMap <- gplots::heatmap.2(orderData, 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram="both", 		# raw and column dendrogram
										 RowSideColors=rowcolors,
										 #                                     lmat=rbind(c(0,5,0,4,0), c(0,3,1,2,0)),
										 #                                     lhei=c(1.5,4.5),
										 #                                     lwid=c(0.10,2,0.25,6.55,0.10),
										 ...) 				# turn off column clustering
		  }
		  else
		  {
			#			rownames(orderData) <- format(orderData[,1],digits=4);
			
			heatMap <- gplots::heatmap.2(orderData, 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram=dendrogram, 		
										 Colv=Colv,
										 Rowv=Rowv,
										 RowSideColors=rowcolors,
										 ...) 					
		  }
		}
		else
		{
		  if (hCluster==FALSE)
		  {
			heatMap <- gplots::heatmap.2(orderData, 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram=dendrogram, 		
										 Colv=Colv,
										 Rowv=Rowv,
										 RowSideColors=rowcolors,
										 ...) 				
		  }
		  else
		  {
			#			rownames(orderData) <- format(orderData[,1],digits=4);
			#        index <- ((orderData[,1]+2)/4)
			#        index <- index*(index>0);
			#        index <- index*(index<1)+0.999*(index>=1);
			#        index <- floor(nrow(orderData)*index)+1;			
			#        rowcolors <- colorRampPalette(theFiveColors)(n = nrow(orderData))[index];
			#        orderData <- orderData[,-1];
			heatMap <- gplots::heatmap.2(orderData, 
										 main = title, 			# heat map title
										 notecol="black",      	# change font color of cell labels to black
										 density.info="none",  	# turns off density plot inside color legend
										 trace="none",         	# turns off trace lines inside the heat map
										 margins =c(7,8),     	# widens margins around plot
										 col=my_palette,       	# use on color palette defined earlier 
										 breaks=col_breaks,    	# enable color transition at specified limits
										 dendrogram="row", 		# 
										 RowSideColors=rowcolors,
										 #                                     lmat=rbind(c(0,5,0,4,0), c(0,3,1,2,0)),
										 #                                     lhei=c(1.5,4.5),
										 #                                     lwid=c(0.10,2,0.25,6.55,0.10),
										 Colv="NA",				# turn off col clustering
										 ...)
		  }
		}
		rowcolors <- colorRampPalette(outcomeColors)(n = min(100,nclass));
		color.bar(rowcolors,min=minout,max=maxout,nticks=nticks,Title=Outcome,omd=c(0.04,0.07,0.07,0.20))
	  }
	  
	  par(opg)
	result <- list(dataMatrix=dataMat,
                 orderMatrix=orderData,
                 heatMap=heatMap)
	}
	else
	{
		cat("No enough features for heatmap\n");
		result <- NULL
	}

  #  par(pty='m',mfrow=c(1,1))
  
  return (result);
}
