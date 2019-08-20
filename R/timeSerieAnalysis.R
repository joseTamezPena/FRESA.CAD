timeSerieAnalysis <-
function(variableList,baseModel,data,timevar="time",contime=".",Outcome=".",...,description=".",Ptoshow = c(1),plegend= c("p"),timesign="-",catgo.names=c("Control", "Case"))
{
# it will do the time analysis for all variables in the var list. It will create a data frame with the results of the analysis
# the analysis will  be done using GLS (nlme package) 


if (!requireNamespace("nlme", quietly = TRUE)) {
   install.packages("nlme", dependencies = TRUE)
} 

	frm1 <- NULL

	vnames <- as.vector(variableList[,1]);
	if (description != ".")
	{
		plottitles <- as.vector(variableList[,description]);
	}
	else
	{
		plottitles <- vnames;
	}
	sigma <- vector();
	size = length(vnames);
	anames <- vector();


	if (contime==".")
	{
	  contime=timevar;
	}

	
	if (timesign=="-")
	{
		timeorder <- data[order(-data[,contime]),];
	}
	else
	{
		timeorder <- data[order(data[,contime]),];
	}

	t <- tapply(timeorder[,contime], timeorder[,timevar], mean,na.rm=TRUE);

	farctionlowess <- min(0.5,2.0/(length(table(data[,timevar]))));

	if (Outcome != ".")
	{
		timeorderCases <- subset(timeorder,get(Outcome) == 1);
		timeorderControl <- subset(timeorder,get(Outcome) == 0);

		tCase <- tapply(timeorderCases[,contime], timeorderCases[,timevar], mean,na.rm=TRUE)
		tControl <- tapply(timeorderControl[,contime], timeorderControl[,timevar], mean,na.rm=TRUE)
	}


	for (j in 1:size)
	{
    
		frm1 <- paste(vnames[j]," ~ ",baseModel);
		cat(frm1,"\n");
		
		obj_s <- eval(parse(text=paste("try(nlme::gls(formula(",frm1,"),data,na.action=na.exclude,...))")))
		if ( inherits(obj_s, "try-error"))
		{
			cat("Getting the gls without parameters\n")
			obj_s <- eval(parse(text=paste("try(nlme::gls(formula(",frm1,"),data,na.action=na.exclude))")))
		}

		predCases <- NULL
		predControl <- NULL
		reg <- NULL
		if ( !inherits(obj_s, "try-error"))
		{
			reg <- summary(obj_s);	
			if (Outcome != ".")
			{
				predCases <- predict(obj_s,timeorderCases,na.action=na.exclude)
				predControl <- predict(obj_s,timeorderControl,na.action=na.exclude)
			}
		}

	mval <- tapply(data[,vnames[j]], data[,timevar], mean,na.rm=TRUE)
	sdval <- tapply(data[,vnames[j]], data[,timevar], sd,na.rm=TRUE)

    mval[is.na(mval)] <- 0.0;
    sdval[is.na(sdval)] <- 0.01;

	size <- tapply(data[,vnames[j]], data[,timevar], length)
	delta <- sdval / sqrt( size)
	miny = min(mval-1.5*sdval);
	maxy = max(mval+1.5*sdval);

		

	if (Outcome != ".")
	{

		mvalc <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], mean,na.rm=TRUE)
		sdvalc <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], sd,na.rm=TRUE)
       
		mvalo <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], mean,na.rm=TRUE)
		sdvalo <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], sd,na.rm=TRUE)

        mvalc[is.na(mvalc)] <- 0.0;
        sdvalc[is.na(sdvalc)] <- 0.01;
        mvalo[is.na(mvalo)] <- 0.0;
        sdvalo[is.na(sdvalo)] <- 0.01;

		miny = min(c(mvalc-1.5*sdvalc,mvalo-1.5*sdvalo));
		maxy = max(c(mvalc+1.5*sdvalc,mvalo+1.5*sdvalo));


	}
	
	if (miny<maxy)
	{

			maxtime = max(t);
			mintime = min(t);
			deltatime = maxtime-mintime;
				
			if (Outcome != ".")
			{

			  
		  
				if (timesign=="-")
				{
					mval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], mean,na.rm=TRUE)
					sdval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], sd,na.rm=TRUE)
					size <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], length)
					delta <- sdval / sqrt( size)
					errbar(-tCase,mval,mval-delta,mval+delta,col="red",ylab=vnames[j],xlab="time",type="b",ylim=c(miny,maxy),xlim=c(mintime,maxtime))
					lines(lowess(-timeorderCases[,contime],predCases,f = farctionlowess),col="red",lty=3)
		   
					mval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], mean,na.rm=TRUE)
					sdval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], sd,na.rm=TRUE)
					size <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], length)
					delta <- sdval / sqrt( size)
					errbar(-tControl,mval,mval-delta,mval+delta,add=TRUE,col="blue",type="b")
		  
					lines(lowess(-timeorderControl[,contime],predControl,f = farctionlowess),col="blue",lty=2)
#					legend(-maxtime+0.1*deltatime,miny + 0.2*(maxy-miny), catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
					legend("bottomleft", catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
				
					if (!is.null(reg))
					{
						for (lg in 1:length(Ptoshow))
						{
							legend(-maxtime+0.65*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),sprintf("t(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],3]),bty="n");
							if (reg$tTable[Ptoshow[lg],4]<0.0001)
							{
								legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"****",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.001)
								{
									legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"***",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.01)
									{
										legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"**",bty="n");
									}
									else
									{
										if (reg$tTable[Ptoshow[lg],4]<0.05)
										{
											legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"*",bty="n");
										}
									}
								}
							}
						}
					}
				}
				else
				{
					mval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], mean,na.rm=TRUE)
					sdval <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], sd,na.rm=TRUE)
					size <- tapply(timeorderCases[,vnames[j]], timeorderCases[,timevar], length)
					delta <- sdval / sqrt( size)

					errbar(tCase,mval,mval-delta,mval+delta,col="red",ylab=vnames[j],xlab="time",type="b",ylim=c(miny,maxy),xlim=c(mintime,maxtime))
					lines(lowess(timeorderCases[,contime],predCases,f = farctionlowess),col="red",lty=3)

					mval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], mean,na.rm=TRUE)
					sdval <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], sd,na.rm=TRUE)
					size <- tapply(timeorderControl[,vnames[j]], timeorderControl[,timevar], length)
					delta <- sdval / sqrt( size)

					errbar(tControl,mval,mval-delta,mval+delta,add=TRUE,col="blue",type="b")
		  
					lines(lowess(timeorderControl[,contime],predControl,f = farctionlowess),col="blue",lty=2)
#					legend(maxtime-0.9*deltatime,miny + 0.2*(maxy-miny), catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
					legend("bottomleft", catgo.names, col=c("blue","red"), lty = 2:3,bty="n")
			
					if (!is.null(reg))
					{
						for (lg in 1:length(Ptoshow))
						{
							legend(maxtime-0.45*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),sprintf("t(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],3]),bty="n");
							if (reg$tTable[Ptoshow[lg],4]<0.0001)
							{
								legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"****",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.001)
								{
									legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"***",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.01)
									{
										legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"**",bty="n");
									}
									else
									{
										if (reg$tTable[Ptoshow[lg],4]<0.05)
										{
											legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"*",bty="n");
										}
									}
								}
							}
						}
					}
				}
				
		  
			}
			else
			{
				if (timesign=="-")
				{			
				  errbar(-t,mval,mval-delta,mval+delta,ylab=vnames[j],xlab="time",type="b",ylim=c(miny,maxy),xlim=c(mintime,maxtime))
				  lines(lowess(-timeorder[,contime],predict(obj_s,timeorder,na.action=na.exclude),f = farctionlowess),col="blue",lty=2)      
				  for (lg in 1:length(Ptoshow))
					{
						legend(-maxtime+0.65*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),sprintf("t(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],3]),bty="n");
							if (reg$tTable[Ptoshow[lg],4]<0.0001)
							{
								legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"****",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.001)
								{
									legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"***",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.01)
									{
										legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"**",bty="n");
									}
									else
									{
										if (reg$tTable[Ptoshow[lg],4]<0.05)
										{
											legend(-maxtime+0.9*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"*",bty="n");
										}
									}
								}
							}
					}
				}
				else
				{
					errbar(t,mval,mval-delta,mval+delta,ylab=vnames[j],xlab="time",type="b",ylim=c(miny,maxy),xlim=c(mintime,maxtime))
					lines(lowess(timeorder[,contime],predict(obj_s,timeorder,na.action=na.exclude),f = farctionlowess),col="blue",lty=2)      
					for (lg in 1:length(Ptoshow))
					{
						legend(maxtime-0.45*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),sprintf("t(%s) = %5.3f",plegend[lg],reg$tTable[Ptoshow[lg],3]),bty="n");
							if (reg$tTable[Ptoshow[lg],4]<0.0001)
							{
								legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"****",bty="n");
							}
							else
							{
								if (reg$tTable[Ptoshow[lg],4]<0.001)
								{
									legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"***",bty="n");
								}
								else
								{
									if (reg$tTable[Ptoshow[lg],4]<0.01)
									{
										legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"**",bty="n");
									}
									else
									{
										if (reg$tTable[Ptoshow[lg],4]<0.05)
										{
											legend(maxtime-0.1*deltatime,miny + (0.35-0.07*lg)*(maxy-miny),"*",bty="n");
										}
									}
								}
							}
					}
				}
			}
			
			title(main=plottitles[j])
			if (!is.null(reg))
			{
				if (j>1) 
				{
					coeff <- rbind (coeff,as.vector(reg$tTable[,1]));
					std.Error <- rbind (std.Error,as.vector(reg$tTable[,2]));
					t.value <- rbind (t.value,as.vector(reg$tTable[,3]));
					p.value <- rbind (p.value,as.vector(reg$tTable[,4]));
				}
				else
				{
					coeff <- rbind (as.vector(reg$tTable[,1]));
					std.Error <- rbind (as.vector(reg$tTable[,2]));
					t.value <- rbind (as.vector(reg$tTable[,3]));
					p.value <- rbind (as.vector(reg$tTable[,4]));
				}
				sigma <- append(sigma,reg$sigma);
				anames <- append(anames,vnames[j]);
			}
		}
		
	}

	rownames(coeff) <- anames;
	rownames(std.Error) <- anames;
	rownames(t.value) <- anames;
	rownames(p.value) <- anames;
	names(sigma) <- rownames(p.value);

	if (!is.null(reg))
	{
		colnames(coeff) <- rownames(reg$tTable);
		colnames(std.Error) <- rownames(reg$tTable);
		colnames(t.value) <- rownames(reg$tTable);
		colnames(p.value) <- rownames(reg$tTable);
	}


	result  <- list(coef=coeff,
	std.Errors=std.Error,
	t.values=t.value,
	p.values=p.value,
	sigmas=sigma,
	lastfit=obj_s);
	return (result);
}

trajectoriesPolyFeatures <- function(data,feature="v1", degree=2, time="t", group="ID",timeOffset=0,strata=NULL,...)
{
  aids <- data[,group]
  ids <- unique(aids)
  coefs <- as.data.frame(matrix(0,nrow = length(ids),ncol = degree + 1));
  rownames(coefs) <- ids;
  miny <- min(data[,feature],na.rm = TRUE);
  maxy <- max(data[,feature],na.rm = TRUE);
  minx <- min(data[,time],na.rm = TRUE);
  maxx <- max(data[,time],na.rm = TRUE);
  plot(1,type="p",xlim=c(minx, maxx), ylim=c(miny, maxy),...)
  range <- (maxx-minx);
  timesamples <- minx + ((0:100)/100.0)*range;
  dataTime <- as.data.frame(cbind(0:100,timesamples));
  for (i in 1:length(ids))
  {
    whoid <- aids == ids[i];
    coefs[i,] <- rep(NA,1 + degree);
    if (sum(whoid) > 0)
    {
      minv = min(data[whoid,time],na.rm = TRUE);
      maxv = max(data[whoid,time],na.rm = TRUE);
      if ((minv < maxv) && (minv <= timeOffset) && (maxv >= timeOffset))
      {
        colnames(dataTime) <- c(ids[i],"t");
        fd <- as.data.frame(cbind(xs = data[whoid,feature],t = data[whoid,time]));
        fitlm <- try(lm(paste("xs ~ poly(I(t -",timeOffset,"),degree = ",degree,", raw=TRUE)"),data = fd,na.action = na.omit));
        if (!inherits(fitlm, "try-error"))
        {
          coefs[i,] <- fitlm$coefficients;
          points(data[whoid,time],data[whoid,feature],type="p",col="red");
          lines(timesamples,predict(fitlm,dataTime,na.action=na.exclude),col="blue",lty=2)
        }
      }
    }
  }
  return(coefs)
}