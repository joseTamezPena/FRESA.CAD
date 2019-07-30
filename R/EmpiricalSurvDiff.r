EmpiricalSurvDiff<-function(times=times,status=status,groups=groups,samples=1000,type=c("SLR","Chi"),plots=FALSE,minAproxSamples=100,computeDist=FALSE,...)
{
  type <- match.arg(type);

LogRankNormal<-function(times=times,status=status,groups=groups)
{
  LR <-.Call("logRank",times,status,groups);
  return (LR)
} 


SLRNullDistribution<-function(times=times,status=status,groups=groups,samples=1000,type=c("Chi","SLR"))
{
  type <- match.arg(type);
  logranksamples<-.Call("SLRNullDistribution",times,status,groups,samples,1+1*(type=="Chi"));
  logranksamples <- logranksamples[order(logranksamples)];
  return (logranksamples)
}

SLRDistribution<-function(times=times,status=status,groups=groups,samples=1000,type=c("Chi","SLR"))
{
  type <- match.arg(type);
  logranksamples<-.Call("SLRDistribution",times,status,groups,samples,1+1*(type=="Chi"));
  logranksamples <- logranksamples[order(logranksamples)];
  return (logranksamples)
}

SLRpvalue <-function(LR=LR,nullLR=NULL,LRDist=NULL,type="",minAproxSamples=100)
{
  
  
  pvalue=1;
  pequal=1;
  psup=1;
  pinf=1;
  if (type=="Chi")
  {
    psup = NA;
    pinf = NA;
    LR <- LR$SLR^2;
    samples=length(nullLR);
	  if (minAproxSamples > (0.20*samples))
	  {
		minAproxSamples <- 0.20*samples;
	  }
	  if (minAproxSamples < 5)
	  {
		minAproxSamples <- 5; # at least five samples
	  }
	topcount <- sum(nullLR>LR)
	pvalue <- topcount/samples;
	if (topcount<minAproxSamples)
	{
		LRATTop <- nullLR[samples-minAproxSamples+1];
		qvalue <-  qchisq(1.0-minAproxSamples/samples,1);
		copvalue <- 1.0 - pchisq(LR*qvalue/LRATTop,1);
#		cpvalue <- 1.0 - pchisq(LR, 1);
#		cat("Observed:",LRATTop,"Expected qvalue:",qvalue," psam:", pvalue," chip", cpvalue," cochip", copvalue,"\n");
		pvalue <- copvalue;
	}

    pequal = pvalue;
  }
  else
  {
    samples=length(nullLR);
	  if (minAproxSamples > (0.20*samples))
	  {
		minAproxSamples <- 0.20*samples;
	  }
	  if (minAproxSamples < 5)
	  {
		minAproxSamples <- 5; # at least five samples
	  }
    LR <- LR$SLR;
	
	

    topcount <- sum(nullLR > abs(LR));
    lowcount <- sum(nullLR <= -abs(LR));

    
    phigh = 0;
    if (topcount>minAproxSamples) 
    {
      phigh <- topcount/samples;
    }
    else
    {
		zvalue <- qnorm(minAproxSamples/samples);
      LRATTop <- nullLR[samples-minAproxSamples+1];
      zhighgain <- abs(LRATTop/zvalue);
      phigh <- pnorm(abs(LR),0,zhighgain,FALSE);
    }

    plow = 0;
    if (lowcount>minAproxSamples) 
    {
      plow <- lowcount/samples;
    }
    else
    {
		zvalue <- qnorm(minAproxSamples/samples);
      LRATbottom <- nullLR[minAproxSamples];
      zlowgain <- LRATbottom/zvalue;
      plow <- pnorm(-abs(LR),0,zlowgain,TRUE);
    }
    pequal = plow+phigh;
    if (2.0*sum(nullLR < LR) > samples)
    {
      pvalue = phigh;
    }
    else
    {
      pvalue = plow;
    }
    if (is.null(LRDist))
    {
      if (LR>0)
      {
        psup <- plow;
        pinf <- 1.0-psup;
      }
      else
      {
        pinf <- phigh;
        psup <- 1.0-pinf;
      }
    }
    else
    {
      LRsamples=length(LRDist);
	  if (minAproxSamples > (0.20*LRsamples))
	  {
		minAproxSamples <- 0.20*LRsamples;
	  }
	  if (minAproxSamples < 5)
	  {
		minAproxSamples <- 5; # at least five samples
	  }

	  topcount <- sum(LRDist > 0);
      lowcount <- sum(LRDist <= 0);

      phigh = 0;
      if (topcount > minAproxSamples) 
      {
        phigh <- topcount/LRsamples;
      }
      else
      {
		  mLR <- median(LRDist)
		  sdLR <- sd(LRDist)
		  zvalue <- abs(qnorm(minAproxSamples/LRsamples));
        LRATTop <- LRDist[LRsamples-minAproxSamples+1];
        sigma <- (LRATTop-mLR)/zvalue;
        phigh <- pnorm(0,mLR,sigma,FALSE);
      }
      plow = 0;
      if (lowcount > minAproxSamples) 
      {
        plow <- lowcount/LRsamples;
      }
      else
      {
		  mLR <- median(LRDist)
		  sdLR <- sd(LRDist)
		  zvalue <- abs(qnorm(minAproxSamples/LRsamples));
        LRATbottom <- LRDist[minAproxSamples];
        sigma <- (mLR-LRATbottom)/zvalue;
        plow <- pnorm(0,mLR,sigma,TRUE);
      }
      if (LR>0)
      {
        psup <- plow;
        pinf <- 1.0-psup;
      }
      else
      {
        pinf <- phigh;
        psup <- 1.0-pinf;
      }
    }
  }
  return (list(LR=LR,pvalue=pvalue,p.equal=pequal,p.sup=psup,p.inf=pinf))
}

  LRDist <- NULL;

  LR <- LogRankNormal(times,status,groups)
  nullLR <- SLRNullDistribution(times,status,groups,samples,type);
  if (computeDist) LRDist <- SLRDistribution(times,status,groups,samples,type);


  pvalues <- SLRpvalue(LR,nullLR,LRDist,type,minAproxSamples);
  if (plots)
  {
    smallnull <- nullLR
	LRp <- LR$SLR
	sLR <- "SLR"
    if (length(nullLR)>1000) smallnull <- nullLR[sample(samples,1000)]
    if (type=="SLR")
    {
      xlim <- c(-6,6)
    }
    else
    {
		sLR <- "Chi"
	 LRp <- LRp*LRp
      xlim <- c(0,36)
    }
    
    plot(ecdf(smallnull),sub=sprintf("%s: %8.2f, 1-p(|SLR|>0)= %8.3e, min(p)= %8.3e",sLR,LRp,pvalues$p.equal,pvalues$pvalue),xlim=xlim,col="black",lwd = 2,verticals = TRUE, do.points = FALSE,...);
    d=density(smallnull,bw="nrd")
    lines(d,col="blue")
    if (!is.null(LRDist))
    {
		smallnull <- LRDist
      if (length(smallnull)>1000) smallnull <- LRDist[sample(samples,1000)]
      lines(ecdf(smallnull),col="red",lwd = 1,verticals = TRUE, do.points = FALSE);  
      d=density(smallnull,bw="nrd")
      lines(d,col="pink")
    }
    if (type=="SLR")
    {
      abline(v=LR$SLR,col = "red");
    }
    else
    {
      abline(v=(LR$SLR)^2,col = "red");
    }
	legend("bottomright", c("CDF", "pdf","SLR"), col=c("black","blue","red"), lwd=c(2,1,1))

  }

  return (list(pvalue=pvalues$pvalue,LR=LR,p.equal=pvalues$p.equal,p.sup=pvalues$p.sup,p.inf=pvalues$p.inf,nullDist=nullLR,LRDist=LRDist))
}


