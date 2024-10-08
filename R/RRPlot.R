RRPlot <-function(riskData=NULL,
                  timetoEvent=NULL,
                  riskTimeInterval=NULL,
                  ExpectedPrevalence=NULL,
                  atRate=c(0.90,0.80),
                  atThr=NULL,
                  plotRR=TRUE,
                  title="",
                  ysurvlim=c(0,1.0))
{

  if (!requireNamespace("corrplot", quietly = TRUE)) {
    install.packages("corrplot", dependencies = TRUE)
  }

  OERatio <- NULL
  OE95ci <- NULL
  OAcum95ci <- NULL
  CumulativeOvs <- NULL
  DCA <- NULL
  OEData <- NULL
  OARatio <- NULL
  netBenefitatP <- NULL
  
  uvalues <- length(unique(riskData[,2]))
  isProbability <- (min(riskData[,2]) >= 0) && (max(riskData[,2]) <= 1.0) && (sd(riskData[,2]) > 0.00001)
  if (is.null(timetoEvent))
  {
    if (ncol(riskData)>2)
    {
      timetoEvent <- riskData[,3]
    }
  }
  steps <- rep(1.0,nrow(riskData));
  if (!is.null(timetoEvent)) # Adjust the censored (no events) at short times
  {
    medianTimeNocens <- median(timetoEvent[riskData[,1] == 1]); # The expected event rate
    cenevents <- riskData[,1] == 0
    # If the event time is short the is a chance of a future event in the near future. 
    steps[cenevents] <- (1.0 - exp(-timetoEvent[cenevents]/(medianTimeNocens))) #1.0 - the probability of zero events.
#    steps[steps > 1.0] <- 1.0; 
#    print(medianTimeNocens)
#    print(summary(steps))
  }
#  print(sum(riskData[,2]))
  if (uvalues < 10)
  {
    warning(paste("Not a continuous variable;",title))
    return (0)
  }
  ## Removing ties
  deltaR <- 0;
  if (uvalues < nrow(riskData))
  {
    deltaR <- IQR(riskData[,2])/nrow(riskData)
    if (deltaR == 0) 
    {
      deltaR <- sd(riskData[,2])/nrow(riskData)
    }
    dupvalues <- table(riskData[,2])
    dupvalues <- as.numeric(names(dupvalues)[dupvalues>1])
    addNoise <- riskData[,2] %in% dupvalues
    if (isProbability)
    {
        addNoise <- addNoise & (riskData[,2] < 0.999) & (riskData[,2] > 0.001)
    }
    riskData[addNoise,2] <- riskData[addNoise,2] + (runif(sum(addNoise))-0.5)*deltaR*1.0e-6

#    print(sum(addNoise))
  }
  
  isRevesed <- 1.0
  
  if (!isProbability)
  {
    if (mean(riskData[riskData[,1]==1,2]) < mean(riskData[riskData[,1]==0,2]))
    {
      riskData[,2] <- -riskData[,2]
      warning("Reversed Sign")
      isRevesed <- -1.0;
      if (!is.null(atThr))  {atThr <- -atThr;}
    }
  }
  
  riskData <- as.data.frame(riskData)
  if (is.null(rownames(riskData)))
  {
    rownames(riskData) <- as.character(c(1:nrow(riskData)))
  }
  totobs <- nrow(riskData)
  
  ## Getting the threshold values for the specified values
  if (is.null(atThr))
  {
    thr_atP <- quantile(riskData[riskData[,1]==0,2],probs=c(atRate)) - deltaR
#   print(thr_atP)
   if (!is.null(timetoEvent))
   {
     noeventdata <- cbind(riskData[riskData[,1]==0,2],steps[riskData[,1]==0])
     noeventdata <- noeventdata[order(noeventdata[,1]),]
     thrstep <- sum(noeventdata[,2])*atRate;
#     print(c(nrow(noeventdata),thrstep))
     acum <- 0;
     for (idx in c(1:nrow(noeventdata)))
     {
       acum <- acum + noeventdata[idx,2];
       if (acum < thrstep[1]) thr_atP[1] <- noeventdata[idx,1];
       if (length(thrstep) > 1)
       {
         if (acum < thrstep[2]) thr_atP[2] <- noeventdata[idx,1];
       }
     }
#     print(thr_atP)
   }
     
    if (atRate[1]<0.5)
    {
      thr_atP <- quantile(riskData[riskData[,1]==1,2],probs=c(atRate)) + deltaR
    }
    if (length(atRate)>1)
    {
      if (atRate[2]>=atRate[1])
      {
        thr_atP <- quantile(riskData[riskData[,1]==1,2],probs=c(1.0-atRate)) + deltaR
      }
      if (atRate[1] < 0.5)
      {
        thr_atP <- quantile(riskData[riskData[,1]==1,2],probs=c(atRate)) + deltaR
      }
      if (thr_atP[1] == thr_atP[2])
      {
        atRate <- atRate[1]
        thr_atP <- thr_atP[1]
      }
    }
  }
  else
  {
    thr_atP <- atThr
    atRate <- atThr
  }
  
#  print(thr_atP)

  
  
        
  numberofEvents <- sum(riskData[,1])
  numberofNoEvents <- sum((riskData[,1]==0)*steps)
  
  samplePrevalence <- numberofEvents/totobs;
  ExpectedNoEventsGain <- 1.0
  pre <- samplePrevalence
  
  if (!is.null(ExpectedPrevalence))
  {
    ExpectedNoEvents <- (1.0-ExpectedPrevalence)*totobs*samplePrevalence/ExpectedPrevalence
    ExpectedNoEventsGain <- ExpectedNoEvents/sum(riskData[,1]==0)
  }
  else
  {
    print(samplePrevalence)
    if (!is.null(riskTimeInterval) && !is.null(timetoEvent))
    {
      samplePrevalence <- sum(timetoEvent < (3*riskTimeInterval) & riskData[,1]==1)/sum(timetoEvent < (3*riskTimeInterval))
      pre <- samplePrevalence
    }
    ExpectedPrevalence <- samplePrevalence
  }
  print(c(samplePrevalence,ExpectedNoEventsGain))
  
  minRiskAtEvent <- min(riskData[riskData[,1]==1,2])
  thrsWhitinEvents <- riskData[riskData[,2] >= minRiskAtEvent,2]
  names(thrsWhitinEvents) <- rownames(riskData[riskData[,2] >= minRiskAtEvent,])
  thrsWhitinEvents <- thrsWhitinEvents[order(thrsWhitinEvents)]
  nobs <- length(thrsWhitinEvents)
  RR <- numeric(nobs)
  SEN <- numeric(nobs)
  SPE <- numeric(nobs)
  isEvent <- riskData[names(thrsWhitinEvents),1]
  PPV <- numeric(nobs)
  netBenefit <- numeric(nobs)
  URCI <- numeric(nobs)
  LRCI <- numeric(nobs)
  BACC <- numeric(nobs)
  names(RR) <- names(thrsWhitinEvents)
  names(SEN) <- names(thrsWhitinEvents)
  names(SPE) <- names(thrsWhitinEvents)
  names(isEvent) <- names(thrsWhitinEvents)
  names(PPV) <- names(thrsWhitinEvents)
  names(URCI) <- names(thrsWhitinEvents)
  names(LRCI) <- names(thrsWhitinEvents)
  names(BACC) <- names(thrsWhitinEvents)
  names(netBenefit) <- names(thrsWhitinEvents)
  
  mxthr <- max(riskData[,2])
  mithr <- min(riskData[,2])
  
  idx <- 1;
  noMoreLowEventsIdx <- minRiskAtEvent
  cenAUC <- 0;
  for (thr in thrsWhitinEvents)
  {
    
    atLowRisk <- riskData[,2] < thr
    atHighRisk <- !atLowRisk
    LowEvents <- sum(riskData[atLowRisk,1])
    if (LowEvents==0) 
    {
      LowEvents <- 0.1;
      if (noMoreLowEventsIdx == minRiskAtEvent)
      {
        noMoreLowEventsIdx <- thr
      }
    }
    HighEvents <- sum(riskData[atHighRisk,1]);
    SEN[idx] <-  HighEvents/numberofEvents;
    SPE[idx] <- sum((riskData[atLowRisk,1]==0)*steps[atLowRisk])/numberofNoEvents;
    BACC[idx] <- (SEN[idx] + SPE[idx])/2
    n1 <- sum(atHighRisk*steps)
    n2 <- sum(atLowRisk*steps)
    PPV[idx] <- HighEvents/n1
    LowFraction <- LowEvents/n2;
    HighFraction <- HighEvents/n1;
    x1 <- HighEvents
    x2 <- LowEvents
    if ((n2/numberofNoEvents) > 0.01)
    {
      RR[idx] <- HighFraction/LowFraction
      cidelta <- 1.96*sqrt((n1-x1)/(x1*n1)+(n2-x2)/(x2*n2))
      URCI[idx] <- exp(log(RR[idx])+cidelta)
      LRCI[idx] <- exp(log(RR[idx])-cidelta)
    }
    else
    {
      if (idx>1)
      {
        RR[idx] <- RR[idx-1]
      }
      else
      {
        RR[idx] <- 1.0;
      }
    }
    if (isProbability)
    {
      netBenefit[idx] <- (HighEvents - ExpectedNoEventsGain*sum(riskData[atHighRisk,1]==0)*thrsWhitinEvents[idx]/(1.0001-thrsWhitinEvents[idx]))/totobs;
    }
    if (idx>1)
    {
      cenAUC <- cenAUC + 0.5*(SEN[idx]+SEN[idx-1])*(SPE[idx]-SPE[idx-1])
    }
#    cat(SEN[idx],",",SPE[idx],"\n");
    idx <- idx + 1;
  }
  idx <- idx - 1;
  cenAUC <- cenAUC + 0.5*(1.0+SEN[1])*SPE[1] + 0.5*SEN[idx]*(1.0-SPE[idx]);
  names(cenAUC) <- NULL
 # print(cenAUC)
  colors <- heat.colors(10)
  rgbcolors <- rainbow(10)

  tmop <- par(no.readonly = TRUE)
  

  
  if (isProbability)
  {
    ## Observed vs Cumulative plot
#    par(pty='s',mfrow=c(1,2),cex=0.5)
    par(pty='s')
    xdata <- riskData[order(-riskData[,2]),]
    tmpObserved <- numeric(nrow(xdata))
    tmpExplected <- numeric(nrow(xdata))
    tmpObserved[1] <- xdata[1,1]
    tmpExplected[1] <- xdata[1,2]
    pgzero <- xdata[,2]
    

    for (idx in c(2:nrow(xdata)))
    {
      if (xdata[idx,1]==1)
      {
        tmpExplected[idx] <- tmpExplected[idx-1] + pgzero[idx];
      }
      else
      {
        tmpExplected[idx] <- tmpExplected[idx-1] + pgzero[idx]*ExpectedNoEventsGain;
      }
      tmpObserved[idx] <- tmpObserved[idx-1] + xdata[idx,1];
    }
    totObserved <- sum(xdata[,1])
    tokeep <- (tmpExplected > 0.10*totObserved)
    tokeep <- tokeep & (tmpObserved > 0.10*totObserved)
    CumulativeOvs <- as.data.frame(cbind(Observed=tmpObserved,Cumulative=tmpExplected,Included=tokeep))
    
    OEratio <- tmpObserved[tokeep]/tmpExplected[tokeep]
    OAcum95ci <- c(mean=mean(OEratio),metric95ci(OEratio))
    
    observedCI <- stats::poisson.test(max(tmpObserved[tokeep]),max(tmpExplected[tokeep]), conf.level = 0.95 )
    OARatio <- observedCI
    OARatio$estimate <- c(OARatio$estimate,OARatio$conf.int,OARatio$p.value)
    names(OARatio$estimate) <- c("O/A","Low","Upper","p.value")
    
    
    rownames(CumulativeOvs) <- rownames(xdata)
    maxobs <- max(c(tmpObserved,tmpExplected))
    
    ## Probability Calibration Curve
    if (plotRR)
    {
     
      plot(tmpExplected,tmpObserved,
           ylab="Observed",
           xlab="Cumulative Probability",
           xlim=c(0,maxobs),
           ylim=c(0,maxobs),
           main=paste("Cumulative vs. Observed:",title),
           pch=c(5+14*xdata[,1]),
           col=c(1+xdata[,1]),
           cex=c(0.2+xdata[,1]))
           se <- 2*sqrt(tmpObserved)
      errbar(tmpExplected,tmpObserved,tmpObserved-se,tmpObserved+se,add=TRUE,pch=0,errbar.col="gray",cex=0.25)
      eventExpected <- tmpExplected[xdata[,1]==1]
      eventObserved <- tmpObserved[xdata[,1]==1]
      coltimes <- rep(1,length(eventExpected));
      if (!is.null(timetoEvent))
      {
        ncolors <- 7
        timecolors  <- heat.colors(6);
        coltimes <- timetoEvent[order(-riskData[,2])];
        coltimes <- coltimes[xdata[,1]==1]
        mintime <- min(coltimes);
        maxtime <- max(coltimes);
        thetimes <- ncolors - 1 - floor((ncolors-1)*(coltimes-mintime)/(maxtime-mintime));
        coltimes <- timecolors[1+thetimes];
        thetimes <- c(0:(ncolors-1))*(maxtime-mintime)/ncolors + mintime;
        legtxt <- sprintf("%3.1f",thetimes)
        text(0.10*maxobs,0.99*maxobs,"Time to Event",cex=0.65)
        corrplot::colorlegend(timecolors, legtxt,xlim=c(0.025*maxobs,0.125*maxobs),ylim=c(c(0.65*maxobs,0.95*maxobs)),cex=0.65)

      }
      lines(x=c(0,maxobs),y=c(0,maxobs),lty=2)
      points(eventExpected,eventObserved,
              pch=c(19),
                 col=coltimes
                 )

      legend("bottomright",legend=c("Observed","Expected"),
             lty=c(-1,2),
             pch=c(19,-1),
             col=c(2,1),cex=0.75)
      par(tmop)
    }

    
    
    ## Decision curve analysis
    pshape <- 4 + 12*isEvent
    xmax <- min(quantile(thrsWhitinEvents,probs=c(0.95),0.95))
    xmin <- min(quantile(thrsWhitinEvents,probs=c(0.01),0.001))
    ymin <- min(quantile(netBenefit,probs=c(0.10)),0)
    DCA <- as.data.frame(cbind(Thrs=thrsWhitinEvents,NetBenefit=netBenefit))
    rownames(DCA) <- names(thrsWhitinEvents)
    
## Decision Curve Analysis
    if (plotRR)
    {
      plot(thrsWhitinEvents,netBenefit,main=paste("Decision Curve Analysis:",title),ylab="Net Benefit",xlab="Threshold",
         ylim=c(ymin,pre),
         xlim=c(xmin,xmax),
         pch=pshape,col=rgbcolors[1+floor(10*(1.0-SEN))],cex=(0.35 + isEvent))
      fiveper <- as.integer(0.05*length(netBenefit)+0.5)
      range <- c(fiveper:length(netBenefit)-fiveper)
      lfit <-try(loess(netBenefit[range]~thrsWhitinEvents[range],span=0.5));
      if (!inherits(lfit,"try-error"))
      {
        plx <- try(predict(lfit,se=TRUE))
        if (!inherits(plx,"try-error"))
        {
          lines(thrsWhitinEvents[range],plx$fit,lty=1)
          lines(thrsWhitinEvents[range],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
          lines(thrsWhitinEvents[range],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
        }
      }
      abline(h=0,col="blue")
      lines(x=c(0,ExpectedPrevalence),y=c(pre,0),col="red")
      legtxt <- sprintf("%3.1f",c(5:0)/5)
      corrplot::colorlegend(rgbcolors, legtxt,xlim=c(0.92*(xmax-xmin)+xmin,0.98*(xmax-xmin)+xmin),ylim=c(pre*0.45,pre*0.1),cex=0.75)
      text(0.92*xmax,pre*0.5,"SEN")
      
      thrlty <- c(1)
      thrLey <- "H. Thr"
      abline(v=thr_atP[1],col="gray",lty=thrlty)
      if (length(thr_atP)>1) 
      {
        abline(v=thr_atP[2],col="gray",lty=2)
        thrlty <- c(thrlty,2)
        thrLey <- c(thrLey,"L. Thr")
      }
      
      legend("topright",legend=c("No Event","Event","loess fit",thrLey),
             pch=c(4,16,-1,-1,-1),
             col=c(1,1,1,"gray","gray"),
             lty=c(-1,-1,1,thrlty)
             )
      legend("bottomleft",legend=c("Treat All","Treat None"),
             lty=c(1,1),
             col=c("red","blue"),cex=0.5)
      par(tmop)
    }
    
  }
  ## Relative Risk plot
  
  ymax <- quantile(RR,probs=c(0.99))
  pshape <- 2 + 14*isEvent
  RRData <- as.data.frame(cbind(thr=thrsWhitinEvents,
                                Sensitivity=SEN,
                                Specificity=SPE,
                                BACC  = BACC,
                                PPV=PPV,
                                RR=RR,
                                LRR=LRCI,
                                URR=URCI,
                                isEvent=isEvent))
  rownames(RRData) <- names(thrsWhitinEvents)
  lfit <- try(loess(RR~SEN,span=0.5));
#  print(thr_atP)
  ## Get the location of the thresholds
  maxBACC <- max(BACC)
  RRval <- RR[thrsWhitinEvents <= thr_atP[1] & SPE > 0.10]
  thrsval <- thrsWhitinEvents[thrsWhitinEvents <= thr_atP[1] & SPE > 0.10]
  maxRR <- max(RRval)
  thr_values <- c(thr_atP,
                  at_max_BACC=min(thrsWhitinEvents[BACC==maxBACC]),
                  at_max_RR=min(thrsval[RRval==maxRR]),
                  atSENE100=noMoreLowEventsIdx)
  if (isProbability)
  {
    thr_values <- c(thr_values,at_0.5=0.5)
  }
  
  thrLoc <- rep(1,length(thr_values))
#  print(thr_values)
  for (tthr in c(1:length(thr_values)))
  {
    thrLoc[tthr] <- which.min(abs(thrsWhitinEvents-thr_values[tthr]))
  }
#  print(thrLoc)
  thrPoints <- as.data.frame(cbind(isRevesed*thrsWhitinEvents[thrLoc],RR[thrLoc],LRCI[thrLoc],URCI[thrLoc],SEN[thrLoc],SPE[thrLoc],BACC[thrLoc]))
  colnames(thrPoints) <- c("Thr","RR","RR_LCI","RR_UCI","SEN","SPE","BACC")
  if (isProbability)
  {
    thrPoints$NetBenefit <- netBenefit[thrLoc]
    rownames(thrPoints) <- c(paste("@",round(atRate,3),sep=":"),"@MAX_BACC","@MAX_RR","@SEN100","p(0.5)")
  }
  else
  {
    rownames(thrPoints) <- c(paste("@",round(atRate,3),sep=":"),"@MAX_BACC","@MAX_RR","@SPE100")
  }
#  print(thrPoints)
  
  ## Sensitivity, Specificity and RR at top threshold
  lowRisk <- (riskData[,2] < thr_atP[1])
  LowEventsFrac <- sum(riskData[lowRisk,1]*steps[lowRisk])/sum(lowRisk*steps)
  HighEventsFrac <- sum(riskData[!lowRisk,1]*steps[!lowRisk])/sum(!lowRisk*steps)
  sensitivity=sum(riskData[!lowRisk,1])/numberofEvents
  who <- names(SEN)[SEN==sensitivity]
  RRAtSen <- c(est=mean(RR[who]),lower=mean(LRCI[who]),upper=mean(URCI[who]))
  

  
  specificity=sum((riskData[lowRisk,1]==0)*steps[lowRisk])/numberofNoEvents
  
  ## Risk Ratio Plot
  if (plotRR)
  {
    ypmax <- max(c(quantile(URCI,probs=c(0.75),na.rm = TRUE),ymax))

    par(mfrow=c(1,1))
  
    plot(SEN,RR,cex=(0.35 + PPV),
       pch=pshape,
       col=colors[1+floor(10*(1.0-SPE))],
       xlim=c(0,1.15),
#       ylim=c(1.0,ymax+0.5),
       ylim=c(0.0,ypmax*1.2),
#       log="y",
       main=paste("Relative Risk:",title))
    atevent <- isEvent==1
    errbar(SEN[atevent],RR[atevent],URCI[atevent],LRCI[atevent],add=TRUE,pch=0,errbar.col="gray",cex=0.25)
    abline(h=1,col="red")
    points(SEN,RR,cex=(0.35 + PPV),
         pch=pshape,
         col=colors[1+floor(10*(1.0-SPE))],
    )  
#    lfit <- try(loess(RR~SEN,span=0.5));
    if (!inherits(lfit,"try-error"))
    {
      plx <- try(predict(lfit,se=TRUE))
      if (!inherits(plx,"try-error"))
      {
        lines(SEN,plx$fit,lty=1)
        lines(SEN,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
        lines(SEN,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
      }
    }
  
    legtxt <- sprintf("%3.1f",c(5:0)/5)
    corrplot::colorlegend(heat.colors(10), legtxt,xlim=c(1.05,1.175),ylim=c((ypmax-0.75)*0.45+1.0,ypmax/10+1.0),cex=0.75)
    text(1.075,1,"SPE")
    sizp <- c(5:1)/5.0 + 0.35
    legend("topright",legend=legtxt[1:5],pch=16,cex=sizp)
    legend("bottomleft",legend=c("No","Yes","Loess Fit"),pch=c(2,16,-1),lty=c(0,0,1),cex=0.75)
    text(0.95,ypmax+0.5,"PPV->")


    abline(v=sensitivity,col="blue")
    if (length(thr_atP)>1) 
    {
       lowRisk2 <- (riskData[,2] < thr_atP[2])
       sensitivity2=sum(riskData[!lowRisk2,1])/numberofEvents
       abline(v=sensitivity2,col="cyan",lty=2)
    }

    text(x=sensitivity,y=ypmax*1.2,sprintf("Index(%3.2f)=%4.3f",specificity,isRevesed*thr_atP[1]),pos=4 - 2*(sensitivity>0.5) ,cex=0.7)
    text(x=sensitivity,y=0.9*(ypmax*1.2-1.0)+1.0,
         sprintf("RR(%3.2f)=%4.3f",
                 sensitivity,
                 RRAtSen[1]),
         pos=4 - 2*(sensitivity>0.5),cex=0.7)
    par(tmop)
  }
  
  ## ROC At threshold
  
  ROCAnalysis <- predictionStats_binary(riskData[,c(1,2)],
                                          thr=thr_atP[1])
                                          
  ## ROC Plot
  if (plotRR)
  {
    par(tmop)
    par(pty='s');
    procdta <- roc(riskData[,1],riskData[,2],
                              grid=c(0.1, 0.1),
                              grid.col=c("gray", "gray"),
                              print.auc=TRUE,
                               quiet = TRUE,
                              main=paste("ROC:",title),
                               plot=TRUE,
                              col="black",
                              lty=1,
                              lwd=3,
                              ) 
    lines(x=SPE,y=SEN,lty=2,lwd=1,col="red");
    text(0.4,0.4,sprintf("AAUC:%4.3f",cenAUC),col="red");
    legend("bottomright",legend=c("Observed","p-Adjusted"),pch=c(-1,-1),lty=c(1,2),col=c("black","red"),cex=0.80)

    par(tmop)
  }
  surfit <- NULL
  surdif <- NULL
  LogRankE <- NULL
  cstat <- NULL
  timetoEventData <- NULL
#  print(sum(riskData[,2]))
  if (!is.null(timetoEvent))
  {
      timetoEventData <- as.data.frame(cbind(eStatus=riskData[,1],
                           class=1*(riskData[,2]>=thr_atP[1]),
                           eTime=timetoEvent,
                           risk=riskData[,2])
                           )
      if (length(thr_atP)>1)
      {
        timetoEventData$class <- 1*(riskData[,2]>=thr_atP[1]) + 1*(riskData[,2]>=thr_atP[2])
      }
      paletteplot <- c("green", "red")
      if (isRevesed > 0)
      {
        labelsplot <- c("Low Risk",sprintf("High Risk > %5.3f",thr_atP[1]));
        if (length(thr_atP)>1)
        {
          labelsplot <- c(sprintf("Low Risk < %5.3f",thr_atP[2]),sprintf("%5.3f <= Risk < %5.3f",thr_atP[2],thr_atP[1]),
                          sprintf("High Risk >= %5.3f",thr_atP[1]));
          paletteplot <- c("green", "blue","red")
        }
      }
      else
      {
        labelsplot <- c("Low Risk",sprintf("High Risk < %5.3f",-thr_atP[1]));
        if (length(thr_atP)>1)
        {
          labelsplot <- c(sprintf("Low Risk > %5.3f",-thr_atP[2]),sprintf("%5.3f >= Risk > %5.3f",-thr_atP[2],-thr_atP[1]),
                          sprintf("High Risk <= %5.3f",-thr_atP[1]));
          paletteplot <- c("green", "blue","red")
        }
        
      }
      if (isProbability)
      {
        ## Time Plot
        prisk <- riskData[,2]
        timetoEventData$lammda <- expectedEventsPerInterval(prisk)
#        print(sum(prisk))
        TOEratio <- tmpObserved/tmpExplected

        aliveEvents <- timetoEventData
        aliveEvents <- aliveEvents[order(-aliveEvents$risk),]
        aliveEvents <- aliveEvents[order(aliveEvents$eTime),]
        
        if (!is.null(riskTimeInterval))
        {
          timeInterval <- riskTimeInterval
        }
        else
        {
          timeInterval <- mean(subset(aliveEvents,
                                aliveEvents$eStatus==1)$eTime);
        }
        timetoEventData$expectedTime <- meanTimeToEvent(prisk,timeInterval)
        attr(timetoEventData,"ClassNames") <- labelsplot
        timed <- unique(aliveEvents[aliveEvents$eStatus==1,"eTime"])
        maxtime <- max(timed)
        Observed <- c(1:length(timed))
        Expected <- Observed
        cClass <- Observed
        aliveEvents[aliveEvents$eStatus == 0,"lammda"] <- aliveEvents[aliveEvents$eStatus == 0,"lammda"]*ExpectedNoEventsGain
#        print(mean(aliveEvents$lammda));
        lastObs <- 0;
        allTimes <- aliveEvents$eTime
        lasttime <- 0
#        print(c(nrow(aliveEvents),sum(aliveEvents$risk),timeInterval))
      
        minTimeCens <- min(aliveEvents[aliveEvents$eStatus==0,"eTime"])
        meanTimeCens <- mean(aliveEvents[aliveEvents$eStatus==0,"eTime"])
        meanCenlammda <- mean(aliveEvents[(aliveEvents$eStatus==0),"lammda"])
        meanEvelammda <- mean(aliveEvents[(aliveEvents$eStatus==1),"lammda"])
        pzeroAtMean <- exp(-meanCenlammda*meanTimeCens/timeInterval);
        pzeroAtMin <- exp(-meanCenlammda*minTimeCens/timeInterval);
        pzEventAtMin <- exp(-meanEvelammda*minTimeCens/timeInterval);
        ExpectCenEventsAtMin <- sum(aliveEvents$eStatus==0)*(1.0 - pzeroAtMin);
        ExpectEventsAtMin <- sum(aliveEvents$eStatus==1)*(1.0 - pzEventAtMin);
        obsAtMin <- sum((aliveEvents$eStatus==1) & (aliveEvents$eTime < minTimeCens))
        fractNoObserb <- (ExpectCenEventsAtMin + ExpectEventsAtMin - obsAtMin)/(ExpectCenEventsAtMin+ExpectEventsAtMin);
        if (fractNoObserb < 0) ## Something is ODD in estimated probabilities
        {
          fractNoObserb <- 0;
        }
        if ((mean(TOEratio) > 1.5) | (mean(TOEratio) < 0.75)) ## Not calibrated probabilities
        {
          fractNoObserb <- 0;
        }
        print(c(mean(TOEratio),pzeroAtMin,minTimeCens,meanTimeCens,ExpectCenEventsAtMin,ExpectEventsAtMin,obsAtMin,fractNoObserb));
        passAcum <- 0;
        acuexpecCen <- 0;
        for (idx in c(1:length(timed)))
        {
          whosum <- (allTimes <= timed[idx]) & (allTimes > lasttime)
          deltatime <- (timed[idx] - lasttime)
          Observed[idx] <- lastObs + sum(aliveEvents[whosum,"eStatus"])
          cClass[idx] <- round(mean(aliveEvents[whosum,"class"]),0)
          lastObs <- Observed[idx]
          eevents <- sum(aliveEvents[allTimes > lasttime,"lammda"])*deltatime/timeInterval;
          Expected[idx] <- passAcum + eevents;
          if ((timed[idx] <= minTimeCens) || (idx==1))
          {
            Noobseevents <- sum(aliveEvents[(allTimes > lasttime) & aliveEvents$eStatus==0,"lammda"])*deltatime/timeInterval;
            missevents <- Noobseevents*fractNoObserb
            Expected[idx] <- Expected[idx] - missevents;
#            print(c(timed[idx],Noobseevents,missevents));
          }
          if (Expected[idx] < 1)
          {
            Expected[idx] <- 1
          }
          passAcum <- Expected[idx]
          lasttime <- timed[idx]
        }
        
        totObserved <- max(Observed)
        maxevents <- max(c(Observed,Expected))
        tokeep <- (Expected > 0.10*totObserved)
        tokeep <- tokeep & (Observed > 0.10*totObserved)
        OEData <- as.data.frame(cbind(time=timed,Observed=Observed,Expected=Expected,Included=tokeep,class=cClass))
        OEratio <- Observed[tokeep]/Expected[tokeep]
        OE95ci <- c(mean=mean(OEratio),metric95ci(OEratio))

        observedCI <- stats::poisson.test(max(Observed[tokeep]),max(Expected[tokeep]), conf.level = 0.95 )
        OERatio <- observedCI
        OERatio$estimate <- c(OERatio$estimate,OERatio$conf.int,OERatio$p.value)
        ## Estimation of O/E for the threshold values
        names(OERatio$estimate) <- c("O/E","Low","Upper","p.value")
        procat <- timetoEventData$class
        values <- unique(procat)
        values <- values[order(values)]
        names(values) <- c("low",names(thr_atP))
        obs <- numeric()
        tmpExplected <- numeric()
        lci <- numeric()
        uci <- numeric()
        lciOE <- numeric()
        uciOE <- numeric()
        pval <- numeric()
        hazards <- timetoEventData$lammda*timetoEventData$eTime/timeInterval
        for (ct in values)
        {
          totobs <- sum(riskData[procat==ct,1])
          obs <-c(obs,totobs)
          expe <- sum(hazards[procat==ct]);
          tmpExplected <- c(tmpExplected,expe)
          pt <- stats::poisson.test(totobs,1)
          lci <- c(lci,pt$conf.int[1])
          uci <- c(uci,pt$conf.int[2])
          pt2 <- stats::poisson.test(totobs,expe)
          lciOE <- c(lciOE,pt2$conf.int[1])
          uciOE <- c(uciOE,pt2$conf.int[2])
          pval <- c(pval,pt2$p.value)
        }
        totobservedCI <- stats::poisson.test(max(Observed),1, conf.level = 0.95 )
        totalEstimated <- c(max(Observed),totobservedCI$conf.int,max(Expected),max(Observed)/max(Expected),OERatio$conf.int,OERatio$p.value)
        OERatio$atThrEstimates <- as.data.frame(cbind(obs,lci,uci,tmpExplected,obs/tmpExplected,lciOE,uciOE,pval))
        OERatio$atThrEstimates <- rbind(totalEstimated,OERatio$atThrEstimates)
        colnames(OERatio$atThrEstimates) <- c("Observed","L.CI","H.CI","Expected","O/E","Low","Upper","pvalue")
        rownames(OERatio$atThrEstimates) <- c("Total",names(values))
        ## Now lets plot
        colorsAtMin <- 1 + (timed < minTimeCens);
        leyGMin <- sprintf("Time > %4.2f",minTimeCens)
        leyLMin <- sprintf("Time < %4.2f",minTimeCens)
        if (plotRR)
        {
          orderplot <- order(cClass)
          pchtypes <- c(1,16,17)
          par(mfrow=c(1,1))
          plot(timed,Expected,
#                pch=4,
                type="b",cex=0.5,
             main=paste("Time vs. Events:",title),
             ylab="Events",
             xlab="Time",
#             col=colorsAtMin,
             col=1,
             pch= 3 + colorsAtMin,
             ylim=c(0,1.05*maxevents),
             xlim=c(0,1.05*maxtime),
             )
          se <- 2*sqrt(Observed)
          errbar(timed,Observed,Observed-se,Observed+se,add=TRUE,pch=0,errbar.col="gray",cex=0.25)
          points(timed,Expected,pch=3+colorsAtMin,type="p",cex=0.5,col=colorsAtMin)
          points(timed[orderplot],Observed[orderplot],pch=pchtypes[1+cClass[orderplot]],col=paletteplot[1+cClass[orderplot]])
          legend("topleft",legend=c(leyGMin,leyLMin,labelsplot),pch=c(4,5,pchtypes),lty=c(1,1,0,0,0),col=c(1,2,paletteplot),cex=0.80)
        }
      }
        ## Survival plot
      LogRankE <- EmpiricalSurvDiff(times=timetoEventData$eTime,
                  status=timetoEventData$eStatus,
                  groups=1*(timetoEventData$class > 0),
                  plots=FALSE,main=paste("Kaplan-Meier:",title))
      
      surfit <- survival::survfit(survival::Surv(eTime,eStatus)~class,data = timetoEventData)
      surdif <- survival::survdiff(survival::Surv(eTime,eStatus)~class,data = timetoEventData)
      cstat <- rcorr.cens(-timetoEventData$risk,survival::Surv(timetoEventData$eTime,timetoEventData$eStatus))
      cstatCI <- c(mean=cstat[1],concordance95ci(as.data.frame(cbind(times=timetoEventData$eTime,
                                       status=timetoEventData$eStatus,
                                       preds=-timetoEventData$risk))))
      cstat <- as.list(cstat)
      cstat$cstatCI <- cstatCI
      
      if (plotRR)
      {
        par(mfrow=c(1,1))
        graph <- survminer::ggsurvplot(surfit,
                                       data = timetoEventData, 
                                       conf.int = TRUE, 
                                       legend.labs = labelsplot,
                                       palette = paletteplot,
                                       ylim = ysurvlim,
                                       ggtheme = ggplot2::theme_bw() + 
                                         ggplot2::theme(plot.title = ggplot2::element_text(
                                           hjust = 0.5, face = "bold", size = 16)),
                        title = paste("Kaplan-Meier:",title),
                        risk.table = TRUE,
                        tables.height = 0.2,
                        tables.theme = survminer::theme_cleantable())
        print(graph)
        par(tmop)
        
      }
    

  }
  
  result <- list(CumulativeOvs=CumulativeOvs,
                 OEData=OEData,
                 DCA=DCA,
                 RRData=RRData,
                 timetoEventData=timetoEventData,
                 keyPoints=thrPoints,
                 OERatio=OERatio,
                 OE95ci=OE95ci,
                 OARatio=OARatio,
                 OAcum95ci=OAcum95ci,
                 fit=lfit,
                 ROCAnalysis=ROCAnalysis,
                 prevalence=pre,
                 thr_atP= isRevesed*thr_atP,
                 c.index=cstat,
                 cenAUC=cenAUC,
                 surfit=surfit,
                 surdif=surdif,
                 LogRankE=LogRankE
                 )
  return (result)
}
