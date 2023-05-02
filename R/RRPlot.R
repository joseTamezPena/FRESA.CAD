RRPlot <-function(riskData=NULL,
                  timetoEvent=NULL,
                  riskTimeInterval=NULL,
                  ExpectedPrevalence=NULL,
                  atProb=c(0.90,0.80),
                  atThr=NULL,
                  plotRR=TRUE,
                  title="",
                  ysurvlim=c(0,1.0))
{
  OERatio <- NULL
  OE95ci <- NULL
  OAcum95ci <- NULL
  CumulativeOvs <- NULL
  DCA <- NULL
  OEData <- NULL
  
  uvalues <- length(unique(riskData[,2]))
  
  if (uvalues < 10)
  {
    warning(paste("Not a continous variable;",title))
    return (0)
  }
  ## Removing ties 
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
    riskData[addNoise,2] <- riskData[addNoise,2] + (runif(sum(addNoise))-0.5)*deltaR*1.0e-6
#    print(sum(addNoise))
  }
  
  if (mean(riskData[riskData[,1]==1,2]) < mean(riskData[riskData[,1]==0,2]))
  {
    riskData[,2] <- -riskData[,2]
    warning("Reversed Sign")
  }
  
  riskData <- as.data.frame(riskData)
  if (is.null(rownames(riskData)))
  {
    rownames(riskData) <- as.character(c(1:nrow(riskData)))
  }
  totobs <- nrow(riskData)
if (!requireNamespace("corrplot", quietly = TRUE)) {
		install.packages("corrplot", dependencies = TRUE)
		}
        
  numberofEvents <- sum(riskData[,1])
  numberofNoEvents <- sum(riskData[,1]==0)
  
  samplePrevalence <- numberofEvents/totobs;
  ExpectedNoEventsGain <- 1.0
  pre <- samplePrevalence
  
  if (!is.null(ExpectedPrevalence))
  {
    ExpectedNoEvents <- (1.0-ExpectedPrevalence)*totobs*samplePrevalence/ExpectedPrevalence
    ExpectedNoEventsGain <- ExpectedNoEvents/numberofNoEvents
  }
  else
  {
    ExpectedPrevalence <- samplePrevalence
  }
  
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
  isProbability <- (mithr>=0) && (mxthr<=1.0)
  
  idx <- 1;
  noMoreLowEventsIdx <- minRiskAtEvent
  for (thr in thrsWhitinEvents)
  {
    
    atLowRisk <- riskData[,2] < thr
    atHighRisk <- !atLowRisk
    LowEvents <- sum(riskData[atLowRisk,1])
    if (LowEvents==0) 
    {
      LowEvents <- 0.1;
      noMoreLowEventsIdx <- thr
    }
    HighEvents <- sum(riskData[atHighRisk,1]);
    SEN[idx] <-  HighEvents/numberofEvents;
    SPE[idx] <- sum(riskData[atLowRisk,1]==0)/numberofNoEvents;
    BACC[idx] <- (SEN[idx] + SPE[idx])/2
    n1 <- sum(atHighRisk)
    n2 <- sum(atLowRisk)
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
    netBenefit[idx] <- (HighEvents - ExpectedNoEventsGain*sum(riskData[atHighRisk,1]==0)*thrsWhitinEvents[idx]/(1.0001-thrsWhitinEvents[idx]))/totobs;
    idx <- idx + 1;
  }
  colors <- heat.colors(10)
  rgbcolors <- rainbow(10)

  tmop <- par(no.readonly = TRUE)
  
  if (isProbability)
  {
    ## Observed vs Cumulative plot
#    par(pty='s',mfrow=c(1,2),cex=0.5)
    par(pty='s')
    xdata <- riskData[order(-riskData[,2]),]
    observed <- numeric(nrow(xdata))
    expected <- numeric(nrow(xdata))
    observed[1] <- xdata[1,1]
    expected[1] <- xdata[1,2]
    for (idx in c(2:nrow(xdata)))
    {
      if (xdata[idx,1]==1)
      {
        expected[idx] <- expected[idx-1]+xdata[idx,2];
      }
      else
      {
        expected[idx] <- expected[idx-1]+xdata[idx,2]*ExpectedNoEventsGain;
      }
      observed[idx] <- observed[idx-1]+xdata[idx,1];
    }
    totObserved <- sum(xdata[,1])
    tokeep <- (expected > 0.10*totObserved)
    tokeep <- tokeep & (observed > 0.10*totObserved)
    CumulativeOvs <- as.data.frame(cbind(Observed=observed,Cumulative=expected,Included=tokeep))
    
    OEratio <- observed[tokeep]/expected[tokeep]
    OAcum95ci <- c(mean=mean(OEratio),metric95ci(OEratio))
    
    observedCI <- stats::poisson.test(max(observed[tokeep]),max(expected[tokeep]), conf.level = 0.95 )
    OARatio <- observedCI
    OARatio$estimate <- c(OARatio$estimate,OARatio$conf.int,OARatio$p.value)
    names(OARatio$estimate) <- c("O/A","Low","Upper","p.value")
    
    
    rownames(CumulativeOvs) <- rownames(xdata)
    maxobs <- max(c(observed,expected))
    
    if (plotRR)
    {
      plot(expected,observed,
           ylab="Observed",
           xlab="Cumulative Probability",
           xlim=c(0,maxobs),
           ylim=c(0,maxobs),
           main=paste("Cumulative vs. Observed:",title),
           pch=c(5+14*xdata[,1]),
           col=c(1+xdata[,1]),
           cex=c(0.2+xdata[,1]))
      lines(x=c(0,maxobs),y=c(0,maxobs),lty=2)
      legend("bottomright",legend=c("Event","Expected"),
             lty=c(-1,2),
             pch=c(19,-1),
             col=c(2,1),cex=0.75)
      par(tmop)
    }

    
    
    ## Decision curve analysis
    pshape <- 4 + 12*isEvent
    xmax <- min(quantile(thrsWhitinEvents,probs=c(0.95),0.95))
    ymin <- min(quantile(netBenefit,probs=c(0.10)),0)
    DCA <- as.data.frame(cbind(Thrs=thrsWhitinEvents,NetBenefit=netBenefit))
    rownames(DCA) <- names(thrsWhitinEvents)
    
    if (plotRR)
    {
      plot(thrsWhitinEvents,netBenefit,main=paste("Decision Curve Analysis:",title),ylab="Net Benefit",xlab="Threshold",
         ylim=c(ymin,pre),
         xlim=c(0,xmax),
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
      corrplot::colorlegend(rgbcolors, legtxt,xlim=c(0.92*xmax,0.98*xmax),ylim=c(pre*0.45,pre*0.1),cex=0.75)
      text(0.92*xmax,pre*0.5,"SEN")
      
      legend("topright",legend=c("No Event","Event","loess fit"),
             pch=c(4,16,-1),
             col=c(1,1,1),
             lty=c(-1,-1,1)
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
  if (is.null(atThr))
  {
    thr_atP <- quantile(riskData[riskData[,1]==0,2],probs=c(atProb))
    if (length(atProb)>1)
    {
      if (atProb[2]>atProb[1])
      {
        thr_atP <- quantile(riskData[riskData[,1]==1,2],probs=c(1.0-atProb))
      }
    }
  }
  else
  {
    thr_atP <- atThr
    atProb <- atThr
  }
  ## Get the location of the thresholds
  maxBACC <- max(BACC)
  RRval <- RR[thrsWhitinEvents <= thr_atP[1] & SPE > 0.10]
  thrsval <- thrsWhitinEvents[thrsWhitinEvents <= thr_atP[1] & SPE > 0.10]
  maxRR <- max(RRval)
  thr_values <- c(thr_atP,
                  at_max_BACC=min(thrsWhitinEvents[BACC==maxBACC]),
                  at_max_RR=min(thrsval[RRval==maxRR]),
                  atSPE100=noMoreLowEventsIdx)
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
  thrPoints <- as.data.frame(cbind(thrsWhitinEvents[thrLoc],RR[thrLoc],SEN[thrLoc],SPE[thrLoc],BACC[thrLoc]))
  colnames(thrPoints) <- c("Thr","RR","SEN","SPE","BACC")
  if (isProbability)
  {
    rownames(thrPoints) <- c(paste("@",atProb,sep=":"),"@MAX_BACC","@MAX_RR","@SPE100","p(0.5)")
  }
  else
  {
    rownames(thrPoints) <- c(paste("@",atProb,sep=":"),"@MAX_BACC","@MAX_RR","@SPE100")
  }
#  print(thrPoints)
  
  ## Sensitivity, Specificity and RR at top threshold
  lowRisk <- riskData[,2] < thr_atP[1]
  LowEventsFrac <- sum(riskData[lowRisk,1])/sum(lowRisk)
  HighEventsFrac <- sum(riskData[!lowRisk,1])/sum(!lowRisk)
  sensitivity=sum(riskData[!lowRisk,1])/numberofEvents
  who <- names(SEN)[SEN==sensitivity]
  RRAtSen <- c(est=mean(RR[who]),lower=mean(LRCI[who]),upper=mean(URCI[who]))
  
  specificity=sum(riskData[lowRisk,1]==0)/numberofNoEvents
  
  if (plotRR)
  {
    par(mfrow=c(1,1))
  
    plot(SEN,RR,cex=(0.35 + PPV),
       pch=pshape,
       col=colors[1+floor(10*(1.0-SPE))],
       xlim=c(0,1.15),
       ylim=c(1.0,ymax+0.5),
#       log="y",
       main=paste("Relative Risk:",title))
    atevent <- isEvent==1
    errbar(SEN[atevent],RR[atevent],URCI[atevent],LRCI[atevent],add=TRUE,pch=0,errbar.col="gray",cex=0.25)
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
    corrplot::colorlegend(heat.colors(10), legtxt,xlim=c(1.05,1.175),ylim=c((ymax-0.75)*0.45+1.0,ymax/10+1.0),cex=0.75)
    text(1.075,1,"SPE")
    sizp <- c(5:1)/5.0 + 0.35
    legend("topright",legend=legtxt[1:5],pch=16,cex=sizp)
    legend("bottomleft",legend=c("No","Yes","Loess Fit"),pch=c(2,16,-1),lty=c(0,0,1),cex=0.75)
    text(0.95,ymax+0.5,"PPV->")


    abline(v=sensitivity,col="blue")
    text(x=sensitivity,y=ymax,sprintf("Index(%3.2f)=%4.3f",specificity,thr_atP[1]),pos=4 - 2*(sensitivity>0.5) ,cex=0.7)
    text(x=sensitivity,y=0.9*(ymax-1.0)+1.0,
         sprintf("RR(%3.2f)=%4.3f",
                 sensitivity,
                 RRAtSen[1]),
         pos=4 - 2*(sensitivity>0.5),cex=0.7)
    par(tmop)
  }
  
  ## ROC At threshold
  
  if (plotRR)
  {
    par(tmop)
    ROCAnalysis <- predictionStats_binary(riskData,
                                        plotname=paste("ROC:",title),
                                        thr=thr_atP[1],cex=0.85)

    par(tmop)
  }
  else
  {
    ROCAnalysis <- predictionStats_binary(riskData,
                                          thr=thr_atP[1])
  }
  surfit <- NULL
  surdif <- NULL
  LogRankE <- NULL
  cstat <- NULL

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
      if (isProbability)
      {
        ## Time Plot
        prisk <- riskData[,2]
        prisk[prisk > 0.999999] <- 0.999999
        timetoEventData$lammda <- -log(1.0-prisk) # From probability to risk
        
        aliveEvents <- timetoEventData
        atEventData <- subset(timetoEventData,timetoEventData$eStatus==1)
        atEventData <- atEventData[order(-atEventData$risk),]
        atEventData <- atEventData[order(atEventData$eTime),]
        timeInterval <- 2.0*mean(atEventData$eTime)
        maxtime <- max(atEventData$eTime)
        
        if (!is.null(riskTimeInterval))
        {
          timeInterval <- riskTimeInterval
        }
        Observed <- numeric(nrow(atEventData))
        Expected <- numeric(nrow(atEventData))
        timed <- numeric(nrow(atEventData))
        passAcum <- 0;
        for (idx in c(1:nrow(atEventData)))
        {
          timed[idx] <- atEventData[idx,"eTime"]
          Observed[idx] <- idx
          roevent <- aliveEvents$lammda*timed[idx]/timeInterval
          roevent[roevent > 1] <- 1
          roevent[aliveEvents$eStatus == 0] <- roevent[aliveEvents$eStatus == 0]*ExpectedNoEventsGain;
          Expected[idx] <- passAcum + sum(roevent);
          pssEvents <- subset(aliveEvents,aliveEvents$eTime <= timed[idx])
          aliveEvents <- subset(aliveEvents,aliveEvents$eTime > timed[idx])
          pnext <- pssEvents$lammda*timed[idx]/timeInterval
          pnext[pnext > 1.0] <- 1.0
          pnext[pssEvents$eStatus == 0] <- pnext[pssEvents$eStatus == 0]*ExpectedNoEventsGain
          passAcum <- passAcum+sum(pnext);
        }
        totObserved <- max(Observed)
        maxevents <- max(c(Observed,Expected))
        tokeep <- (Expected > 0.10*totObserved)
        tokeep <- tokeep & (Observed > 0.10*totObserved)
        OEData <- as.data.frame(cbind(time=timed,Observed=Observed,Expected=Expected,Included=tokeep))
        OEratio <- Observed[tokeep]/Expected[tokeep]
        OE95ci <- c(mean=mean(OEratio),metric95ci(OEratio))
        rownames(OEData) <- rownames(atEventData)

        observedCI <- stats::poisson.test(max(Observed[tokeep]),max(Expected[tokeep]), conf.level = 0.95 )
        OERatio <- observedCI
        OERatio$estimate <- c(OERatio$estimate,OERatio$conf.int,OERatio$p.value)
        names(OERatio$estimate) <- c("O/E","Low","Upper","p.value")
        
        if (plotRR)
        {
          par(mfrow=c(1,1))
          plot(timed,Expected,pch=4,type="b",cex=0.5,
             main=paste("Time vs. Events:",title),
             ylab="Events",
             xlab="Time",
             ylim=c(0,1.05*maxevents),
             xlim=c(0,1.05*maxtime))
          points(timed,Observed,pch=1,col="red")
          legend("topleft",legend=c("Expected","Observed"),pch=c(4,1),lty=c(1,0),col=c(1,"red"))
        }
      }
        ## Survival plot
      labelsplot <- c("Low",sprintf("At Risk > %5.3f",thr_atP[1]));
      paletteplot <- c("green", "red")
      if (length(thr_atP)>1)
      {
        labelsplot <- c(sprintf("Low Risk < %5.3f",thr_atP[2]),sprintf("%5.3f <= Risk < %5.3f",thr_atP[2],thr_atP[1]),
                        sprintf("High Risk >= %5.3f",thr_atP[1]));
        paletteplot <- c("green", "cyan","red")
      }
      LogRankE <- EmpiricalSurvDiff(times=timetoEventData$eTime,
                  status=timetoEventData$eStatus,
                  groups=timetoEventData$class,
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
                 keyPoints=thrPoints,
                 OERatio=OERatio,
                 OE95ci=OE95ci,
                 OARatio=OARatio,
                 OAcum95ci=OAcum95ci,
                 fit=lfit,
                 ROCAnalysis=ROCAnalysis,
                 prevalence=pre,
                 thr_atP=thr_values,
                 SEN_atP=sensitivity,
                 LowEventsFrac_atP=LowEventsFrac,
                 HighEventsFrac_atP=HighEventsFrac,
                 RR_atP=RRAtSen,
                 c.index=cstat,
                 surfit=surfit,
                 surdif=surdif,
                 LogRankE=LogRankE
                 )
  return (result)
}
