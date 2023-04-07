RRPlot <-function(riskData=NULL,atProb=c(0.90,0.80),atThr=NULL,title="",timetoEvent=NULL,ysurvlim=c(0,1.0))
{
  totobs <- nrow(riskData)
if (!requireNamespace("corrplot", quietly = TRUE)) {
		install.packages("corrplot", dependencies = TRUE)
		}
        
  numberofEvents <- sum(riskData[,1])
  numberofNoEvents <- sum(riskData[,1]==0)
  
  minRiskAtEvent <- min(riskData[riskData[,1]==1,2])
  risksGreaterThanM <- riskData[riskData[,2] >= minRiskAtEvent,2]
  names(risksGreaterThanM) <- rownames(riskData[riskData[,2] >= minRiskAtEvent,])
  risksGreaterThanM <- risksGreaterThanM[order(risksGreaterThanM)]
  nobs <- length(risksGreaterThanM)
  RR <- numeric(nobs)
  SEN <- numeric(nobs)
  SPE <- numeric(nobs)
  isEvent <- riskData[names(risksGreaterThanM),1]
#  print(head(isEvent))
  PPV <- numeric(nobs)
  idx <- 1;
  netBenefit <- numeric(nobs)
  pthr <- numeric(nobs)
  thrs <- numeric(nobs)
  mxthr <- max(riskData[,2])
  mithr <- min(riskData[,2])
  for (thr in risksGreaterThanM)
  {
    
    atLowRisk <- riskData[,2] < thr
    atHighRisk <- !atLowRisk
    LowEvents <- sum(riskData[atLowRisk,1])
    if (LowEvents==0) LowEvents <- 0.5;
    HighEvents <- sum(riskData[atHighRisk,1]);
    SEN[idx] <-  HighEvents/numberofEvents;
    SPE[idx] <- sum(riskData[atLowRisk,1]==0)/numberofNoEvents;
    PPV[idx] <- HighEvents/sum(atHighRisk)
    LowFraction <- LowEvents/sum(atLowRisk);
    HighFraction <- HighEvents/sum(atHighRisk);
    if ((sum(atLowRisk)/numberofNoEvents) > 0.01)
    {
      RR[idx] <- HighFraction/LowFraction
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
    pthr[idx] <- thr
    thrs[idx] <- thr
    netBenefit[idx] <- (HighEvents - sum(riskData[atHighRisk,1]==0)*pthr[idx]/(1.0001-pthr[idx]))/totobs;
    idx <- idx + 1;
  }
  colors <- heat.colors(10)
  
  if ((mithr>=0) && (mxthr<=1.0))
  {
    ## Observed vs Expected plot
    tmop <- par(no.readonly = TRUE)
    par(pty='s')
    xdata <- riskData[order(-riskData[,2]),]
    observed <- numeric(nrow(xdata))
    expected <- numeric(nrow(xdata))
    observed[1] <- xdata[1,1]
    expected[1] <- xdata[1,2]
    for (idx in c(2:nrow(xdata)))
    {
      expected[idx] <- expected[idx-1]+xdata[idx,2]
      observed[idx] <- observed[idx-1]+xdata[idx,1]
    }
    maxobs <- max(c(observed,expected))
    plot(expected,observed,
         ylab="Observed",
         xlab="Expected",
         xlim=c(0,maxobs),
         ylim=c(0,maxobs),
         main=paste("Expected vs. Observed:",title),
         pch=c(5+14*xdata[,1]),
         col=c(1+xdata[,1]),
         cex=c(0.2+xdata[,1]))
    lines(x=c(0,maxobs),y=c(0,maxobs),lty=2)
    par(tmop)
    
    ## Decision curve analysis
    pshape <- 4 + 12*isEvent
    pre=numberofEvents/totobs
    xmax <- quantile(thrs,probs=c(0.99))
    plot(thrs,netBenefit,main=paste("Decision Curve Analysis:",title),ylab="Net Benefit",xlab="Threshold",
         ylim=c(min(netBenefit),pre),
         xlim=c(0,xmax),
         pch=pshape,col=colors[1+floor(10*(1.0-SEN))],cex=(0.35 + isEvent))
    lfit <-loess(netBenefit~thrs,span=0.25);
    plx <- predict(lfit,se=TRUE)
    lines(thrs,plx$fit,lty=1)
    lines(thrs,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
    lines(thrs,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
    abline(h=0,col="blue")
    lines(x=c(0,pre),y=c(pre,0),col="red")
    legtxt <- sprintf("%3.1f",c(5:0)/5)
    colorlegend(heat.colors(10), legtxt,xlim=c(0.92*xmax,0.98*xmax),ylim=c(pre*0.45,pre*0.1),cex=0.75)
    text(0.92*xmax,pre*0.5,"SEN")
    
    legend("topright",legend=c("No Event","Event"),
           pch=c(4,16),
           col=c(1,1))
    legend("bottomleft",legend=c("Treat All","Treat None"),
           lty=c(1,1),
           col=c("red","blue"),cex=0.5)
  }
  ## Relative Risk plot
  
  ymax <- quantile(RR,probs=c(0.99))
  pshape <- 2 + 14*isEvent
  plot(SEN,RR,cex=(0.35 + PPV),
       pch=pshape,
       col=colors[1+floor(10*(1.0-SPE))],
       xlim=c(0,1.15),
       ylim=c(1.0,ymax+0.5),
#       log="y",
       main=paste("Relative Risk:",title))
  lfit <-loess(RR~SEN,span=0.5);
  plx <- predict(lfit,se=TRUE)
  lines(SEN,plx$fit,lty=1)
  lines(SEN,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
  lines(SEN,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

  legtxt <- sprintf("%3.1f",c(5:0)/5)
  colorlegend(heat.colors(10), legtxt,xlim=c(1.05,1.175),ylim=c((ymax-0.75)*0.45+1.0,ymax/10+1.0),cex=0.75)
  text(1.075,1,"SPE")
  sizp <- c(5:1)/5.0 + 0.35
  legend("topright",legend=legtxt[1:5],pch=16,cex=sizp)
  legend("bottomleft",legend=c("No","Yes"),pch=c(2,16),cex=0.5)
  text(0.95,ymax+0.5,"PPV->")


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
  }
  
  lowRisk <- riskData[,2] < thr_atP[1]
  LowEventsFrac <- sum(riskData[lowRisk,1])/sum(lowRisk)
  HighEventsFrac <- sum(riskData[!lowRisk,1])/sum(!lowRisk)
  sensitivity=sum(riskData[!lowRisk,1])/numberofEvents
  specificity=sum(riskData[lowRisk,1]==0)/numberofNoEvents
  abline(v=sensitivity,col="blue")
  text(x=sensitivity,y=ymax,sprintf("Index(%3.2f)=%4.3f",specificity,thr_atP[1]),pos=4 - 2*(sensitivity>0.5) ,cex=0.7)
  text(x=sensitivity,y=0.9*(ymax-1.0)+1.0,
       sprintf("RR(%3.2f)=%4.3f",
               sensitivity,
               predict(lfit,sensitivity)),
       pos=4 - 2*(sensitivity>0.5),cex=0.7)
  minthr <- min(risksGreaterThanM)
  maxthr <- max(risksGreaterThanM)
  rangethr <- minthr + (10:0)/10*(maxthr-minthr);
  surfit <- NULL
  LogRankE <- NULL

  if (!is.null(timetoEvent))
  {
    ## Survival Plot
      timetoEventData <- as.data.frame(cbind(event=riskData[,1],
                           class=1*(riskData[,2]>=thr_atP[1]),
                           time=timetoEvent,
                           risk=riskData[,2])
                           )
      if (length(thr_atP)>1)
      {
        timetoEventData$class <- 1*(riskData[,2]>=thr_atP[1]) + 1*(riskData[,2]>=thr_atP[2])
      }
    
 
      ## Survival plot
      labelsplot <- c("Low",sprintf("At Risk > %5.3f",thr_atP[1]));
      paletteplot <- c("green", "red")
      if (length(thr_atP)>1)
      {
        labelsplot <- c("Low",sprintf("%5.3f <= Risk < %5.3f",thr_atP[2],thr_atP[1]),
                        sprintf("Risk >= %5.3f",thr_atP[1]));
        paletteplot <- c("green", "cyan","red")
      }
      LogRankE <- EmpiricalSurvDiff(times=timetoEventData$time,
                  status=timetoEventData$event,
                  groups=timetoEventData$class,
                  plots=FALSE,main=paste("Kaplan-Meier:",title))
      
      surfit <- survival::survfit(Surv(time,event)~class,data = timetoEventData)
      surdif <- survival::survdiff(Surv(time,event)~class,data = timetoEventData)
      
      graph <- survminer::ggsurvplot(surfit,
                                     data=timetoEventData, 
                                     conf.int = TRUE, 
                                     legend.labs=labelsplot,
                                     palette = paletteplot,
                                     ylim = ysurvlim,
                                     ggtheme = ggplot2::theme_bw() + 
                                       ggplot2::theme(plot.title = ggplot2::element_text(
                                         hjust = 0.5, face = "bold", size = 20)),
                      title = paste("Kaplan-Meier:",title),
                      risk.table = TRUE,
                      tables.height = 0.2,
                      tables.theme = survminer::theme_cleantable())
      print(graph)
    

  }
  
  result <- list(EventsThr=risksGreaterThanM,
                 isEvent=isEvent,
                 RR=RR,
                 SEN=SEN,
                 SPE=SPE,
                 PPV=PPV,
                 fit=lfit,
                 thr_atP=thr_atP,
                 SEN_atP=sensitivity,
                 LowEventsFrac_atP=LowEventsFrac,
                 HighEventsFrac_atP=HighEventsFrac,
                 RR_atP=predict(lfit,sensitivity),
                 surfit=surfit,
                 sufdif=surdif,
                 LogRankE=LogRankE
                 )
  return (result)
}
