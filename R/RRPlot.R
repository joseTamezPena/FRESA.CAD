RRPlot <-function(riskData=NULL,atProb=c(0.90,0.80),atThr=NULL,title="Relative Risk",timetoEvent=NULL,titleS="Kaplan-Meier",ysurvlim=c(0,1.0))
{
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
  PPV <- numeric(nobs)
  NPV <- numeric(nobs)
  isEvent <- riskData[names(risksGreaterThanM),1]
#  print(head(isEvent))
  Sensitivity <- numeric(nobs)
  idx <- 1;
  for (thr in risksGreaterThanM)
  {
    
    atLowRisk <- riskData[,2] < thr
    atHighRisk <- !atLowRisk
    LowEvents <- sum(riskData[atLowRisk,1])
    if (LowEvents==0) LowEvents <- 0.5;
    HighEvents <- sum(riskData[atHighRisk,1]);
    PPV[idx] <-  HighEvents/numberofEvents;
    NPV[idx] <- sum(riskData[atLowRisk,1]==0)/numberofNoEvents;
    Sensitivity[idx] <- HighEvents/sum(atHighRisk)
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
    idx <- idx + 1;
  }
  ymax <- quantile(RR,probs=c(0.99))
  colors <- heat.colors(10)
  pshape <- 2 + 14*isEvent
  
  plot(PPV,RR,cex=(0.35 + Sensitivity),
       pch=pshape,
       col=colors[1+floor(10*(1.0-NPV))],
       xlim=c(0,1.15),
       ylim=c(1.0,ymax+0.5),
#       log="y",
       main=title)
  lfit <-loess(RR~PPV,span=0.5);
  plx <- predict(lfit,se=TRUE)
  lines(PPV,plx$fit,lty=1)
  lines(PPV,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
  lines(PPV,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

  legtxt <- sprintf("%3.1f",c(5:0)/5)
  colorlegend(heat.colors(10), legtxt,xlim=c(1.05,1.175),ylim=c((ymax-0.75)*0.45+1.0,ymax/10+1.0),cex=0.75)
  text(1.075,1,"NPV")
  sizp <- c(5:1)/5.0 + 0.35
  legend("topright",legend=legtxt[1:5],pch=16,cex=sizp)
  legend("bottomleft",legend=c("No","Yes"),pch=c(2,16),cex=0.5)
  text(0.95,ymax+0.5,"Sen->")


  if (is.null(atThr))
  {
    thr_atP <- quantile(riskData[riskData[,1]==0,2],probs=c(atProb))
    if (length(atProb)>1)
    {
      if (atProb[2]>atProb[1])
      {
        thr_atP <- quantile(riskData[riskData[,1]==1,2],probs=c(1.0-atProb))
        print(thr_atP)
  #      atProb <- atProb[2]
      }
  #    atProb <- atProb[1]
    }
  }
  else
  {
    thr_atP <- atThr
  }
  print(thr_atP)
  
  lowRisk <- riskData[,2] < thr_atP[1]
  LowEventsFrac <- sum(riskData[lowRisk,1])/sum(lowRisk)
  HighEventsFrac <- sum(riskData[!lowRisk,1])/sum(!lowRisk)
  precision=sum(riskData[!lowRisk,1])/numberofEvents
  abline(v=precision,col="blue")
  text(x=precision,y=ymax,sprintf("Index(%3.2f)=%4.3f",atProb[1],thr_atP[1]),pos=4 - 2*(precision>0.5) ,cex=0.7)
  text(x=precision,y=0.9*(ymax-1.0)+1.0,
       sprintf("RR(%3.2f)=%4.3f",
               precision,
               predict(lfit,precision)),
       pos=4 - 2*(precision>0.5),cex=0.7)
  minthr <- min(risksGreaterThanM)
  maxthr <- max(risksGreaterThanM)
  rangethr <- minthr + (10:0)/10*(maxthr-minthr);
#  axis(3, at=c(0:10)/10,labels=round(rangethr,digits=2),line=0,outer=FALSE,col.axis="gray", las=1, cex.axis=0.7, tck=-.01)
#  mtext("Index", side=3, line=0, cex.lab=0.5,las=1, col="gray")
  surfit <- NULL
  LogRankE <- NULL

  if (!is.null(timetoEvent))
  {
      timetoEventData <- as.data.frame(cbind(event=riskData[,1],
                           class=1*(riskData[,2]>=thr_atP[1]),
                           time=timetoEvent))
      if (length(thr_atP)>1)
      {
        timetoEventData$class <- 1*(riskData[,2]>=thr_atP[1]) + 1*(riskData[,2]>=thr_atP[2])
      }
    
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
                  plots=FALSE,main=titleS)
      
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
                      title = titleS,
                      risk.table = TRUE,
                      tables.height = 0.2,
                      tables.theme = survminer::theme_cleantable())
      print(graph)
    

  }
  
  result <- list(EventsThr=risksGreaterThanM,
                 isEvent=isEvent,
                 RR=RR,
                 PPV=PPV,
                 NPV=NPV,
                 Sensitivity=Sensitivity,
                 fit=lfit,
                 thr_atP=thr_atP,
                 PPV_atP=precision,
                 LowEventsFrac_atP=LowEventsFrac,
                 HighEventsFrac_atP=HighEventsFrac,
                 RR_atP=predict(lfit,precision),
                 surfit=surfit,
                 sufdif=surdif,
                 LogRankE=LogRankE
                 )
  return (result)
}
