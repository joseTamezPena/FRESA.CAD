RRPlot <-function(riskData=NULL,atProb=0.90,title="Relative Risk",timetoEvent=NULL,titleS="Kaplan-Meier")
{
if (!requireNamespace("corrplot", quietly = TRUE)) {
		install.packages("corrplot", dependencies = TRUE)
		}
        
  numberofEvents <- sum(riskData[,1])
  numberofNoEvents <- sum(riskData[,1]==0)
  
  minRiskAtEvent <- min(riskData[riskData[,1]==1,2])
  risksGreaterThanM <- riskData[riskData[,2] >= minRiskAtEvent,2]
  risksGreaterThanM <- risksGreaterThanM[order(risksGreaterThanM)]
  nobs <- length(risksGreaterThanM)
  RR <- numeric(nobs)
  PPV <- numeric(nobs)
  NPV <- numeric(nobs)
  isEvent <- riskData[names(risksGreaterThanM),1]
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
    RR[idx] <- HighFraction/LowFraction
    idx <- idx + 1;
  }
  ymax <- max(RR)
  colors <- heat.colors(10)
  pshape <- 1 + 15*isEvent
  plot(PPV,RR,cex=(0.3 + Sensitivity),
       pch=pshape,
       col=colors[1+floor(10*(1.0-NPV))],
       xlim=c(0,1.15),
       ylim=c(1.0,ymax+0.5),
       main=title)
  lfit <-loess(RR~PPV);
  plx <- predict(lfit,se=TRUE)
  lines(PPV,plx$fit,lty=1)
  lines(PPV,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
  lines(PPV,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

  legtxt <- sprintf("%3.1f",c(5:0)/5)
  colorlegend(heat.colors(10), legtxt,xlim=c(1.05,1.175),ylim=c(ymax/2,1),cex=0.75)
  text(1.075,ymax/2+0.25,"NPV")
  sizp <- c(5:1)/5.0 + 0.3
  legend("topright",legend=legtxt[1:5],pch=16,cex=sizp)
  text(0.95,ymax,"Sen->")


  thr_atP <- quantile(riskData[riskData[,1]==0,2],probs=c(atProb))
  lowRisk <- riskData[,2] < thr_atP
  LowEventsFrac <- sum(riskData[lowRisk,1])/sum(lowRisk)
  HighEventsFrac <- sum(riskData[!lowRisk,1])/sum(!lowRisk)
  precision=sum(riskData[!lowRisk,1])/numberofEvents
  abline(v=precision,col="blue")
  text(x=precision,y=ymax,sprintf("Index(%3.2f)=%4.3f",atProb,thr_atP),pos=4,cex=0.7)
  surfit <- NULL
  LogRankE <- NULL

  if (!is.null(timetoEvent))
  {
      timetoEventData <- as.data.frame(cbind(event=riskData[,1],
                           class=1*(riskData[,2]>=thr_atP),
                           time=timetoEvent))
    
      labelsplot <- c("Low",sprintf("At Risk > %5.2f",thr_atP));
      paletteplot <- c("#00bbff", "#ff0000")
      
      LogRankE <- EmpiricalSurvDiff(times=timetoEventData$time,
                  status=timetoEventData$event,
                  groups=timetoEventData$class,
                  plots=FALSE,main=titleS)
      
      surfit <- survival::survfit(Surv(time,event)~class,data = timetoEventData)
      
      graph <- survminer::ggsurvplot(surfit,
                                     data=timetoEventData, 
                                     conf.int = TRUE, 
                                     legend.labs=labelsplot,
                                     palette = paletteplot,
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
                 LogRankE=LogRankE
                 )
  return (result)
}
