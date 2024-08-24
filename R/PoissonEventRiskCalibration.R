adjustProb <- function(probGZero,gain)
{
  probGZero[probGZero > 0.999999] <- 0.999999
  hazard <- -log(1.0-probGZero)*gain
  probGZero <- 1.0-exp(-hazard)
  return (probGZero)
}

ppoisGzero <- function(index,h0)
{
  hazard <- h0*exp(index)
  probGZero <- 1.0-exp(-hazard)
  return (probGZero)
}

expectedEventsPerInterval <- function(probGZero)
{
  probGZero[probGZero > 0.999999] <- 0.999999
  return (-log(1.0-probGZero))
}

meanTimeToEvent <- function(probGZero,timeInterval)
{
  meanEvents <- expectedEventsPerInterval(probGZero)
  return (0.5*timeInterval/meanEvents)
}


CalibrationProbPoissonRisk <- function(Riskdata,trim=0.10)
### Riskdata is a matrix of a Poisson event with Event, Probability of Event>0, and Time to Event]
{

  Riskdata <- as.data.frame(Riskdata)
  colnames(Riskdata) <- c("Event","pGZ","Time")
  observed <- sum(Riskdata$Event)
  
  meaninterval <- mean(subset(Riskdata,Riskdata$Event==1)$Time);
  timeInterval <- meaninterval;
  
  
  probGZero <- Riskdata$pGZ
  probGZero[probGZero > 0.999999] <- 0.999999

  hazard <- -log(1.00-probGZero)
  h0 <- mean(Riskdata$Event)

  index <- log(-log(1.00-probGZero)/h0)
  hazard <- -log(1.00-probGZero)
  
  
  expected <- sum(probGZero)
  delta <- abs(observed-expected)/observed
  totgain <- 1.0;
  eindex <- exp(index)
#  print(c(expected,observed,h0,totgain,delta))
  while (delta>0.001)
  {
    gain <- expected/observed
    h0 <- h0/gain
    hazard <- h0*eindex
    probGZero <- 1.0-exp(-hazard)
    expected <- sum(probGZero)
    delta <- abs(observed-expected)/observed
    totgain <- totgain/gain
#    print(c(h0,gain,totgain,delta))
  }
  probGZero <- adjustProb(Riskdata$pGZ,totgain)
  expected <- sum(probGZero)
  delta <- abs(observed-expected)/observed
#  print(c(expected,observed,h0,totgain,delta))

  timeSorted <- Riskdata
  timeSorted$pGZ <- probGZero
  probGZero[probGZero > 0.999999] <- 0.999999
  timeSorted$hazard <- -log(1.0-probGZero)
  timeSorted <- timeSorted[order(-timeSorted$hazard),]
  timeSorted <- timeSorted[order(timeSorted$Time),]
#  print(head(timeSorted))
  touse <- c(1:nrow(timeSorted))
  firstrimObserved <- trim*observed
  lasttrimObserved <- (1.0 - trim)*observed 
  timed <- unique(timeSorted[timeSorted$Event==1,"Time"])
  allTimes <- timeSorted$Time
#  for (lp in 1:2)
  {
    lasttime <- 0;
    lastObs <- 0;
    totObs <- 0;
    acuHazard <- 0;
    passHazard <- 0;
    gainAdded <- 0;
    meanGain <- 0;
    for (idx in c(1:length(timed)))
    {
          whosum <- (allTimes <= timed[idx]) & (allTimes > lasttime)
          whohazard <- allTimes > lasttime
          totObs <- lastObs + sum(timeSorted[whosum,"Event"])
          lastObs <- totObs
          deltatime <- pmin(allTimes[whohazard],rep(timed[idx],sum(whohazard))) - lasttime
          eevents <- sum(timeSorted[whohazard,"hazard"]*deltatime)/timeInterval
          acuHazard <- passHazard + eevents;
          passHazard <- acuHazard;
          lasttime <- timed[idx]
          if ((totObs >= firstrimObserved) & (totObs <= lasttrimObserved))
          {
            wt <- log(totObs)
            gainAdded <- gainAdded + wt
            meanGain <- meanGain + wt*(acuHazard/totObs)
          }
#          cat(totObs,",",eevents,",",acuHazard,",",min(deltatime),",",max(deltatime),"\n")
    }
    if (gainAdded > 0)
    {
      meanGain <- meanGain/gainAdded
    }
    else
    {
      meanGain <- acuHazard/totObs
    }
    cat("(",timeInterval,",",gainAdded,",",meanGain,",",totObs,",",acuHazard,")\n")
  }
  timeInterval <- timeInterval*meanGain;
  nevent <- timeSorted$hazard*timeSorted$Time/timeInterval
  Ahazard <- sum(nevent)
  cat("(Delta Time:",timeInterval,",Tot:",totObs,", Tot Hazard:",Ahazard,")\n")

  result <- list(index=index,
                 probGZero=probGZero,
                 hazard=hazard,
                 h0=h0,
                 hazardGain=totgain,
                 timeInterval = timeInterval,
                 meaninterval = meaninterval,
                 Ahazard=Ahazard,
                 delta=delta)
  return (result)
}

CoxRiskCalibration <- function(ml,data,outcome,time,trim=0.10,timeInterval=NULL)
{
  index <- predict(ml,data)
  

  if (is.null(timeInterval))
  {
    meaninterval <- mean(data[data[,outcome]==1,time]);
    timeInterval <- meaninterval;
  }
#  h0 <- sum(data[,outcome]==1 & data[,time] <= timeInterval)
#  h0 <- h0/sum((data[,time] > timeInterval) | (data[,outcome]==1))
  
  h0 <- 1.0

  
  hazard <- h0*exp(index)
  probGZero <- 1.0-exp(-hazard)

  Riskdata <- cbind(data[,outcome],probGZero,data[,time])
  result <- CalibrationProbPoissonRisk(Riskdata,trim)
  

  return (result)
}



