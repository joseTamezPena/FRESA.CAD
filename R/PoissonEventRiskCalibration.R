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

meanTimeToEvent <- function(probGZero,timeInterval)
{
  meanEvents <- -log(1.0-probGZero)
  return (0.5*timeInterval/meanEvents)
}


CalibrationProbPoissonRisk <- function(Riskdata,trim=0.10,timeInterval=NULL)
### Riskdata is a matrix of a Poisson event with Event, Probability of Event>0, and Time to Event]
{

  Riskdata <- as.data.frame(Riskdata)
  colnames(Riskdata) <- c("Event","pGZ","Time")
  observed <- sum(Riskdata$Event)
  
  meaninterval <- mean(subset(Riskdata,Riskdata$Event==1)$Time);
  if (is.null(timeInterval))
  {
    timeInterval <- 2*meaninterval;
  }
  
  h0 <- sum(Riskdata$Event & Riskdata$Time <= timeInterval)
  h0 <- h0/sum((Riskdata$Time > timeInterval) | (Riskdata$Event==1))
  
#  print(c(h0,timeInterval))
  
  probGZero <- Riskdata$pGZ
  probGZero[probGZero > 0.999999] <- 0.999999

  hazard <- -log(1.00-probGZero)

   ## Adjust probabilites of no-event that share similar time than events 
  atRisktime <- 3.0*meaninterval
  isnoevent <- (Riskdata$Time < atRisktime) & (Riskdata$Event==0) # Adjust only if short time to event on censored events
  probGZero[isnoevent] <- 1.0 - exp(-hazard[isnoevent]*Riskdata[isnoevent,"Time"]/timeInterval)

  index <- log(-log(1.00-probGZero)/h0)
  hazard <- -log(1.00-probGZero)
  
  
  expected <- sum(probGZero)
  delta <- abs(observed-expected)/observed
  totgain <- 1.0;
  eindex <- exp(index)
  while (delta>0.005)
  {
    gain <- expected/observed
#    print(c(gain,totgain))
    h0 <- h0/gain
    hazard <- h0*eindex
    probGZero <- 1.0-exp(-hazard)
    expected <- sum(probGZero)
    delta <- abs(observed-expected)/observed
    totgain <- totgain/gain
  }
  probGZero <- adjustProb(Riskdata$pGZ,totgain)
  timeSorted <- Riskdata
  timeSorted$pGZ <- probGZero
  probGZero[probGZero > 0.999999] <- 0.999999
  timeSorted$hazard <- -log(1.0-probGZero)
  timeSorted <- timeSorted[order(-timeSorted$hazard),]
  timeSorted <- timeSorted[order(timeSorted$Time),]
#  print(head(timeSorted))
  touse <- c(1:nrow(timeSorted))
  firstrimObserved <- trim*observed
  lasttrimObserved <- (1.0 - 0.5*trim)*observed 
  timed <- unique(timeSorted[timeSorted$Event==1,"Time"])
  allTimes <- timeSorted$Time
#  for (lp in 1:3)
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
          deltatime <- (timed[idx] - lasttime)
          totObs <- lastObs + sum(timeSorted[whosum,"Event"])
          lastObs <- totObs
          eevents <- sum(timeSorted[allTimes > lasttime,"hazard"])*deltatime/timeInterval
          acuHazard <- passHazard + eevents;
          passHazard <- acuHazard;
          lasttime <- timed[idx]
          if ((totObs >= firstrimObserved) & (totObs <= lasttrimObserved))
          {
            wt <- sqrt(totObs)
            gainAdded <- gainAdded + wt
            meanGain <- meanGain + wt*(acuHazard/totObs)
          }
#          cat(totObs,",",eevents,",",acuHazard,"\n")
    }
    if (gainAdded > 0)
    {
      meanGain <- meanGain/gainAdded
    }
    else
    {
      meanGain <- acuHazard/totObs
    }
    timeInterval <- timeInterval*meanGain
#    cat("(",timeInterval,",",gainAdded,",",meanGain,",",totObs,",",acuHazard,")")
  }
  nevent <- timeSorted$"hazard"*timeSorted$Time/timeInterval
  nevent[nevent > 1.0]  <- 1.0
  Ahazard <- sum(nevent)

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
    timeInterval <- 2*meaninterval;
  }
  h0 <- sum(data[,outcome]==1 & data[,time] <= timeInterval)
  h0 <- h0/sum((data[,time] > timeInterval) | (data[,outcome]==1))
  
  
  hazard <- h0*exp(index)
  probGZero <- 1.0-exp(-hazard)

  Riskdata <- cbind(data[,outcome],probGZero,data[,time])
  result <- CalibrationProbPoissonRisk(Riskdata,trim,timeInterval)
  

  return (result)
}



