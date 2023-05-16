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
#  startTimeInterval <- timeInterval
  timeSorted <- Riskdata
  timeSorted$pGZ <- probGZero
  probGZero[probGZero > 0.999999] <- 0.999999
  timeSorted$hazard <- -log(1.0-probGZero)
  timeSorted <- timeSorted[order(-timeSorted$hazard),]
  timeSorted <- timeSorted[order(timeSorted$Time),]
#  print(head(timeSorted))
  touse <- c(1:nrow(timeSorted))
  firstrimObserved <- trim*observed + 1
  lasttrimObserved <- (1.0-trim)*observed + 1
  for (lp in 1:3)
  {
    Ahazard <- 0
    cobserved <- 0
    deltaObs <- 0;
    deltaAdded <- 0;
    alive <- timeSorted;
    pastEvents <- 0;
    idx <- 0
#    cat("Rows",nrow(alive),"Time:",idx,">")
    for (idx in touse)
    {
      oevent <- alive$hazard*(timeSorted[idx,"Time"]/timeInterval)
      oevent[oevent > 1.0] <- 1.0
      Ahazard <- pastEvents + sum(oevent)
      cobserved <- cobserved + timeSorted[idx,"Event"]
      if ((cobserved > firstrimObserved) && (cobserved < lasttrimObserved))
      {
        deltaObs <- deltaObs + Ahazard/cobserved;
        deltaAdded <- deltaAdded + 1.0;
#        cat("Rows",nrow(alive),"Time:",idx,">")
      }
      donealive <- alive[alive$Time <= timeSorted[idx,"Time"],]
      pevent <- donealive$hazard*timeSorted[idx,"Time"]/timeInterval
      pevent[pevent > 1.0] <- 1.0
      pastEvents <- pastEvents + sum(pevent)
#      cat(":",timeSorted[idx,"Time"],",",sum(alive$Time > timeSorted[idx,"Time"]),"\n")
      alive <- alive[alive$Time > timeSorted[idx,"Time"],]
    }
    if (deltaAdded > 0)
    {
      gainTime <- deltaObs/deltaAdded
    }
    else
    {
      gainTime <- Ahazard/observed
    }
#    print(c(Ahazard/observed,gainTime))
    timeInterval <- timeInterval*gainTime
#    cat(deltaObs,":",Ahazard,":",cobserved,":",timeInterval,"\n")
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



