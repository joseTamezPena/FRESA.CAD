adjustProb <- function(probGZero,gain)
{
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

CalibrationProbPoissonRisk <- function(Riskdata,trim=0.10)
### Riskdata is a matrix of a Poisson event with Event, Probability of Event>0, and Time to Event]
{

  Riskdata <- as.data.frame(Riskdata)
  colnames(Riskdata) <- c("Event","pGZ","Time")
  observed <- sum(Riskdata$Event)
  meaninterval <- mean(Riskdata[Riskdata$Event==1,"Time"]);
  h0 <- 1.0/meaninterval
  probGZero <- Riskdata$pGZ
  index <- log(-log(1.0-probGZero)/h0)
  hazard <- -log(1.0-probGZero)
  expected <- sum(probGZero)
  delta <- abs(observed-expected)/observed
  totgain <- 1.0;
  while (delta>0.005)
  {
    gain <- expected/observed
    h0 <- h0/gain
    hazard <- h0*exp(index)
    probGZero <- 1.0-exp(-hazard)
    expected <- sum(probGZero)
    delta <- abs(observed-expected)/observed
    totgain <- totgain/gain
#    print(totgain)
  }
  timeInterval <- 2.0*meaninterval
  timeSorted <- Riskdata
  timeSorted$pGZ <- probGZero
  timeSorted$hazard <- -log(1.0-probGZero)
  timeSorted <- timeSorted[order(-timeSorted$hazard),]
  timeSorted <- timeSorted[order(timeSorted$Time),]
#  print(head(timeSorted))
  touse <- c(1:nrow(timeSorted))
  firstrimObserved <- trim*observed
  lasttrimObserved <- (1.0-trim)*observed
  for (lp in 1:2)
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
      if ((cobserved >= firstrimObserved) && (cobserved <= lasttrimObserved))
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

CoxRiskCalibration <- function(ml,data,outcome,time)
{
  index <- predict(ml,data)
  meaninterval <- mean(data[data[,outcome]==1,time]);
  h0 <- 1.0/meaninterval
  hazard <- h0*exp(index)
  probGZero <- 1.0-exp(-hazard)

  Riskdata <- cbind(data[,outcome],probGZero,data[,time])
  result <- CalibrationProbPoissonRisk(Riskdata)
  

  return (result)
}



