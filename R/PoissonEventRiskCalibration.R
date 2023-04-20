CalibrationProbPoissonRisk <- function(Riskdata)
### Riskdata is a matrix of a Poisson event with Event, Probability of Event>0, and Time to Event]
{

  Riskdata <- as.data.frame(Riskdata)
  colnames(Riskdata) <- c("Event","pGZ","Time")
  observed <- sum(Riskdata$Event)
  meaninterval <- mean(Riskdata[Riskdata$Event==1,"Time"]);
  h0 <- 1.0/meaninterval
  index <- log(-log(1.0-Riskdata$pGZ)/h0)
  hazard <- h0*exp(index)
  probGZero <- 1.0-exp(-hazard)
  expected <- sum(probGZero)
  delta <- abs(observed-expected)/observed
  while (delta>0.005)
  {
    gain <- expected/observed
    h0 <- h0/gain
    hazard <- h0*exp(index)
    probGZero <- 1.0-exp(-hazard)
    expected <- sum(probGZero)
    delta <- abs(observed-expected)/observed
  }
  timeInterval <- 2.0*meaninterval
  timeSorted <- Riskdata
  timeSorted$hazard <- hazard
  timeSorted <- timeSorted[order(Riskdata$Time),]
  Ahazard <- 0
  for (idx in c(1:nrow(timeSorted)))
  {
    nevent <- timeSorted[idx,"hazard"]*timeSorted[idx,"Time"]/timeInterval
    if (nevent > 1.0) nevent <- 1.0
    Ahazard <- Ahazard + nevent
  }
  timeInterval <- timeInterval*Ahazard/observed
  Ahazard <- 0
  for (idx in c(1:nrow(timeSorted)))
  {
    nevent <- timeSorted[idx,"hazard"]*timeSorted[idx,"Time"]/timeInterval
    if (nevent > 1.0) nevent <- 1.0
    Ahazard <- Ahazard + nevent
  }
  
  result <- list(index=index,
                 probGZero=probGZero,
                 hazard=hazard,
                 h0=h0,
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



