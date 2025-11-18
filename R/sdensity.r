# This function will estimate the probability density of a random variable based on the sampled data.
# x : sample observations
# lambda: degree of smoothing (0,1]
# at return it will provide the spline function of the pdf
sdensity <- function(x,lambda=0.5,npts=512,posSdef=FALSE,from=NULL,to=NULL)
{
  IQRRange <- IQR(x)
  maxval <- max(x)
  minval <- min(x)
  sfilt <- c(1,1,1)/3 # Filter for smoothing the pdf and estimation of the CDF
  ddpts <- 11 # Extra points
  hhd <- x + 1.0e-10*rnorm(length(x))/(max(x)-min(x)) # Small perturbation to break ties
  xval <- unique(hhd)
  xval <- xval[order(xval)]

  # Add tail samples to the left and right of the distribution
  xtail <- c(1:ddpts)*IQRRange/ddpts
  llim <- minval - xtail[1]
  uplim <- maxval + xtail[1]
  xval <- c(llim - rev(xtail),xval,uplim + xtail)
  if (posSdef)
  {
    if (min(xval)<0)
    {
      xval <- c(-xtail[1],xval[xval >= 0])
    }
  }
  # Add an dummy sample at the end
  hhd <- c(hhd,maxval + 0.5*xtail[1])
  # The CDF
  CDFh <- numeric(length(xval))
  for (i in c(1:length(xval)))
  {
    CDFh[i] <- sum(1*(hhd<=xval[i]))
  }
  rawX <- xval
  CDFh <- CDFh/max(CDFh)
  rawCDF <- CDFh
  # The interpolation function of the CDF for monotonic functions
  sCDFh <- splinefun(xval,CDFh,"hyman")

  # The regular interval for the pdf
  smin <- min(xval)
  smax <- max(xval)
  # Estimating the number of step of the pyramidal estimation of the CDF
  fres <- 1 + round(log(1.0 + npts/length(hhd))/log(2))
  nxval <- rawX
  CDFhff <- CDFh
  pxval <- rawX
  tfilt <- 3
  snpts <- length(nxval)
  # Smoothing the empirical CDF from the original resolution to the specified resolution
  for (sc in c(1:fres))
  {
    if (snpts < npts)
    {
      nxval <- c(nxval,pxval,x) # The possible X observations of the CDF
      nxval <- unique(nxval)
      nxval <- nxval[order(nxval)]
      pxval <- nxval
      CDFhff <- sCDFh(nxval) # Get the CDF
      ndpts <- length(CDFhff)
      for (n in c(1:tfilt)) # Smoothing the CDF and the X support
      {
        CDFhf <- as.numeric(filter(CDFhff,sfilt))
        CDFhff[2:(ndpts-1)] <- CDFhf[2:(ndpts-1)]
        nxvalf <- as.numeric(filter(nxval,sfilt))
        nxval[2:(ndpts-1)] <- nxvalf[2:(ndpts-1)]
      }
      sCDFh <- splinefun(nxval,CDFhff,"hyman") # The new CDF function
    }
    snpts <- length(nxval)
  }
  rxval <- nxval # Store the used X support of the CDF

  smin <- min(rxval)
  smax <- max(rxval)
  sCDFh <- splinefun(rxval,CDFhff,"hyman") # The final CDF function
  nxval <- smin+c(1:npts)*(smax-smin)/npts
  snpts <- npts

  # the pdf is the derivative of the CDF
  dxf <- 0.5*(nxval[2]-nxval[1])
  pdfh <- (sCDFh(nxval+dxf) - sCDFh(nxval-dxf))/(2*dxf) # Numerical differentiation
  # the tails of the pdf are zeros
  pdfh <- c(rep(0,ddpts),pdfh,rep(0,ddpts))
  ndpts <- length(nxval)
  nxval <- c(nxval[1] - c(ddpts:1)*dxf,
             nxval,
             nxval[ndpts] + c(1:ddpts)*dxf)
  masslowzero <- 0
  atzero <- nxval[1]
  if (posSdef) # Correct for positive definitive
  {
    masslowzero <- sum(pdfh[nxval<0])
    atzero <- nxval[nxval>=0]
    atzero <- nxval == atzero[1]
    if (masslowzero > 0)
    {
      pdfh[nxval<0] <- 0
      pdfh[atzero] <- pdfh[atzero] + masslowzero
    }
  }
  # Smoothing the pdf using recursive filtering
  ndpts <- length(pdfh)
  tfilt <- 1+round(2*IQRRange/dxf*lambda) # The number of steps depends on the IQR of the observations
  for (n in c(1:tfilt))
  {
    pdfhf <- as.numeric(filter(pdfh,sfilt))
    pdfh[2:(ndpts-1)] <- pdfhf[2:(ndpts-1)]
    if (posSdef) # Adjust for positive definitive 
    {
      masslowzero <- sum(pdfh[nxval<0])
      if (masslowzero>0)
      {
        pdfh[nxval<0] <- 0
        pdfh[atzero] <- pdfh[atzero] + masslowzero
      }
    }
  }

  # Lets get the interpolation function of the pdf
  deltax <- c(1:round(0.25*ndpts))*IQRRange/round(0.25*ndpts);
  frangex <- c(min(nxval)-rev(deltax),nxval,max(nxval)+deltax)
  pdfext <- c(rep(0,length(deltax)),pdfh,rep(0,length(deltax)))

  pdf=splinefun(frangex,pdfext,"natural")

  # The pdf for the input specified sampling
  if (is.null(from))
  {
    from <- smin
  }
  if (is.null(to))
  {
    to <- smax
  }
  nxval <- from+c(1:npts)*(to-from)/npts
  pdfh <- pdf(nxval)
  pdfh[nxval<smin] <- 0;
  pdfh[nxval>smax] <- 0;
  CDF.f <-sapply(nxval,CDFint,pdf,ll=(min(nxval)-dxf))
  CDF.f[CDF.f<0] <- 0
  
  
  # the density analysis
  result <- list(x=nxval,
                 y=pdfh,
                 d.f=CDF.f,
                 pdf.f=pdf,
                 eCDF.f=sCDFh,
                 x.eCDF = rawX,
                 eCDF = rawCDF,
                 res.x = rxval
                 )
  return(result)
  
}