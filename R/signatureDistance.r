signatureDistance <- 
function (template, data=NULL, method = c("pearson","spearman","kendall","RSS","MAN"),fwts=NULL)
{

#given the template: mean,median,sample, etc....;signatureDistance it will return the distance between the template to each row of the dataframe
#the template is a named numeric vector
#the data is a colnamed data frame
#methods:
# RSS: Normalized Root Sum Square
# MAN: Normalized Manhattan distance
# pearson: 2*(1-Pearson correlation coefficient)
# spearman: 2*(1-spearman correlation coefficient)
# kendall: 2*(1-kendall correlation coefficient)

	method <- match.arg(method)
	theQuant <- c(0.025,0.16,0.25,0.5,0.75,0.84,0.975);
	samplesize <- 2.0;
	if (class(template)[1] == "list")
	{
		samplesize <- template$samples;
		theQuant <- template$quant;
		meant <- template$meanv;
		sdt <- template$sdv;
		template <- template$template;
	}
	if (is.null(fwts))
	{
		fwts <- rep(1,ncol(template));
	}

	wvalues <- 1.0/abs(qt(theQuant,df=samplesize-1));
	
	if (class(template)[1]=="matrix")
	{
		vnames <- colnames(template);
	}
	else
	{
		vnames <- names(template);
	}
	datasubset <- as.matrix(data[,vnames]);
	
	medianv <- as.integer((length(theQuant) + 1)/2);
	ld <- NULL;
	ud <- NULL;
	qld <- NULL;
	qud <- NULL;
	if (class(template)[1] == "matrix")
	{
		tem <- meant;
		
#		cat("median:")
#		print(tem)
		wts <- numeric(length(tem));
		ld <- numeric(length(tem));
		for (i in 1:(medianv - 1))
		{
			tdis <- tem - template[i,];
			w <- (tdis >= 0)*theQuant[i];
			wts <- wts + w;
			tdis[tdis < 0] <- 0;
			ld <- ld + w*wvalues[i]*tdis;
		}
		wts[wts == 0] <- 1.0e-10;
		ld <- ld/wts;
		tdis <- tem - template[medianv - 2,];
		tdis[tdis < 0] <- 0;
		qld <- tdis*wvalues[medianv - 2];

		wts <- numeric(length(tem));
		ud <- numeric(length(tem));
		for (i in (medianv + 1):length(wvalues))
		{
			tdis <- template[i,] - tem;
			w <- (tdis >= 0)*(1.0-theQuant[i]);
			wts <- wts + w;
			tdis[tdis < 0] <- 0;
			ud <- ud + w*wvalues[i]*tdis;
		}
		wts[wts == 0] <- 1.0e-10;
		ud <- ud/wts;
		tdis <- template[medianv + 2,] - tem;
		tdis[tdis < 0] <- 0;
		qud <- tdis*wvalues[medianv + 2];

		ld[ld == 0] <- 0.5*ud[ld == 0];
		ld[ld == 0] <- sdt[ld == 0];
		ld[ld == 0] <- 0.25;
		qld[qld == 0] <- ld[qld == 0];

		ud[ud == 0] <- 0.5*ld[ud == 0];
		ud[ud == 0] <- sdt[ud == 0];
		ud[ud == 0] <- 0.25;
		qud[qud == 0] <- ud[qud == 0];

#		cat("ld:")
#		print(ld)
#		cat("qld:")
#		print(qld)
#		cat("ud:")
#		print(ud)
#		cat("qud:")
#		print(qud)
	}
	else
	{
		tem <- template;
		ld <- sd(template);
		ld[ld == 0] <- 0.25;

		ud <- ld;
		qld <- IQR(template)/abs(qnorm(0.25))/2;
		qld[qld == 0] <- ld[qld == 0];
		qud <- qld;

	}
	switch(method, 
		RSS = 
		{ 
			RSSDistance <- function (x,template,ld,ud,wts) 
			{
				md <- (x-template);
				tsum = sum(wts);
				md <- sqrt(sum(wts*(pmax(md/ud,-md/ld)^2),na.rm=TRUE)/tsum);
				return (md)
			}
			metric <- apply(datasubset,1,RSSDistance,tem,ld,ud,fwts);
		},
		MAN = 
		{ 
			manDistance <- function (x,template,ld,ud,wts) 
			{
				md <- (x-template)*wts;
				tsum = sum(wts);
				md <- sum(pmax(md/ud,-md/ld),na.rm=TRUE)/tsum;
				return (md)
			}
			metric <- apply(datasubset,1,manDistance,tem,qld,qud,fwts);
	  },
		{
			corDistance <- function (x,template,method) {md <- 3.0*(1.0-cor(x,template,method=method,use="pairwise.complete.obs")); return (md)}
			if (class(template)[1]=="matrix")
			{
				metric <- numeric(nrow(datasubset));
				swts <- 0;
				for (i in 1:length(theQuant))
				{
					tem <- template[i,];
					wts <- theQuant[i];
					if (wts > 0.5)
					{
						wts <- 1.0-wts;
					}
					metric <- metric + wts*(apply(datasubset,1,corDistance,template=tem,method=method));
					swts <- swts + wts;
				}
				metric <- metric/swts;
			}
			else
			{
				tem <- template;
				metric <- apply(datasubset,1,corDistance,template=tem,method=method);
			}
		}
	)
	names(metric) <- rownames(data);
	attr(metric,"ld") <- ld;
	attr(metric,"ud") <- ud;
	attr(metric,"qld") <- qld;
	attr(metric,"qud") <- qud;
	metric[is.na(metric)] <- 1.0e10;
	
	result <- metric
	return (result);
}
