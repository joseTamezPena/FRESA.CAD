signatureDistance <- 
function (template, data=NULL, method = c("pearson","spearman","kendall","RSS","MAN"))
{

#given the template: mean,median,sample, etc....;signatureDistance it will return the distance between the template to each row of the dataframe
#the template is a named numeric vector
#the data is a colnamed data frame
#methods:
# RSS: Root Sum Square
# MAN: Manhattan distance
# pearson: 1-Pearson correlation coefficient
# spearman: 1-spearman correlation coefficient
# kendall: 1-kendall correlation coefficient

	method <- match.arg(method)
	theQuant <- c(0.025,0.16,0.25,0.5,0.75,0.84,0.975);
	
	if (class(template) == "list")
	{
		theQuant <- template$quant;
		template <- template$template;
	}

	wvalues <- 1.0/abs(qnorm(theQuant));
	
	if (class(template)=="matrix")
	{
		vnames <- colnames(template);
	}
	else
	{
		vnames <- names(template);
	}
	datasubset <- as.matrix(data[,vnames]);
	
	medianv <- as.integer((length(theQuant) + 1)/2);
	if (class(template) == "matrix")
	{
		tem <- template[medianv,];
		wts <- 0.0;
		ld <- numeric(length(tem));
		for (i in 1:(medianv - 1))
		{
			wts <- wts + 1.0;
			ld <- ld + wvalues[i]*(tem - template[i,]);
		}
		ld <- ld/wts;
		mv <- min(ld[ld > 0]);
		if (is.numeric(mv))	{ ld[ld == 0] <- mv }
		else { ld[ld == 0] <- 1.0; }
		qld <- tem - template[medianv - 1,];
		qld[qld == 0] <- ld[qld == 0];

		wts <- 0;
		ud <- numeric(length(tem));
		for (i in (medianv + 1):length(wvalues))
		{
			wts <- wts + 1.0;
			ud <- ud + wvalues[i]*(template[i,] - tem);
		}
		ud <- ud/wts;
		mv <- min(ud[ud > 0]);
		if (is.numeric(mv))	{ ud[ud == 0] <- mv; }
		else { ud[ud == 0] <- 1.0; }
		qud <- template[medianv + 1,] - tem;
		qud[qud == 0] <- ud[qud == 0];
	}
	else
	{
		tem <- template;
		ld <- sd(template);
		ud <- ld;
		qld <- IQR(template)/2;
		qud <- qld;
	}
	switch(method, 
		RSS = 
		{ 
			RSSDistance <- function (x,template,ld,ud) 
			{
				md <- x-template
				md <- sqrt(sum(pmax(md/ud,-md/ld)^2,na.rm=TRUE));
				return (md)
			}
			metric <- apply(datasubset,1,RSSDistance,tem,ld,ud);
		},
		MAN = 
		{ 
			manDistance <- function (x,template,ld,ud) 
			{
				md <- x-template
				md <- sum(pmax(md/ud,-md/ld),na.rm=TRUE);
				return (md)
			}
			metric <- apply(datasubset,1,manDistance,tem,qld,qud);
	  },
		{
			corDistance <- function (x,template,method) {md <- 1.0-cor(x,template,method=method,use="pairwise.complete.obs"); return (md)}
			if (class(template)=="matrix")
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
	metric[is.na(metric)] <- 1.0e10;
	
	result <- metric
	return (result);
}
