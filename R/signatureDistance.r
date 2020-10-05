signatureDistance <- 
function (template, data=NULL, method = c("pearson","spearman","kendall","RSS","MAN"))
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
	
	if (class(template)[1] == "list")
	{
		theQuant <- template$quant;
		template <- template$template;
	}

	wvalues <- 1.0/abs(qnorm(theQuant));
	
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
	if (class(template)[1] == "matrix")
	{
		tem <- template[medianv,];
		
#		cat("median:")
#		print(tem)
		wts <- 0.0;
		ld <- numeric(length(tem));
		for (i in 1:(medianv - 1))
		{
			wts <- wts + theQuant[i];
			ld <- ld + theQuant[i]*wvalues[i]*(tem - template[i,]);
		}
		ld <- ld/wts;
		qld <- (tem - template[medianv - 1,])*wvalues[medianv - 1];

		wts <- 0;
		ud <- numeric(length(tem));
		for (i in (medianv + 1):length(wvalues))
		{
			wts <- wts + (1.0-theQuant[i]);
			ud <- ud + (1.0-theQuant[i])*wvalues[i]*(template[i,] - tem);
		}
		ud <- ud/wts;
		qud <- (template[medianv + 1,] - tem)*wvalues[medianv + 1];

		ld[ld == 0] <- 0.33*ud[ld == 0];
		ld[ld == 0] <- 0.33;
		qld[qld == 0] <- ld[qld == 0];

		ud[ud == 0] <- 0.33*ld[ud == 0];
		ud[ud == 0] <- 0.33;
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
		ld[ld == 0] <- 0.33;

		ud <- ld;
		qld <- IQR(template)/abs(qnorm(0.25))/2;
		qld[qld == 0] <- ld[qld == 0];
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
			metric <- apply(datasubset,1,RSSDistance,tem,ld,ud)/sqrt(ncol(template));
		},
		MAN = 
		{ 
			manDistance <- function (x,template,ld,ud) 
			{
				md <- x-template
				md <- sum(pmax(md/ud,-md/ld),na.rm=TRUE);
				return (md)
			}
			metric <- apply(datasubset,1,manDistance,tem,qld,qud)/ncol(template);
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
	metric[is.na(metric)] <- 1.0e10;
	
	result <- metric
	return (result);
}
