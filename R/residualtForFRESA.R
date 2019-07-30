residualForFRESA <-
function (object,testData,Outcome,eta=0.05) 
{
#Set eta to zero to get prediction residuals for cox models. eta to 1 get the Martingale residuals
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
#	cat(classlen," ",cobj,"\n");
	if ((cobj!="co") && !is.null(object$family[1])) cobj = "fi";
	switch(cobj,
		co =
		{
			if (classlen==1)
			{
				lpp <- 1.0/(1.0+exp(-predict(object,testData,'lp',na.action=na.omit)));
				s <- is.na(lpp);
				if (any(s)) 
				{
					lpp[s] = 1.0e10; # set to a large residual
				}
				lppres <- lpp - testData[,Outcome];
				matingale <- testData[,Outcome]-predict(object,testData,'expected',na.action=na.omit);
				s <- is.na(matingale);
				if (any(s))
				{
					matingale[s] = 1.0e10; # set to a large residual
				}
				out <- (1-eta)*lppres - eta*matingale;			
			}
			else
			{
				out <- (1 - 2*testData[,Outcome]);
			}
		},
		lm =
		{
			out <- predict.fitFRESA(object,testData,'linear') - testData[,Outcome];
#			out <- predict(object,testData) - testData[,Outcome];
		},
		fi =
		{
			if (object$type == "LM")
			{
				out <- predict.fitFRESA(object,testData,'linear') - testData[,Outcome];
			}
			else
			{
				out <- predict.fitFRESA(object,testData,'prob') - testData[,Outcome];
			}
		},
		tr =
		{
			if (object$type == "LM")
			{
				out <- predict.fitFRESA(object,testData,'linear') - testData[,Outcome];
			}
			else
			{
				out <- predict.fitFRESA(object,testData,'prob') - testData[,Outcome];
			}
		},
		{
			if (object$family[1] == "binomial")
			{
				out <- predict.fitFRESA(object,testData,'prob') - testData[,Outcome];
			}
			else
			{
				out <- predict.fitFRESA(object,testData,'linear') - testData[,Outcome];
			}
		}
	)	
	s <- is.na(out);
	if (any(s)) 
	{
		warning("Warning NA predictFor NeRIs \n");
		out[s] <- 1.23456789e10; # Set a large residual
	}
    return (out)
}
