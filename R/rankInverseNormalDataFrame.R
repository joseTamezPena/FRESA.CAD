rankInverseNormalDataFrame <-
function(variableList,data,referenceframe,strata=NA) 
{
  
	if (class(variableList)=="character")
	{
		colnamesList <- variableList;
	}
	else
	{
		colnamesList <- as.vector(variableList[,1]);
	}
	size = length(colnamesList);

	if (!is.na(strata)) 
	{
		maxStrata = max(referenceframe[,strata],na.rm = TRUE);
		minStrata = min(referenceframe[,strata],na.rm = TRUE);
	}
	else 
	{
		maxStrata=1;
		minStrata=1;
	}
#	cat ("Min Strata:",minStrata,"Max Strata:",maxStrata,"\n");
	created=0;
	for (sta in minStrata:maxStrata)
	{
		if (!is.na(strata))
		{
			stracondition = paste (strata,paste('==',sta));
			strastatement = paste ("subset(referenceframe,",paste(stracondition,")"));
#			cat ("Strata:",stracondition,"\n");
			cstrataref <- eval(parse(text=strastatement));
			strastatement = paste ("subset(data,",paste(stracondition,")"));
			cstrata <- eval(parse(text=strastatement));
#			cat ("Rows:",nrow(cstrataref),"Rows 2",nrow(cstrata)," \n");
		}
		else
		{
			cstrataref = referenceframe;
			cstrata = data;
		}
		if ((nrow(cstrata)>1) && ( nrow(cstrataref)>1))
		{

			SortedCtr = cstrataref;
			nrowsCtr =  nrow(cstrataref);
			nrows = nrow(cstrata);
			InverseFrame = cstrata;
			idxs <- as.integer(0.01*nrowsCtr); #for outlies smoothing
			for (i in 1:size)
			{ 
				if (length(table(cstrata[,colnamesList[i]]))>2)
				{
					SortedCtr[,colnamesList[i]]<-SortedCtr[order(cstrataref[,colnamesList[i]]),colnamesList[i]];
					minvalue <- SortedCtr[1,colnamesList[i]];
					maxvalue <- SortedCtr[nrowsCtr,colnamesList[i]];
#					cat(" Variable: ",colnamesList[i],"Min: ",minvalue," Max: ",maxvalue)      
					if (idxs>0)
					{
						for (n in idxs:1) SortedCtr[n,colnamesList[i]] <- 0.5*(SortedCtr[n,colnamesList[i]]+SortedCtr[n+1,colnamesList[i]]);
						for (n in (nrowsCtr-idxs):nrowsCtr) SortedCtr[n,colnamesList[i]] <- 0.5*(SortedCtr[n,colnamesList[i]]+SortedCtr[n-1,colnamesList[i]]);
						minvalue <- SortedCtr[1,colnamesList[i]];
						maxvalue <- SortedCtr[nrowsCtr,colnamesList[i]];
#						cat("-> NMin: ",minvalue," NMax: ",maxvalue)      
					}
#					cat("\n")      
					InverseFrame[,colnamesList[i]]<-.Call("rankInverseNormalCpp",nrows,cstrata[,colnamesList[i]],minvalue,maxvalue,SortedCtr[,colnamesList[i]]);
				}
			}
		}
		if (created == 1) 
		{
			zRankInverseFrame <- rbind(zRankInverseFrame,InverseFrame);
		}
		else
		{
			created = 1;
			zRankInverseFrame = InverseFrame;
		}
	}
	return (zRankInverseFrame);
}
