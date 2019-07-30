orderFeatures <- function(modelFormulas,univariate=NULL,baseForm="",useFreq=TRUE,n_bootstrap=1)
{
	VarFrequencyTable <- NULL;
	theterms <- list();
	features <- vector();
	forder <- numeric();
	fcount <- numeric();
	inserted <- 1;
	loops <- length(modelFormulas);
	formulaLoops <- 1;
	if (is.numeric(useFreq))
	{
		formulaLoops <- useFreq;
		useFreq <- TRUE;
	}
	iscompletef <- (gregexpr(pattern ='~',modelFormulas[1])[1]>0)
	if (baseForm=="")
	{
		if (!iscompletef)
		{
			baseForm <- "Dummy ~ 1";
		}
		else
		{
			baseForm <- paste(unlist(strsplit(as.character(modelFormulas[1]),"[~]"))[1],"~ 1");
		}
	}
#	print(baseForm);
	removed <- numeric();
	nformula <- 0;
#	cat("orderFeatures\n");
	for (n in 1:loops)
	{
#		print(modelFormulas[n]);
		if (useFreq) 
		{
			inserted <- 1;
		}
		if (grepl("=-=",as.character(modelFormulas[n]),fixed = TRUE))
		{
			inserted <- 1;
			modelFormulas[n] <- baseForm;
			removed <- c(removed,n);
		}
		else
		{
			nformula <- nformula + 1
			if (iscompletef)	
			{
				modelFormulas[n] <- unlist(strsplit(as.character(modelFormulas[n]),"[~]"))[2];
			}
			if (nchar(modelFormulas[n])>0)
			{
				modelFormulas[n] <- paste(baseForm,"+",modelFormulas[n]);
			}
			else
			{
				modelFormulas[n] <- baseForm;
			}
	#		print(modelFormulas[n]);
			termList <- str_replace_all(attr(terms(formula(modelFormulas[n])),"term.labels"),":","\\*");
	#		print(termList);
			if (length(termList)>0)
			{
				features <- append(features,termList);
				rtermlist <- termList[!(termList %in% names(forder))]
				sforder <- numeric(length(rtermlist));
				names(sforder) <- rtermlist
				forder <- append(forder,sforder)
				fcount <- append(fcount,sforder)
				for (f in 1:length(termList))
				{
					forder[termList[f]] <- forder[termList[f]]+inserted;
					fcount[termList[f]] <- fcount[termList[f]]+1;
					inserted <- inserted+1;
				}
			}
			theterms[[nformula]] <- termList;
		}
	}
	numberofBreaks <- length(removed);
	if (length(removed) > 0)
	{
#		print(removed)
#		print(length(modelFormulas));
		modelFormulas <- modelFormulas[-removed];
		loops <- length(modelFormulas);
#		print(modelFormulas);
	}
	if (length(features)>0)
	{
		VarFrequencyTable <- table(features);
#		print(VarFrequencyTable);
	}
	if (length(VarFrequencyTable)>1) 
	{
		forder=forder/fcount;
		oVarFrequencyTable <- VarFrequencyTable;
		if (formulaLoops > 1 )
		{
			cycles <- loops/formulaLoops;
			VarFrequencyTable <- 0*VarFrequencyTable;
			vnames <- rownames(VarFrequencyTable);
			for (lo in 1:cycles)
			{
				features <- vector();
				offs <- (lo-1)*formulaLoops;
				for (m in 1:formulaLoops)
				{
					features <- append(features,theterms[[m+offs]]);
				}
				tmpVarFrequencyTable <- table(features);
				if (length(tmpVarFrequencyTable)>0)
				{
					thetnames <- vnames %in% rownames(tmpVarFrequencyTable)
					if (!is.null(thetnames))
					{
						otablenames <- rownames(VarFrequencyTable[thetnames]);
						if (!is.null(otablenames))
						{
							freq1 <- as.numeric(VarFrequencyTable[thetnames]);
							freq2 <- as.numeric(tmpVarFrequencyTable[otablenames]);
							VarFrequencyTable[thetnames] <- pmax(freq1,freq2);
						}
					}
				}
			}
			fmax <- max(oVarFrequencyTable);
			VarFrequencyTable <- VarFrequencyTable+(oVarFrequencyTable-1.0)/fmax;
		}
		if (!is.null(univariate))
		{
			if (!is.null(univariate$ZUni))
			{
				VarFrequencyTable <- VarFrequencyTable[order(-univariate[rownames(VarFrequencyTable),"ZUni"])];		
			}				
		}
		VarFrequencyTable <- VarFrequencyTable[order(forder[rownames(VarFrequencyTable)])];
		if (useFreq)
		{
			VarFrequencyTable <- VarFrequencyTable[order(-VarFrequencyTable)]
		}
#		print(VarFrequencyTable);
#		print(forder);
	}
#		cat(" End orderFeatures\n");

	result <- list(VarFrequencyTable=VarFrequencyTable,modelFormulas=modelFormulas,forder=forder,theterms=theterms,features=features,numberofBreaks=numberofBreaks)
	return (result);
}
