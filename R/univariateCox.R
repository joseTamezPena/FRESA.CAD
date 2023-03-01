univariate_cox <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", limit=0,...,n = 0)
{
  # Outcome <- "pgstat"
  baseformula <- as.character(Outcome);
  featuresOnSurvival <- baseformula[2]
  featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
  featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
  
  varlist <- colnames(data);
  varlist <- varlist[!varlist %in% featuresOnSurvivalObject[[1]]];
  
  # univ <- MultUnivariateCox(varlist,data,Outcome)
  # univ <- univ[order(univ)];
  # unitPvalues <- p.adjust(univ,adjustMethod);
  # top <- unitPvalues[1];
  # unitPvalues <- unitPvalues[unitPvalues <= pvalue];
  # if (length(unitPvalues) > 1)
  # {
    # unitPvalues <- try(correlated_RemoveToLimit(data,unitPvalues,limit,...));
    # if (!inherits(unitPvalues, "try-error"))
		# {
      # return(unitPvalues);
    # }
    # else{
      # return(top);
    # }
  # }
  # else
  # {
    # return(top);
  # }
  
    unitPvalues <- MultUnivariateCox(varlist,data,Outcome)
    unitPvalues <- unitPvalues[order(unitPvalues)];
    unadjusted <- unitPvalues;
	unitPvalues <- p.adjust(unitPvalues,adjustMethod,n=max(n,length(unitPvalues)));
	top <- unitPvalues[1];
	if ((top < 0.45) && (unadjusted[1] < 0.05))
	{
		top <- unitPvalues[unitPvalues <= 1.01*top];
	}
	unitPvalues <- unitPvalues[unitPvalues <= pvalue];
	if (length(unitPvalues) > 1) 
	{
		unitPvalues <- correlated_RemoveToLimit(data,unitPvalues,limit,...);
	}
	else 
	{
		unitPvalues <- top;
	}
	attr(unitPvalues,"Unadjusted") <- unadjusted;
    return(unitPvalues);
}

MultUnivariateCox <- function(varlist = NULL,data = NULL,Outcome = NULL,...)
{
  baseformula <- as.character(Outcome);
  featuresOnSurvival <- baseformula[2]
  featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
  featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
  time <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][1]]))
  status <- as.numeric(unlist(data[featuresOnSurvivalObject[[1]][2]]))
  srv <- survival::Surv(time,status)
  
#  univ_formulas <- sapply(varlist, function(x) as.formula(paste('Surv(time,status) ~ ', x)))
  
#  univ_models <- lapply(univ_formulas, function(x){survival::coxph(x, data = data)})
  
  #signif(summary(univ_models$age)$wald["pvalue"], digits=2)
  # Extract data 
#  univ_results <- lapply(univ_models,
#                         function(x){ 
#                           x <- summary(x)
#                           p.value<-signif(x$wald["pvalue"], digits=2)
#                           return(p.value)
#                         })
#  pvalues <- as.matrix(univ_results)
#  names(pvalues)<-rownames(pvalues)

  coxpvalue <- function(x)
  {
    s <- summary(survival::coxph(srv ~  x))
    return(s$wald["pvalue"])
  }
  pvalues <- apply(data[,varlist],2,coxpvalue);
  names(pvalues) <- varlist;
  
  return(pvalues)
}
