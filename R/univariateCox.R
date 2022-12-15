univariate_cox <- function(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", limit=0,...)
{
  # Outcome <- "pgstat"
  baseformula <- as.character(Outcome);
  featuresOnSurvival <- baseformula[2]
  featuresOnSurvival <- gsub(" ", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("Surv\\(", "", featuresOnSurvival)
  featuresOnSurvival <- gsub("\\)", "", featuresOnSurvival)
  featuresOnSurvivalObject <- strsplit(featuresOnSurvival, ",")
  
  varlist <-colnames(data);
  varlist <- varlist[!varlist %in% featuresOnSurvivalObject[[1]]];
  
  univ <- MultUnivariateCox(varlist,data,Outcome) 
  univ <- univ[order(univ)];
  unitPvalues <- p.adjust(univ,adjustMethod);
  unitPvalues <- unitPvalues[order(unitPvalues)];
  top <- unitPvalues[1];
  unitPvalues <- unitPvalues[unitPvalues <= pvalue];
  #	print(unitPvalues)
  if (length(unitPvalues) > 1)
  {
    #		unitPvalues <- unitPvalues[correlated_Remove(data,names(unitPvalues))];
    unitPvalues <- try(correlated_RemoveToLimit(data,unitPvalues,limit,...));
    if (!inherits(unitPvalues, "try-error"))
		{
      return(unitPvalues);
    }
    else{
      return(top);
    }
  }
  else
  {
    return(top);
  }
  
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
  
  univ_formulas <- sapply(varlist, function(x) as.formula(paste('Surv(time,status) ~ ', x)))
  
  univ_models <- lapply(univ_formulas, function(x){survival::coxph(x, data = data)})
  
  #signif(summary(univ_models$age)$wald["pvalue"], digits=2)
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           return(p.value)
                         })
  pvalues <- as.matrix(univ_results)
  names(pvalues)<-rownames(pvalues)
  
  return(pvalues)
}
