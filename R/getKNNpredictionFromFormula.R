getKNNpredictionFromFormula <-
function (model.formula,trainData,testData,Outcome="CLASS",nk=3) 
{

if (!requireNamespace("class", quietly = TRUE)) {
   install.packages("class", dependencies = TRUE)
} 

#	print(as.character(model.formula))
#	cat(Outcome,"\n")
	varlist <- NULL;
	if (class(model.formula) == "character")
	{
		if (length(model.formula)==1)
		{
			model.formula <- formula(model.formula);
		}
		else
		{
			varlist <- model.formula;
		}
	}
	if (class(model.formula) == "formula")
	{
		varlist <- as.character(model.formula)
		varlist <- all.vars(formula(paste(Outcome,"~",varlist[3])));
		names(varlist) <- varlist;
		varlist <- varlist[!(varlist %in% Outcome)];
	}

	nrows <- nrow(trainData);
	ncolms <- length(varlist);
	if (ncolms>0)
	{
#		cat("Rows:",nrows," Columns:",ncolms,"\n");

		trainframe <- 	as.data.frame(trainData[,varlist]+matrix(rnorm(nrows*ncolms,0,1e-10),nrow=nrows))
		testframe <- as.data.frame(testData[,varlist])
	}
	else
	{
		cat("Warning: no predictor for KNN\n");

		trainframe <- 	as.data.frame(matrix(rnorm(nrows,0,1e-10),nrow=nrows))
		tnrows <- nrow(testData);
		testframe <- as.data.frame(matrix(rnorm(tnrows,0,1e-10),nrow=tnrows))
	}
	knnclass <- try(class::knn(trainframe,testframe,factor(trainData[,Outcome]),nk,prob=TRUE))
	if (!inherits(knnclass, "try-error"))
	{
		prop <- attributes(knnclass);
		binProp <-abs(prop$prob-1*(knnclass=="0"))
	}
	else
	{
		prop <- NULL
		binProp <- vector(nrow(testData));
	}
	result <- list(prediction=knnclass,
	prob=prop,
	binProb=binProp,
	featureList=varlist)
    return (result)
}
