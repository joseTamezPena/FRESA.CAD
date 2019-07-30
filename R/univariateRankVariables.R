univariateRankVariables <-
function(variableList,formula,Outcome,data,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail","Tail","RawRaw"),type=c("LOGIT","LM","COX"),rankingTest=c("zIDI","zNRI","IDI","NRI","NeRI","Ztest","AUC","CStat","Kendall"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,description=".",uniType=c("Binary","Regression"),FullAnalysis=TRUE,acovariates=NULL,timeOutcome=NULL) 
{
	orderframe <- uniRankVar(variableList,formula,Outcome,data,categorizationType,type,rankingTest,cateGroups,raw.dataFrame,description,uniType,FullAnalysis,acovariates,timeOutcome)$orderframe;
	return (orderframe)
}