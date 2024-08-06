\name{FRESA.CAD-package}
\alias{FRESA.CAD-package}
\alias{FRESA.CAD}
\docType{package}
\title{FeatuRE Selection Algorithms for Computer-Aided Diagnosis (FRESA.CAD)}
\description{
Contains a set of utilities for building and testing formula-based models for Computer Aided Diagnosis/prognosis applications via feature selection.
Bootstrapped Stage Wise Model Selection (B:SWiMS) controls the false selection (FS) for linear, logistic, or Cox proportional hazards regression models.
Utilities include functions for: univariate/longitudinal analysis, data conditioning (i.e. covariate adjustment and normalization), model validation and visualization.
}


\details{
    \tabular{ll}{
        Package: \tab FRESA.CAD\cr
        Type: \tab Package\cr
        Version: \tab 3.4.9\cr
        Date: \tab 2024-08-01\cr
        License: \tab LGPL (>= 2)\cr
    }
    Purpose: The design of diagnostic or prognostic multivariate models via the selection of significantly discriminant features.
    The models are selected via the bootstrapped step-wise selection of model features that offer a significant improvement in subject classification/error.
    The false selection control is achieved by train-test partitions, where train sets are used to select variables and test sets used to evaluate model performance.
    Variables that do not improve subject classification/error on the blind test are not included in the models.
    The main function of this package is the selection and cross-validation of diagnostic/prognostic linear, logistic, or Cox proportional hazards regression model constructed from a large set of candidate features.
    The variable selection may start by conditioning all variables via a covariate-adjustment and a \emph{z}-inverse-rank-transformation. 
    In order to integrate features with partial discriminant power, the package can be used to categorize the continuous variables and rank their discriminant power.
    Once ranked, each feature is bootstrap-tested in a multivariate model, and its blind performance is evaluated.
    Variables with a statistical significant improvement in classification/error are stored and finally inserted into the final model according to their relative store frequency. 
    A cross-validation procedure may be used to diagnose the amount of model shrinkage produced by the selection scheme.
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{
Jose Gerardo Tamez-Pena, Antonio Martinez-Torteya, Israel Alanis and Jorge Orozco 
Maintainer: <jose.tamezpena@tec.mx>
}
\examples{
    \dontrun{
    ### Fresa Package Examples ####
    library("epiR")
    library("FRESA.CAD")
    library(network)
    library(GGally)
    library("e1071")


    # Start the graphics device driver to save all plots in a pdf format
    pdf(file = "Fresa.Package.Example.pdf",width = 8, height = 6)


    # Get the stage C prostate cancer data from the rpart package

    data(stagec,package = "rpart")
    options(na.action = 'na.pass')
    dataCancer <- cbind(pgstat = stagec$pgstat,
                        pgtime = stagec$pgtime,
                        as.data.frame(model.matrix(Surv(pgtime,pgstat) ~ .,stagec))[-1])

    #Impute missing values
    dataCancerImputed <- nearestNeighborImpute(dataCancer)

    # Remove the incomplete cases
    dataCancer <- dataCancer[complete.cases(dataCancer),]


    # Load a pre-stablished data frame with the names and descriptions of all variables
    data(cancerVarNames)
    # the Heat Map
    hm <- heatMaps(cancerVarNames,varRank=NULL,Outcome="pgstat",
                   data=dataCancer,title="Heat Map",hCluster=FALSE
                   ,prediction=NULL,Scale=TRUE,
                   theFiveColors=c("blue","cyan","black","yellow","red"),
                   outcomeColors = 
                     c("blue","lightgreen","yellow","orangered","red"),
                   transpose=FALSE,cexRow=0.50,cexCol=0.80,srtCol=35)

    # The univariate analysis
    UniRankFeaturesRaw <- univariateRankVariables(variableList = cancerVarNames,
                                                  formula = "pgstat ~ 1+pgtime",
                                                  Outcome = "pgstat",
                                                  data = dataCancer, 
                                                  categorizationType = "Raw", 
                                                  type = "LOGIT", 
                                                  rankingTest = "zIDI",
                                                  description = "Description",
                                                  uniType="Binary")

    print(UniRankFeaturesRaw)

    # A simple BSIWMS Model

    BSWiMSModel <- BSWiMS.model(formula = Surv(pgtime, pgstat) ~ 1, dataCancerImputed)

    # The Log-Rank Analysis using survdiff

    lrsurvdiff <- survdiff(Surv(pgtime,pgstat)~
                  BSWiMSModel$BSWiMS.model$back.model$linear.predictors > 0,
                  data=dataCancerImputed)

    # The Log-Rank Analysis EmpiricalSurvDiff and permutations of the null Chi distribution
    lrp <- EmpiricalSurvDiff(dataCancerImputed$pgtime,dataCancerImputed$pgstat,
                             BSWiMSModel$BSWiMS.model$back.model$linear.predictors > 0,
                             type="Chi",plots=TRUE,samples=10000)

    # The Log-Rank Analysis EmpiricalSurvDiff and permutations of the null SLR distribution
    lrp <- EmpiricalSurvDiff(dataCancerImputed$pgtime,dataCancerImputed$pgstat,
                             BSWiMSModel$BSWiMS.model$back.model$linear.predictors > 0,
                             type="SLR",plots=TRUE,samples=10000)

    # The Log-Rank Analysis EmpiricalSurvDiff and bootstrapping the SLR distribution
    lrp <- EmpiricalSurvDiff(dataCancerImputed$pgtime,dataCancerImputed$pgstat,
                             BSWiMSModel$BSWiMS.model$back.model$linear.predictors > 0,
                             computeDist=TRUE,plots=TRUE)

    #The performance of the final model using the summary function
    sm <- summary(BSWiMSModel$BSWiMS.model$back.model)
    print(sm$coefficients)
    pv <- plot(sm$bootstrap)

    # The equivalent model
    eq <- reportEquivalentVariables(BSWiMSModel$BSWiMS.model$back.model,data=dataCancer,
                                    variableList=cancerVarNames,Outcome = "pgstat",
                                    timeOutcome="pgtime",
                                    type = "COX");

    print(eq$equivalentMatrix)

    #The list of all models of the bootstrap forward selection 
    print(BSWiMSModel$forward.selection.list)

    #With FRESA.CAD we can do a leave-one-out using the list of models
    pm <- ensemblePredict(BSWiMSModel$forward.selection.list,
                          dataCancer,predictType = "linear",type="LOGIT",Outcome="pgstat")

    #Ploting the ROC with 95% CI of the models
    pm <- plotModels.ROC(cbind(dataCancer$pgstat,
                               pm$ensemblePredict),main=("LOO Forward Selection Median Predict"))

    #The plotModels.ROC provides the diagnosis confusion matrix.
    summary(epi.tests(pm$predictionTable))



    #FRESA.CAD can be used to create a bagged model using the forward selection formulas
    bagging <- baggedModel(BSWiMSModel$forward.selection.list,dataCancer,useFreq=32)
    pm <- predict(bagging$bagged.model)
    pm <- plotModels.ROC(cbind(dataCancer$pgstat,pm),main=("Bagged"))

    #Let's check the performance of the model 
    sm <- summary(bagging$bagged.model)
    print(sm$coefficients)

    #Using bootstrapping object I can check the Jaccard Index
    print(bagging$Jaccard.SM)

    #Ploting the evolution of the coefficient value
    plot(bagging$coefEvolution$grade,main="Evolution of grade")


    gplots::heatmap.2(bagging$formulaNetwork,trace="none",
                      mar=c(10,10),main="eB:SWIMS Formula Network")
    barplot(bagging$frequencyTable,las = 2,cex.axis=1.0,
            cex.names=0.75,main="Feature Frequency")
    n <- network::network(bagging$formulaNetwork, directed = FALSE,
                          ignore.eval = FALSE,names.eval = "weights")
    ggnet2(n, label = TRUE, size = "degree",size.cut = 3,size.min = 1, 
           mode = "circle",edge.label = "weights",edge.label.size=4)


    # Get a Cox proportional hazards model using:
    # - The default parameters

    mdCOXs <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,data = dataCancer)
    sm <- summary(mdCOXs$BSWiMS.model)
    print(sm$coefficients)

    # The model with singificant improvement in the residual error
    mdCOXs <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
                          data = dataCancer,OptType = "Residual" )
    sm <- summary(mdCOXs$BSWiMS.model)
    print(sm$coefficients)

    # Get a Cox proportional hazards model using second order models:
    mdCOX <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
                         data = dataCancer,categorizationType="RawRaw")
    sm <- summary(mdCOX$BSWiMS.model)
    print(sm$coefficients)

    namesc <- names(mdCOX$BSWiMS.model$coefficients)[-1]
    hm <- heatMaps(mdCOX$univariateAnalysis[namesc,],varRank=NULL,
                   Outcome="pgstat",data=dataCancer,
                   title="Heat Map",hCluster=FALSE,prediction=NULL,Scale=TRUE,
                   theFiveColors=c("blue","cyan","black","yellow","red"),
                   outcomeColors = c("blue","lightgreen","yellow","orangered","red"),
                   transpose=FALSE,cexRow=0.50,cexCol=0.80,srtCol=35)

    # The LOO estimation
    pm <- ensemblePredict(mdCOX$BSWiMS.models$formula.list,dataCancer,
                          predictType = "linear",type="LOGIT",Outcome="pgstat")
    pm <- plotModels.ROC(cbind(dataCancer$pgstat,pm$ensemblePredict),main=("LOO Median Predict"))
    #Let us check the diagnosis performance
    summary(epi.tests(pm$predictionTable))

    # Get a Logistic model using FRESA.Model
    # - The default parameters
    dataCancer2 <-dataCancer 
    dataCancer2$pgtime <-NULL
    mdLOGIT <- FRESA.Model(formula = pgstat ~ 1,data = dataCancer2)
    if (!is.null(mdLOGIT$bootstrappedModel)) pv <- plot(mdLOGIT$bootstrappedModel)
    sm <- summary(mdLOGIT$BSWiMS.model)
    print(sm$coefficients)


    ## FRESA.Model with Cross Validation and Recursive Partitioning and Regression Trees


    md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,data = dataCancer,
                      CVfolds = 10,repeats = 5,equivalent = TRUE,usrFitFun=rpart::rpart)

    colnames(md$cvObject$Models.testPrediction)

    pm <- plotModels.ROC(md$cvObject$LASSO.testPredictions,theCVfolds=10,main="CV LASSO",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$KNN.testPrediction,theCVfolds=10,main="KNN",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="Prediction",main="B:SWiMS Bagging",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="Ensemble.B.SWiMS"
                         ,main="Forward Selection Median Ensemble",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="Ensemble.Forward",main="Forward Selection Bagging",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="eB.SWiMS",main="Equivalent Model",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="Forward.Selection.Bagged",main="The Forward Bagging",cex=0.90)

    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=20,
                         predictor="usrFitFunction",
                         main="Recursive Partitioning and Regression Trees",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=20,
                         predictor="usrFitFunction_Sel",
                         main="Recursive Partitioning and Regression Trees with FS",cex=0.90)


    ## FRESA.Model with Cross Validation, LOGISTIC and Support Vector Machine


    md <- FRESA.Model(formula = pgstat ~ 1,data = dataCancer2,
                      CVfolds = 10,repeats = 5,equivalent = TRUE,usrFitFun=svm)

    pm <- plotModels.ROC(md$cvObject$LASSO.testPredictions,theCVfolds=10,main="CV LASSO",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$KNN.testPrediction,theCVfolds=10,main="KNN",cex=0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="Prediction",main="B:SWiMS Bagging",cex=0.90)

    md$cvObject$Models.testPrediction[,"usrFitFunction"] <- 
                      md$cvObject$Models.testPrediction[,"usrFitFunction"] - 0.5
    md$cvObject$Models.testPrediction[,"usrFitFunction_Sel"] <- 
                      md$cvObject$Models.testPrediction[,"usrFitFunction_Sel"] - 0.5
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="usrFitFunction",
                         main="SVM",cex = 0.90)
    pm <- plotModels.ROC(md$cvObject$Models.testPrediction,theCVfolds=10,
                         predictor="usrFitFunction_Sel",
                         main="SVM with FS",cex = 0.90)


    # Shut down the graphics device driver
    dev.off()

   }
}
\keyword{package}
