\name{randomCV}
\alias{randomCV}
\title{Cross Validation of Prediction Models}
\description{
	The data set will be divided into a random train set and a test sets. The train set will be modeled by the user provided fitting 
	method. Each fitting method must have a prediction function that will be used to predict the outcome of the test set.
}

\usage{
	randomCV(theData = NULL,
             theOutcome = "Class",
             fittingFunction=NULL, 
             trainFraction = 0.5,
             repetitions = 100,
             trainSampleSets=NULL,
             featureSelectionFunction=NULL,
             featureSelection.control=NULL,
             asFactor=FALSE,
             addNoise=FALSE,
             classSamplingType=c("Augmented","NoAugmented","Proportional","Balanced"),
             ...
             )
}
\arguments{
	\item{theData}{
	The data-frame for cross-validation
	}
	\item{theOutcome}{
	The name of the outcome
	}
	\item{fittingFunction}{
	The fitting function used to model the data
	}
	\item{trainFraction}{
	The percentage of the data to be used for training
	}
	\item{repetitions}{
	The number of times that the CV process will be repeated
	}
	\item{trainSampleSets}{
	A set of train samples
	}
	\item{featureSelectionFunction}{
	The feature selection function to be used to filter out irrelevant features
	}
	\item{featureSelection.control}{
	The parameters to control the feature selection function
	}
	\item{asFactor}{
	Set theOutcome as factor
	}
	\item{addNoise}{
	if TRUE will add 0.1% of random noise to the predictors
	}
	\item{classSamplingType}{
	if "Proportional": proportional to the data classes. 
	"Augmented": Augment samples to balance training class
	"NoAugmented": Avoid Data augmentation hece training may not be balanced
	"Balanced": All class in training set have the same samples
	}
	\item{...}{
	Parameters to be passed to the fitting function
	}
}

\value{
	\item{testPredictions}{
	All the predicted outcomes. Is a data matrix with three columns c("Outcome","Model","Prediction"). Each row has a prediction for a given test subject
	}
	\item{trainPredictions}{
	All the predicted outcomes in the train data set. Is a data matrix with three columns c("Outcome","Model","Prediction"). Each row has a prediction for a given test subject
	}
	\item{medianTest}{
	The median of the test prediction for each subject
	}
	\item{medianTrain}{
	The median of the prediction for each train subject
	}
	\item{boxstaTest}{
	The statistics of the boxplot for test data
	}
	\item{boxstaTrain}{
	The statistics of the boxplot for train data
	}
	\item{trainSamplesSets}{
	The id of the subjects used for training
	}
	\item{selectedFeaturesSet}{
	A list with all the features used at each training cycle
	}
	\item{featureFrequency}{
	A order table object that describes how many times a feature was selected.
	}
	\item{jaccard}{
	The jaccard index of the features as well as the average number of features used for prediction
	}
	\item{theTimes}{
	The CPU time analysis
	}
}
\author{Jose G. Tamez-Pena}
\examples{
	\dontrun{

        ### Cross Validation Example ####
        # Start the graphics device driver to save all plots in a pdf format
        pdf(file = "CrossValidationExample.pdf",width = 8, height = 6)

        # Get the stage C prostate cancer data from the rpart package
        data(stagec,package = "rpart")

        # Prepare the data. Create a model matrix with interactions but no event time 
        stagec$pgtime <- NULL
        stagec$eet <- as.factor(stagec$eet)
        options(na.action = 'na.pass')
        stagec_mat <- cbind(pgstat = stagec$pgstat,
                         as.data.frame(model.matrix(pgstat ~ .*.,stagec))[-1])
        fnames <- colnames(stagec_mat)
        fnames <- str_replace_all(fnames,":","__")
        colnames(stagec_mat) <- fnames

        # Impute the missing data
        dataCancerImputed <- nearestNeighborImpute(stagec_mat)
        dataCancerImputed[,1:ncol(dataCancerImputed)] <- sapply(dataCancerImputed,as.numeric)

        # Cross validating a Random Forest classifier
        cvRF <- randomCV(dataCancerImputed,"pgstat",
                         randomForest::randomForest,
                         trainFraction = 0.8, 
                         repetitions = 10,
                         asFactor = TRUE);

        # Evaluate the prediction performance of the Random Forest classifier
        RFStats <- predictionStats_binary(cvRF$medianTest,
        plotname = "Random Forest",cex = 0.9);

        # Cross validating a BSWiMS with the same train/test set
        cvBSWiMS <- randomCV(fittingFunction = BSWiMS.model,
        	trainSampleSets = cvRF$trainSamplesSets);

        # Evaluate the prediction performance of the BSWiMS classifier
        BSWiMSStats <- predictionStats_binary(cvBSWiMS$medianTest,
        	plotname = "BSWiMS",cex = 0.9);

        # Cross validating a LDA classifier with a t-student filter
        cvLDA <- randomCV(dataCancerImputed,"pgstat",MASS::lda,
                          trainSampleSets = cvRF$trainSamplesSets,
                          featureSelectionFunction = univariate_tstudent,
                          featureSelection.control = list(limit = 0.5,thr = 0.975));

        # Evaluate the prediction performance of the LDA classifier
        LDAStats <- predictionStats_binary(cvLDA$medianTest,plotname = "LDA",cex = 0.9);

        # Cross validating a QDA classifier with LDA t-student features and RF train/test set
        cvQDA <- randomCV(fittingFunction = MASS::qda,
                          trainSampleSets = cvRF$trainSamplesSets,
                          featureSelectionFunction = cvLDA$selectedFeaturesSet);

        # Evaluate the prediction performance of the QDA classifier
        QDAStats <- predictionStats_binary(cvQDA$medianTest,plotname = "QDA",cex = 0.9);

        #Create a barplot with 95%CI that compare the balance error of the classifiers
        errorciTable <- rbind(RFStats$berror,
        	BSWiMSStats$berror,
        	LDAStats$berror,
        	QDAStats$berror)

        bpCI <- barPlotCiError(as.matrix(errorciTable),metricname = "Balanced Error",
                        	   thesets =  c("Classifier Method"),
                        	   themethod = c("RF","BSWiMS","LDA","QDA"),
                        	   main = "Balanced Error",
                        	   offsets = c(0.5,0.15),
                        	   scoreDirection = "<",
                        	   ho = 0.5,
                        	   args.legend = list(bg = "white",x = "topright"),
                        	   col = terrain.colors(4));



        dev.off()
	}
}



\keyword{Model_CV}
