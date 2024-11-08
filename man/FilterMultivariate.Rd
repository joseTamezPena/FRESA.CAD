\name{multivariate_BinEnsemble}
\alias{multivariate_BinEnsemble}

\title{Multivariate Filters}

\description{
    Returns the top set of features that are associated with the outcome based on Multivariate logistic models: LASSO and BSWiMS
}
\usage{
multivariate_BinEnsemble(data,Outcome,pvalue=0.2,limit=-1,adjustMethod="BH",...)
}
\arguments{
    \item{data}{
    The data frame
    }
    \item{Outcome}{
    The outcome feature
    }
    \item{adjustMethod}{
    The method used by the p.adjust method
    }
    \item{pvalue}{
    The maximum adjusted p-value
    }
    \item{limit}{
    The samples-wise fraction of features to return. 
    }
    \item{...}{
    Parameters to be passed to the correlated_Remove function
    }
}

\value{
    Named vector with the adjusted p-values of the associted features
}
\author{Jose G. Tamez-Pena}
\examples{
    \dontrun{

        library("FRESA.CAD")

        ### Univariate Filter Examples ####

        # Get the stage C prostate cancer data from the rpart package
        data(stagec,package = "rpart")

        # Prepare the data. Create a model matrix without the event time and interactions
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

        # Get the top Features associated to pgstat

        q_values <- multivariate_BinEnsemble(data=dataCancerImputed, 
                                    Outcome="pgstat")



    }
}
\keyword{Feature_Selection}

