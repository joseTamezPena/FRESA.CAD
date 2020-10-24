\name{FilterUnivariate}
\alias{FilterUnivariate}
\alias{univariate_Logit}
\alias{univariate_residual}
\alias{univariate_Wilcoxon}
\alias{univariate_tstudent}
\alias{univariate_correlation}
\alias{correlated_Remove}
\alias{univariate_cox}
\title{Univariate Filters}

\description{
    Returns the top set of features that are statistically associated with the outcome.
}
\usage{

univariate_Logit(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH", 
                 uniTest=c("zIDI","zNRI"),limit=0,...,n=0)
univariate_residual(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                    uniTest=c("Ftest","Binomial","Wilcox","tStudent"),
                    type=c("LM","LOGIT"),limit=0,...,n=0)
univariate_tstudent(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                    limit=0,...,n=0)
univariate_Wilcoxon(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                     limit=0,...,n=0)
univariate_KS(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                     limit=0,...,n=0)
univariate_correlation(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                       method = "kendall",limit=0,...,n=0)
univariate_cox(data=NULL, Outcome=NULL, pvalue=0.2, adjustMethod="BH",
                     limit=0,...)
univariate_BinEnsemble(data=NULL, Outcome=NULL,...)
correlated_Remove(data=NULL,fnames=NULL,thr=0.999)

}
\arguments{
    \item{data}{
    The data frame
    }
    \item{Outcome}{
    The outcome feature
    }
    \item{pvalue}{
    The threshold pvalue used after the p.adjust method
    }
    \item{adjustMethod}{
    The method used by the p.adjust method
    }
    \item{uniTest}{
    The unitTest to be performed by the linear fitting model
    }
    \item{type}{
    The type of linear model: LM or LOGIT
    }
    \item{method}{
    The correlation method: pearson,spearman or kendall.
    }
    \item{limit}{
    The samples-wise fraction of features to return. 
    }
    \item{fnames}{
    The list of features to test inside the correlated_Remove function
    }
    \item{thr}{
    The maximum correlation to allow between features
    }
    \item{...}{
    Parameters to be passed to the correlated_Remove function
    }
    \item{n}{
	the number of original features passed to p.adjust
    }
}

\value{
Named vector with the adjusted p-values or the list of no-correlated features for the correlated_Remove 
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

        q_values <- univariate_Logit(data=dataCancerImputed, 
                                    Outcome="pgstat",
                                    pvalue = 0.05)

        qValueMatrix <- q_values
        idiqValueMatrix <- q_values
        barplot(-log(q_values),las=2,cex.names=0.4,ylab="-log(Q)",
        main="Association with PGStat: IDI Test")

        q_values <- univariate_Logit(data=dataCancerImputed, 
                                    Outcome="pgstat", 
                                    uniTest="zNRI",pvalue = 0.05)
        qValueMatrix <- cbind(idiqValueMatrix,q_values[names(idiqValueMatrix)])

        q_values <- univariate_residual(data=dataCancerImputed, 
                                    Outcome="pgstat", 
                                    pvalue = 0.05,type="LOGIT")
        qValueMatrix <- cbind(qValueMatrix,q_values[names(idiqValueMatrix)])

        q_values <- univariate_tstudent(data=dataCancerImputed, 
                                       Outcome="pgstat", 
                                       pvalue = 0.05)
        qValueMatrix <- cbind(qValueMatrix,q_values[names(idiqValueMatrix)])

        q_values <- univariate_Wilcoxon(data=dataCancerImputed, 
                                       Outcome="pgstat", 
                                       pvalue = 0.05)
        qValueMatrix <- cbind(qValueMatrix,q_values[names(idiqValueMatrix)])

        q_values <- univariate_correlation(data=dataCancerImputed, 
                                       Outcome="pgstat", 
                                       pvalue = 0.05)
        qValueMatrix <- cbind(qValueMatrix,q_values[names(idiqValueMatrix)])

        q_values <- univariate_correlation(data=dataCancerImputed, 
                                          Outcome="pgstat", 
                                          pvalue = 0.05,
                                          method = "pearson")

        #The qValueMatrix has the qValues of all filter methods.  
        qValueMatrix <- cbind(qValueMatrix,q_values[names(idiqValueMatrix)])
        colnames(qValueMatrix) <- c("IDI","NRI","F","t","W","K","P")
        #Do the log transform to display the heatmap
        qValueMatrix <- -log10(qValueMatrix)
        #the Heatmap of the q-values
        gplots::heatmap.2(qValueMatrix,Rowv = FALSE,dendrogram = "col",
        main = "Method q.values",cexRow = 0.4)

    }
}
\keyword{Feature_Filtering}
