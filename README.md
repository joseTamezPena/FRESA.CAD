# FRESA.CAD

Feature Selection Algorithms for Computer Aided Diagnosis.

Set of functions for: Conditioning, Feature Selection, Machine Learning, Cross-Validation, and Visual Evaluation

## Table of Contents

-   [Overview](#Overview)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Contributing](#contributing)
-   [License](#license)
-   [Contact](#contact)

## Overview

The design of diagnostic or prognostic multivariate models via the selection of significantly discriminant features is complex.

FRESA.CAD provides a series of functions for: Data conditioning, Feature Selection, Machine Learning, Benchmarking, Visualization and Reporting.

| Category                       | Function(s)                | Purpose                                   |
|---------------------|---------------------|------------------------------|
| **Conditioning/Preprocessing** | nearestNeighborImpute()    | Impute missing values                     |
| **Conditioning/Preprocessing** | FRESA.Scale()              | Data Scale/Normalization                  |
| **Conditioning/Preprocessing** | featureAdjustment()        | Adjust variables removing collinearity    |
| **Conditioning/Preprocessing** | IDeA()/ILAA()              | Multicollinearity Mitigation              |
|                                |                            |                                           |
| **Feature Selection**          | uniRankVar()               | Univariate Analysis                       |
| **Feature Selection**          | BSWiMS.model()             | Linear Model Subset Selection             |
| **Feature Selection**          | univariate_BinEnsemble()   | Ensemble Select Top Features              |
| **Feature Selection**          | univariate...              | Filter Select Top Features ...            |
|                                |                            |                                           |
| **Machine Learning**           | BSWiMS.model()             | Bootstrap Modeling                        |
| **Machine Learning**           | filteredFit()              | Pipeline ML: Scale/Filter/Transform/Learn |
| **Machine Learning**           | HLCM()/HLCM_EM()           | Latent-Class Based Modeling               |
| **Machine Learning**           | GMVECluster()              | Unsupervised Clustering via GMVE          |
|                                |                            |                                           |
| **Benchmarking / Evaluation**  | RandomCV()                 | Random Holdout Validation                 |
| **Benchmarking / Evaluation**  | BinaryBenchmark()          | Binary Model Evaluation                   |
| **Benchmarking / Evaluation**  | OrdinalBenchmark()         | Ordinal Model Evaluation                  |
| **Benchmarking / Evaluation**  | CoxBenchmark()             | Cox-based Model Evaluation                |
|                                |                            |                                           |
| **Visualization / Reporting**  | RRPlot()                   | Survival Model Evaluation                 |
| **Visualization / Reporting**  | predictionStats_binary()   | Report Cross Validation Results Binary    |
| **Visualization / Reporting**  | predictionStats_Ordinal()  | Report Cross Validation Results Ordinal   |
| **Visualization / Reporting**  | predictionStats_survival() | Report Cross Validation Results Survival  |

Besides the above listed functions the library provides predictors and wrappers of common machine learning methods, and many other auxiliary functions.

## Installation

You can install the official release of the package from CRAN using:

``` r
install.packages("FRESA.CAD")
```

To install the development version from GitHub, use:

``` r
# Install 'devtools' package if you haven't already
install.packages("devtools")

# Install the package from GitHub
devtools::install_github("https://github.com/joseTamezPena/FRESA.CAD")
```

## Usage

``` r
#Load the package
library(FRESA.CAD)

#For comprehensive evaluaiton of confusion tables
library("epiR")

# Example usage

data(stagec,package = "rpart")
options(na.action = 'na.pass')
dataCancer <- cbind(pgstat = stagec$pgstat,
                        pgtime = stagec$pgtime,
                        as.data.frame(
                          model.matrix(Surv(pgtime,pgstat) ~ .,stagec))[-1])

#Impute missing values
dataCancerImputed <- nearestNeighborImpute(dataCancer)
data(cancerVarNames)

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
#The list of all models of the bootstrap forward selection 
print(BSWiMSModel$forward.selection.list)

#With FRESA.CAD we can do a leave-one-out using the list of models
pm <- ensemblePredict(BSWiMSModel$forward.selection.list,
                          dataCancer,predictType = "linear",type="LOGIT",Outcome="pgstat")

#Ploting the ROC with 95
pm <- plotModels.ROC(cbind(dataCancer$pgstat,
                               pm$ensemblePredict),
                     main=("LOO Forward Selection Median Predict"))

#The plotModels.ROC provides the diagnosis confusion matrix.
summary(epi.tests(pm$predictionTable))
    
```

More examples of FRESA.CAD usage can be found at: <https://rpubs.com/J_Tamez>

## Contributing

Contributions are welcome! If you'd like to contribute to this project, please follow these guidelines:

\- Fork the repository.

\- Create a new branch: `git checkout -b feature/new-feature`.

\- Make your changes and commit them: `git commit -m 'Add new feature'`.

\- Push to the branch: `git push origin feature/new-feature`.

\- Submit a pull request.

## License

This project is licensed under the LGPL\>=2.0.

## Contact

For any questions or feedback, feel free to contact us at:

Email: jose.tamezpena\@tec.mx

Twitter: [\@tamezpena](https://twitter.com/jtamezpena)
