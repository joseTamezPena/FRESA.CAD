#' @method plot FRESA_benchmark
plot.FRESA_benchmark <-
function(x,...) 
{
	prefix <- "";
	parameters <- list(...);
	if (!is.null(parameters$prefix))
	{
		prefix <- parameters$prefix;
	}

	op <- par(no.readonly=TRUE);
	result = NULL;

	if (class(x)[2] == "Binary")
	{

		args.legend = list(bg = "white",x="bottomright")
		mar = par("mar")
		if (mar[4] >= 8)
		{
			args.legend = list(bg = "white",x="bottomright",inset=c(-0.25,0),cex=0.75)
		}

		x$errorciTable[is.na(x$errorciTable)] <- 0;
		bpBER <- barPlotCiError(as.matrix(x$errorciTable),metricname = "Balanced Error",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Balanced Error"),offsets = c(0.5,1),scoreDirection = "<",ho=0.5,args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testBalancedError <- bpBER$ciTable$mean;
		testBalancedErrormin <- min(bpBER$ciTable$low95);
		testBalancedErrormax <- max(bpBER$ciTable$top95);

		x$accciTable[is.na(x$accciTable)] <- 0;
		bpACC <- barPlotCiError(as.matrix(x$accciTable),metricname = "Accuracy",thesets = x$thesets,themethod = x$theMethod,main =  paste(prefix,"Accuracy"),offsets = c(0.5,1),args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testACC <- bpACC$ciTable$mean;
		testACCmin <- min(bpACC$ciTable$low95);
		testACCmax <- max(bpACC$ciTable$top95);

		x$aucTable[is.na(x$aucTable)] <- 0;
		bpAUC <- barPlotCiError(as.matrix(x$aucTable),metricname = "AUC",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"ROC AUC"),offsets = c(0.5,1),ho=0.5,args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testAUC <- bpAUC$ciTable$mean;
		testAUCmin <- min(bpAUC$ciTable$low95);
		testAUCmax <- max(bpAUC$ciTable$top95);

		x$cidxTable[is.na(x$cidxTable)] <- 0;
		bpCIDX <- barPlotCiError(as.matrix(x$cidxTable),metricname = "CIndex",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Concordance"),offsets = c(0.5,1),ho=0.5,args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testCIDX <- bpCIDX$ciTable$mean;
		testCIDXmin <- min(bpCIDX$ciTable$low95);
		testCIDXmax <- max(bpCIDX$ciTable$top95);

		x$senTable[is.na(x$senTable)] <- 0;
		bpSEN <- barPlotCiError(as.matrix(x$senTable),metricname = "Sensitivity",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Sensitivity"),offsets = c(0.5,1),args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testSEN <- bpSEN$ciTable$mean;
		testSENmin <- min(bpSEN$ciTable$low95);
		testSENmax <- max(bpSEN$ciTable$top95);

		x$speTable[is.na(x$speTable)] <- 0;
		bpSPE <- barPlotCiError(as.matrix(x$speTable),metricname = "Specificity",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Specificity"),offsets = c(0.5,1),args.legend = args.legend,col = terrain.colors(length(x$theMethod)),...);
		testSPE <- bpSPE$ciTable$mean;
		testSPEmin <- min(bpSPE$ciTable$low95);
		testSPEmax <- max(bpSPE$ciTable$top95);
		brnames <- names(testBalancedError);
		

		
		metrics <- rbind(BER = testBalancedError,
							ACC = testACC[brnames],
							AUC = testAUC[brnames],
							SEN = testSEN[brnames],
							SPE = testSPE[brnames],
							CIDX = testCIDX[brnames]
						);
		barPlotsCI <- list(BER=bpBER,
							ACC=bpACC,
							AUC=bpAUC,
							SEN=bpSEN,
							SPE=bpSPE,
							CIDX=bpCIDX
							);

		metrics_filter <- NULL;
		barPlotsCI_filter <- NULL;
		minMaxMetrics <- NULL;
							
		if (!is.null(x$errorciTable_filter))
		{
			par(op)
			barplot(x$jaccard[order(-x$jaccard)],las = 2,cex.axis = 1,cex.names = 0.7,main = paste(prefix,"Jaccard Index"),ylab = "Jaccard")
			unsize <- x$featsize
			unsize[unsize == 0] <- 1;
			barplot(unsize[order(unsize)],las = 2,cex.axis = 1,cex.names = 0.7,log = "y",main = paste(prefix,"Number of Features"),ylab = "# of Features")
			par(op)
			x$errorciTable_filter[is.na(x$errorciTable_filter)] <- 0;
			bpBER_filter <- barPlotCiError(as.matrix(x$errorciTable_filter),metricname = "Balanced Error",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"Balanced Error"),scoreDirection = "<",ho=0.5,args.legend = list(bg = "white",x = "topleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilBalancedError <- apply(bpBER_filter$ciTable$mean,2,median);
			testFilBalancedErrormax <- max(apply(bpBER_filter$ciTable$top95,2,max));
			testFilBalancedErrormin <- min(apply(bpBER_filter$ciTable$low95,2,min))

			x$accciTable_filter[is.na(x$accciTable_filter)] <- 0;
			bpACC_filter <- barPlotCiError(as.matrix(x$accciTable_filter),metricname = "Accuracy",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"Accuracy"),offsets = c(0.5,1),args.legend = list(bg = "white",x = "bottomleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilACC <- apply(bpACC_filter$ciTable$mean,2,median);
			testFilACCmin <- min(apply(bpACC_filter$ciTable$low95,2,min));
			testFilACCmax <- max(apply(bpACC_filter$ciTable$top95,2,max));

			x$aucTable_filter[is.na(x$aucTable_filter)] <- 0;
			bpAUC_filter <- barPlotCiError(as.matrix(x$aucTable_filter),metricname = "AUC",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"ROC AUC"),offsets = c(0.5,1),ho=0.5,args.legend = list(bg = "white",x = "bottomleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilAUC <- apply(bpAUC_filter$ciTable$mean,2,median);
			testFilAUCmin <- min(apply(bpAUC_filter$ciTable$low95,2,min));
			testFilAUCmax <- max(apply(bpAUC_filter$ciTable$top95,2,max));

			x$cindexTable_filter[is.na(x$cindexTable_filter)] <- 0;
			bpCIDX_filter <- barPlotCiError(as.matrix(x$cindexTable_filter),metricname = "CIDX",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"Concordance"),offsets = c(0.5,1),ho=0.5,args.legend = list(bg = "white",x = "bottomleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilCIDX <- apply(bpCIDX_filter$ciTable$mean,2,median);
			testFilCIDXmin <- min(apply(bpCIDX_filter$ciTable$low95,2,min));
			testFilCIDXmax <- max(apply(bpCIDX_filter$ciTable$top95,2,max));
			
			x$senciTable_filter[is.na(x$senciTable_filter)] <- 0;
			bpSEN_filter <- barPlotCiError(as.matrix(x$senciTable_filter),metricname = "Sensitivity",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"Sensitivity"),offsets = c(0.5,1),args.legend = list(bg = "white",x = "bottomleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilSEN <- apply(bpSEN_filter$ciTable$mean,2,median);
			testFilSENmin <- min(apply(bpSEN_filter$ciTable$low95,2,min));
			testFilSENmax <- max(apply(bpSEN_filter$ciTable$top95,2,max));

			x$speciTable_filter[is.na(x$speciTable_filter)] <- 0;
			bpSPE_filter <- barPlotCiError(as.matrix(x$speciTable_filter),metricname = "Specificity",thesets = x$theFiltersets,themethod = x$theClassMethod,main = paste(prefix,"Specificity"),offsets = c(0.5,1),args.legend = list(bg = "white",x = "bottomleft"),col = terrain.colors(length(x$theClassMethod)),...)
			testFilSPE <- apply(bpSPE_filter$ciTable$mean,2,median);
			testFilSPEmin <- min(apply(bpSPE_filter$ciTable$low95,2,min));
			testFilSPEmax <- max(apply(bpSPE_filter$ciTable$top95,2,max));
			brnames <- names(testFilBalancedError);
			metrics_filter <- rbind(BER = testFilBalancedError,
								ACC = testFilACC[brnames],
								AUC = testFilAUC[brnames],
								SEN = testFilSEN[brnames],
								SPE = testFilSPE[brnames],
								CIDX = testFilCIDX[brnames]
								);

			barPlotsCI_filter <- list(BER=bpBER_filter,
									ACC=bpACC_filter,
									SEN=bpSEN_filter,
									AUC=bpAUC_filter,
									SPE=bpSPE_filter,
									CIDX=bpCIDX_filter
								);
			minMaxMetrics <- list(BER = c(min(testBalancedErrormin,testFilBalancedErrormin),max(testBalancedErrormax,testFilBalancedErrormax)),
						ACC = c(min(testACCmin,testFilACCmin),max(testACCmax,testFilACCmax)),
						AUC = c(min(testAUCmin,testFilAUCmin),max(testAUCmax,testFilAUCmax)),
						SEN = c(min(testSENmin,testFilSENmin),max(testSENmax,testFilSENmax)),
						SPE = c(min(testSPEmin,testFilSPEmin),max(testSPEmax,testFilSPEmax)),
						CIDX = c(min(testFilCIDXmin,testFilCIDXmin),max(testFilCIDXmax,testFilCIDXmax))
						);
		}
		mcnemar <- matrix(0,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		pmcnemar <- matrix(1,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		if ((ncol(x$testPredictions)-1) > 2)
		{
			for (i in 2:(ncol(x$testPredictions)-1))
			{
				for (j in (i+1):ncol(x$testPredictions))
				{
					tb <- table(x$testPredictions[,i] > 0.5,x$testPredictions[,j] > 0.5)
					pmcnemar[i-1,j-1] <- epiR::epi.kappa(tb)$mcnemar$p.value;
					pmcnemar[j-1,i-1] <- pmcnemar[i-1,j-1];
					mcnemar[i-1,j-1] <- -log10(max(pmcnemar[i-1,j-1],0.0001));
					mcnemar[j-1,i-1] <- mcnemar[i-1,j-1];
				}
			}
			mcnemar[is.nan(mcnemar)] <- 0;
			colnames(mcnemar) <- colnames(x$testPredictions)[-1]
			rownames(mcnemar) <- colnames(x$testPredictions)[-1]
			colnames(pmcnemar) <- colnames(x$testPredictions)[-1]
			rownames(pmcnemar) <- colnames(x$testPredictions)[-1]
			par(op);
			par(mfrow = c(1,1),mar = c(2,2,2,2));
			gplots::heatmap.2(mcnemar,trace = "none",mar = c(5,10),col=rev(heat.colors(8)),main = "p(Method A = Method B)",cexRow = 0.65,cexCol = 0.75,srtCol = 25,key.xlab="-log(p)",xlab="Method B", ylab="Method A")
			par(op);
		}

		AUCtable <- matrix(0,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		pAUCtable <- matrix(1,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		if ((ncol(x$testPredictions)-1) > 2)
		{
			for (i in 2:ncol(x$testPredictions))
			{
				for (j in 2:ncol(x$testPredictions))
				{
					roct <- roc.test(roc(x$testPredictions$Outcome,x$testPredictions[,i],quiet = TRUE),
									roc(x$testPredictions$Outcome,x$testPredictions[,j],quiet = TRUE),
									alternative="less")
					pAUCtable[i-1,j-1] <-  roct$p.value;
					AUCtable[i-1,j-1] <- -log10(max(pAUCtable[i-1,j-1],0.0001));
				}
			}
			AUCtable[is.nan(AUCtable)] <- 0.0;
			colnames(AUCtable) <- colnames(x$testPredictions)[-1]
			rownames(AUCtable) <- colnames(x$testPredictions)[-1]
			colnames(pAUCtable) <- colnames(x$testPredictions)[-1]
			rownames(pAUCtable) <- colnames(x$testPredictions)[-1]
			par(op);
			par(mfrow = c(1,1),mar = c(2,2,2,2));
			topf <- apply(AUCtable,1,mean)
			gplots::heatmap.2(AUCtable[order(topf),],trace = "none",mar = c(5,10),col=rev(heat.colors(8)),Rowv=FALSE,dendrogram = "column",cexRow = 0.65,cexCol = 0.75,srtCol = 25,key.xlab="-log(p)",main = "p(ROC_AUC A > ROC_AUC B)",xlab="Method B",ylab="Method A")
			par(op);
		}
		
		result <- list(metrics = metrics, barPlotsCI = barPlotsCI,metrics_filter=metrics_filter,barPlotsCI_filter=barPlotsCI_filter, minMaxMetrics = minMaxMetrics,mcnemar=pmcnemar,AUCtable=pAUCtable);
	}
	if (class(x)[2] == "Ordinal")
	{
		x$KappaTable[is.na(x$KappaTable)] <- 0;
		x$BMAETable[is.na(x$x$BMAETable)] <- 0;
		x$KendallTable[is.na(x$KendallTable)] <- 0;
		x$BiasTable[is.na(x$BiasTable)] <- 0;
		x$ACCTable[is.na(x$ACCTable)] <- 0;
		x$SENTable[is.na(x$SENTable)] <- 0;
		x$AUCTable[is.na(x$AUCTable)] <- 0;
		x$KappaTable_filter[is.na(x$KappaTable_filter)] <- 0;
		x$KendallTable_filter[is.na(x$KendallTable_filter)] <- 0;
		x$BMAETable_filter[is.na(x$BMAETable_filter)] <- 0;
		x$ACCTable_filter[is.na(x$ACCTable_filter)] <- 0;
		x$SENTable_filter[is.na(x$SENTable_filter)] <- 0;
		x$AUCTable_filter[is.na(x$AUCTable_filter)] <- 0;

		bpKAPPA <- barPlotCiError(as.matrix(x$KappaTable),metricname = "KAPPA",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"KAPPA"),offsets = c(0.5,0.05),ho=0.0,args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testKAPPA <- bpKAPPA$ciTable$mean;
		testKAPPAmin <- min(bpKAPPA$ciTable$low95)
		testKAPPAmax <- max(bpKAPPA$ciTable$top95)

		bpKendall <- barPlotCiError(as.matrix(x$KendallTable),metricname = "Kendall Correlation",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Kendall Correlation"),ho=0.0,args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testKendall <- bpKendall$ciTable$mean;
		testKendallmin <- min(bpKendall$ciTable$low95)
		testKendallmax <- max(bpKendall$ciTable$top95)

		bpBMAE <- barPlotCiError(as.matrix(x$BMAETable),metricname = "BMAE",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"BMAE"),scoreDirection = "<",args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testBMAE <- bpBMAE$ciTable$mean;
		testBMAEmin <- min(bpBMAE$ciTable$low95)
		testBMAEmax <- max(bpBMAE$ciTable$top95)

		bpBias <- barPlotCiError(as.matrix(x$BiasTable),metricname = "Bias",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Bias"),offsets = c(0.5,0.5),scoreDirection = "==",args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testBias <- bpBias$ciTable$mean;
		testBiasmin <- min(bpBias$ciTable$low95)
		testBiasmax <- max(bpBias$ciTable$top95)

		bpACC <- barPlotCiError(as.matrix(x$ACCTable),metricname = "Accuracy",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Accuracy"),offsets = c(0.5,0.5),args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testACC <- bpACC$ciTable$mean;
		testACCmin <- min(bpACC$ciTable$low95)
		testACCmax <- max(bpACC$ciTable$top95)

		bpSEN <- barPlotCiError(as.matrix(x$SENTable),metricname = "Sensitivity",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Sensitivity"),offsets = c(0.5,0.5),args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testSEN <- bpSEN$ciTable$mean;
		testSENmin <- min(bpSEN$ciTable$low95)
		testSENmax <- max(bpSEN$ciTable$top95)

		bpAUC <- barPlotCiError(as.matrix(x$AUCTable),metricname = "ROC AUC",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"ROC AUC"),offsets = c(0.5,0.5),ho=0.5,args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testAUC <- bpAUC$ciTable$mean;
		testAUCmin <- min(bpAUC$ciTable$low95)
		testAUCmax <- max(bpAUC$ciTable$top95)
		
		brnames <- names(testKAPPA);
		metrics <- rbind(KAPPA = testKAPPA,
							BMAE = testBMAE[brnames],
							Kendall = testKendall[brnames],
							Bias = testBias[brnames],
							ACC = testACC[brnames],
							SEN = testSEN[brnames], 
							AUC = testAUC[brnames]
						);
		barPlotsCI <- list(KAPPA=bpKAPPA,
							BMAE=bpBMAE,
							Kendall=bpKendall,
							Bias=bpBias,
							ACC=bpACC,
							SEN=bpSEN,
							AUC=bpAUC
						);
		
		bpKAPPA_filter <- barPlotCiError(as.matrix(x$KappaTable_filter),metricname = "Kappa",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"Kappa"),offsets = c(0.0,0.1),ho=0.0,args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilKappa <- apply(bpKAPPA_filter$ciTable$mean,2,median);
		testFilKappamin <- min(apply(bpKAPPA_filter$ciTable$low95,2,min));
		testFilKappamax <- max(apply(bpKAPPA_filter$ciTable$top95,2,max));

		bpKendall_filter <- barPlotCiError(as.matrix(x$KendallTable_filter),metricname = "Kendall Correlation",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"Kendall Correlation"),offsets = c(0.0,0.1),ho=0.0,args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilKendall <- apply(bpKendall_filter$ciTable$mean,2,median);
		testFilKendallmin <- min(apply(bpKendall_filter$ciTable$low95,2,min));
		testFilKendallmax <- max(apply(bpKendall_filter$ciTable$top95,2,max));

		bpBMAE_filter <- barPlotCiError(as.matrix(x$BMAETable_filter),metricname = "BMAE",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"BMAE"),scoreDirection = "<",offsets = c(0.0,0.1),args.legend = list(bg = "white",x = "bottomright"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilBMAE <- apply(bpBMAE_filter$ciTable$mean,2,median);
		testFilBMAEmin <- min(apply(bpBMAE_filter$ciTable$low95,2,min));
		testFilBMAEmax <- max(apply(bpBMAE_filter$ciTable$top95,2,max));

		bpACC_filter <- barPlotCiError(as.matrix(x$ACCTable_filter),metricname = "Accuracy ",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"Accuracy "),offsets = c(0.0,0.1),args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilACC <- apply(bpACC_filter$ciTable$mean,2,median);
		testFilACCmin <- min(apply(bpACC_filter$ciTable$low95,2,min));
		testFilACCmax <- max(apply(bpACC_filter$ciTable$top95,2,max));

		bpSEN_filter <- barPlotCiError(as.matrix(x$SENTable_filter),metricname = "Sensitivity",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"Sensitivity"),offsets = c(0.0,0.1),args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilSEN <- apply(bpSEN_filter$ciTable$mean,2,median);
		testFilSENmin <- min(apply(bpSEN_filter$ciTable$low95,2,min));
		testFilSENmax <- max(apply(bpSEN_filter$ciTable$top95,2,max));

		bpAUC_filter <- barPlotCiError(as.matrix(x$AUCTable_filter),metricname = "ROC AUC ",thesets = x$theFiltersets,themethod = x$theOrdinalMethod,main = paste(prefix,"ROC AUC "),offsets = c(0.0,0.1),ho=0.5,args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theOrdinalMethod)),...)
		testFilAUC <- apply(bpAUC_filter$ciTable$mean,2,median);
		testFilAUCmin <- min(apply(bpAUC_filter$ciTable$low95,2,min));
		testFilAUCmax <- max(apply(bpAUC_filter$ciTable$top95,2,max));
		brnames <- names(testFilKappa);

		metrics_filter <- rbind(KAPPA = testFilKappa,
								BMAE = testFilBMAE[brnames],
								Kendall = testFilKendall[brnames],
								ACC = testFilACC[brnames],
								SEN = testFilSEN[brnames], 
								AUC = testFilAUC[brnames]
								);
		barPlotsCI_filter <- list(KAPPA=bpKAPPA_filter,
									BMAE=bpBMAE_filter,
									Kendall=bpKendall_filter,
									ACC=bpACC_filter,
									SEN=bpSEN_filter,
									AUC=bpAUC_filter
									);
		minMaxMetrics <- list(KAPPA = c(min(testKAPPAmin,testFilKappamin),max(testKAPPAmax,testFilKappamax)),
						BMAE = c(min(testBMAEmin,testFilBMAEmin),max(testBMAEmax,testFilBMAEmax)),
						Kendall = c(min(testKendallmin,testFilKendallmin),max(testKendallmax,testFilKendallmax)),
						Bias = c(testBiasmin,testBiasmax),
						ACC = c(min(testACCmin,testFilACCmin),max(testACCmax,testFilACCmax)),
						SEN =c(min(testSENmin,testFilSENmin),max(testSENmax,testFilSENmax)),
						AUC =c(min(testAUCmin,testFilAUCmin),max(testAUCmax,testFilAUCmax))
						);
		fttable <- matrix(0,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		pfttable <- matrix(1,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		if ((ncol(x$testPredictions)-1) > 2)
		{
			for (i in 2:ncol(x$testPredictions))
			{
				for (j in 2:ncol(x$testPredictions))
				{
					rss1 <- max(sum((x$testPredictions[,j]-x$testPredictions$Outcome)^2),1e-10);
					rss2 <- max(sum((x$testPredictions[,i]-x$testPredictions$Outcome)^2),1e-10);
					rss2 <- rss2/rss1;
					pfttable[i-1,j-1] <-  pf(nrow(cp$testPredictions)*(rss2-1.0),1.0,nrow(cp$testPredictions),lower.tail = FALSE);
					fttable[i-1,j-1] <- -log10(max(pfttable[i-1,j-1],0.0001));
				}
			}
			fttable[is.nan(fttable)] <- 0.0;
			colnames(fttable) <- colnames(x$testPredictions)[-1]
			rownames(fttable) <- colnames(x$testPredictions)[-1]
			colnames(pfttable) <- colnames(x$testPredictions)[-1]
			rownames(pfttable) <- colnames(x$testPredictions)[-1]
			par(op);
			par(mfrow = c(1,1),mar = c(2,2,2,2));
			topf <- apply(fttable,1,mean)
			gplots::heatmap.2(fttable[order(topf),],trace = "none",mar = c(5,10),col=rev(heat.colors(8)),Rowv=FALSE,dendrogram = "column",cexRow = 0.65,cexCol = 0.75,srtCol = 25,key.xlab="-log(p)",main = "p(Method A > Method B)",xlab="Method B",ylab="Method A")
			par(op);
		}


		result <- list(metrics = metrics, barPlotsCI = barPlotsCI,metrics_filter=metrics_filter,barPlotsCI_filter=barPlotsCI_filter, minMaxMetrics = minMaxMetrics,fttable=pfttable);
	}
	if (class(x)[2] == "Regression")
	{
		x$MAETable[is.na(x$MAETable)] <- 0;
		bpMAE <- barPlotCiError(as.matrix(x$MAETable),metricname = "MAE",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"MAE"),scoreDirection = "<",args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testMAE <- bpMAE$ciTable$mean;
		testMAEmax <- max(bpMAE$ciTable$top95)
		testMAEmin <- min(bpMAE$ciTable$low95)

		x$CorSpearman[is.na(x$CorSpearman)] <- 0;
		bpSpearman <- barPlotCiError(as.matrix(x$CorSpearman),metricname = "Spearman Correlation",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Spearman Correlation"),offsets = c(0.5,0.05),ho=0.0,args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testSpearman <- bpSpearman$ciTable$mean;
		testSpearmanmin <- min(bpSpearman$ciTable$low95)
		testSpearmanmax <- max(bpSpearman$ciTable$top95)

		x$RMSETable[is.na(x$RMSETable)] <- 0;
		bpRMSE <- barPlotCiError(as.matrix(x$RMSETable),metricname = "RMSE",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"RMSE"),offsets = c(0.5,0.5),scoreDirection = "<",args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testRMSE <- bpRMSE$ciTable$mean;
		testRMSEmax <- max(bpRMSE$ciTable$top95)
		testRMSEmin <- min(bpRMSE$ciTable$low95)

		x$CorTable[is.na(x$CorTable)] <- 0;
		bpPearson <- barPlotCiError(as.matrix(x$CorTable),metricname = "Pearson Correlation",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Pearson Correlation"),ho=0.0,args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testPearson <- bpPearson$ciTable$mean;
		testPearsonmin <- min(bpPearson$ciTable$low95)
		testPearsonmax <- max(bpPearson$ciTable$top95)
		  
		x$BiasTable[is.na(x$BiasTable)] <- 0;
		bpBias <- barPlotCiError(as.matrix(x$BiasTable),metricname = "Bias",thesets = x$thesets,themethod = x$theMethod,main = paste(prefix,"Bias"),offsets = c(0.5,0.5),scoreDirection = "==",args.legend = list(bg = "white",x = "bottomright"),col = terrain.colors(length(x$theMethod)),...)
		testBias <- bpBias$ciTable$mean;
		testBiasmax <- max(bpBias$ciTable$top95)
		testBiasmin <- min(bpBias$ciTable$low95)

		brnames <- names(testSpearman);
		metrics <- rbind(Spearman = testSpearman,
							MAE = testMAE[brnames],
							Pearson = testPearson[brnames],
							RMSE = testRMSE[brnames],
							Bias = testBias[brnames]
							);
		barPlotsCI <- list(Spearman=bpSpearman,
							MAE=bpMAE,
							Pearson=bpPearson,
							RMSE=bpRMSE,
							Bias=bpBias
							);
		
		x$MAETable_filter[is.na(x$MAETable_filter)] <- 0;
		bpMAE_filter <- barPlotCiError(as.matrix(x$MAETable_filter),metricname = "MAE",thesets = x$theFiltersets,themethod = x$theRegressMethod,main = paste(prefix,"MAE"),scoreDirection = "<",offsets = c(0.0,1.0),args.legend = list(bg = "white",x = "bottomright"),angle = 45,col = terrain.colors(length(x$theRegressMethod)),...)
		testFilMAE <- apply(bpMAE_filter$ciTable$mean,2,median);
		testFilMAEmax <- max(apply(bpMAE_filter$ciTable$top95,2,max));
		testFilMAEmin <- min(apply(bpMAE_filter$ciTable$low95,2,min));

		x$CorSpearman_filter[is.na(x$CorSpearman_filter)] <- 0;
		bpSpearman_filter <- barPlotCiError(as.matrix(x$CorSpearman_filter),metricname = "Correlation",thesets = x$theFiltersets,themethod = x$theRegressMethod,main = paste(prefix,"Spearman Correlation"),offsets = c(0.0,0.1),ho=0.0,args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theRegressMethod)),...)
		testFilSpearman <- apply(bpSpearman_filter$ciTable$mean,2,median);
		testFilSpearmanmin <- min(apply(bpSpearman_filter$ciTable$low95,2,min));
		testFilSpearmanmax <- max(apply(bpSpearman_filter$ciTable$top95,2,max));

		x$RMSETable_filter[is.na(x$RMSETable_filter)] <- 0;
		bpRMSE_filter <- barPlotCiError(as.matrix(x$RMSETable_filter),metricname = "RMSE",thesets = x$theFiltersets,themethod = x$theRegressMethod,main = paste(prefix,"RMSE"),scoreDirection = "<",offsets = c(0.0,1.0),args.legend = list(bg = "white",x = "bottomright"),angle = 45,col = terrain.colors(length(x$theRegressMethod)),...)
		testFilRMSE <- apply(bpRMSE_filter$ciTable$mean,2,median);
		testFilRMSEmax <- max(apply(bpRMSE_filter$ciTable$top95,2,max));
		testFilRMSEmin <- min(apply(bpRMSE_filter$ciTable$low95,2,min));

		x$CorTable_filter[is.na(x$CorTable_filter)] <- 0;
		bpPearson_filter <- barPlotCiError(as.matrix(x$CorTable_filter),metricname = "Pearson Correlation",thesets = x$theFiltersets,themethod = x$theRegressMethod,main = paste(prefix,"Pearson Correlation"),offsets = c(0.0,0.1),ho=0.0,args.legend = list(bg = "white",x = "bottomleft"),angle = 45,col = terrain.colors(length(x$theRegressMethod)),...)
		testFilPearson <- apply(bpPearson_filter$ciTable$mean,2,median);
		testFilPearsonmin <- min(apply(bpPearson_filter$ciTable$low95,2,min));
		testFilPearsonmax <- max(apply(bpPearson_filter$ciTable$top95,2,max));
		  
		x$BiasTable_filter[is.na(x$BiasTable_filter)] <- 0;
		bpBias_filter <- barPlotCiError(as.matrix(x$BiasTable_filter),metricname = "Bias",thesets = x$theFiltersets,themethod = x$theRegressMethod,main = paste(prefix,"Bias"),offsets = c(0.5,1),args.legend = list(bg = "white",x = "bottomleft",cex = 0.75),angle = 45,scoreDirection = "==",col = terrain.colors(length(x$theRegressMethod)),...)
		testFilBias <- apply(bpBias_filter$ciTable$mean,2,median);
		testFilBiasmin <- min(apply(bpBias_filter$ciTable$low95,2,min));
		testFilBiasmax <- max(apply(bpBias_filter$ciTable$top95,2,max));

		brnames <- names(testFilSpearman);
		metrics_filter <- rbind(Spearman = testFilSpearman,
								MAE = testFilMAE[brnames],
								Pearson = testFilPearson[brnames],
								RMSE = testFilRMSE[brnames],
								Bias = testFilBias[brnames]
								);
		barPlotsCI_filter <- list(Spearman=bpSpearman_filter,
									MAE=bpMAE_filter,
									Pearson=bpPearson_filter,
									RMSE=bpRMSE_filter,
									Bias=bpBias_filter
									);
		minMaxMetrics <- list(MAE = c(min(c(testMAEmin,testFilMAEmin)),max(c(testMAEmax,testFilMAEmax))),
						Spearman = c(min(c(testSpearmanmin,testFilSpearmanmin)),max(testSpearmanmax,testFilSpearmanmax)),
						RMSE = c(min(c(testRMSEmin,testFilRMSEmin)),max(c(testRMSEmax,testFilRMSEmax))),
						Pearson = c(min(c(testPearsonmin,testFilPearsonmin)),max(c(testPearsonmax,testFilPearsonmax))),
						Bias = c(min(c(testBiasmin,testFilBiasmin)),max(c(testBiasmax,testFilBiasmax)))
						);
						

		fttable <- matrix(0,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		pfttable <- matrix(1,ncol = ncol(x$testPredictions)-1,nrow=ncol(x$testPredictions)-1)
		if ((ncol(x$testPredictions)-1) > 2)
		{
			for (i in 2:ncol(x$testPredictions))
			{
				for (j in 2:ncol(x$testPredictions))
				{
					rss1 <- max(sum((x$testPredictions[,j]-x$testPredictions$Outcome)^2),1e-10);
					rss2 <- max(sum((x$testPredictions[,i]-x$testPredictions$Outcome)^2),1e-10);
					rss2 <- rss2/rss1;
					pfttable[i-1,j-1] <-  pf(nrow(cp$testPredictions)*(rss2-1.0),1.0,nrow(cp$testPredictions),lower.tail = FALSE);
					fttable[i-1,j-1] <- -log10(max(pfttable[i-1,j-1],0.0001));
				}
			}
			fttable[is.nan(fttable)] <- 0.0;
			colnames(fttable) <- colnames(x$testPredictions)[-1]
			rownames(fttable) <- colnames(x$testPredictions)[-1]
			colnames(pfttable) <- colnames(x$testPredictions)[-1]
			rownames(pfttable) <- colnames(x$testPredictions)[-1]
			par(op);
			par(mfrow = c(1,1),mar = c(2,2,2,2));
			topf <- apply(fttable,1,mean)
			gplots::heatmap.2(fttable[order(topf),],trace = "none",mar = c(5,10),col=rev(heat.colors(8)),Rowv=FALSE,dendrogram = "column",cexRow = 0.65,cexCol = 0.75,srtCol = 25,key.xlab="-log(p)",main = "p(Method A > Method B)",xlab="Method B",ylab="Method A")
			par(op);
		}

		result <- list(metrics = metrics, barPlotsCI = barPlotsCI,metrics_filter=metrics_filter,barPlotsCI_filter=barPlotsCI_filter, minMaxMetrics = minMaxMetrics,fttable=pfttable);


	}
	par(op);
	return (result);
}