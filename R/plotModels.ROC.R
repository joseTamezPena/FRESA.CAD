plotModels.ROC <-
function(modelPredictions,number.of.models=0,specificities=c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),theCVfolds=1,predictor="Prediction",cex=1.0,thr=NULL,...) 
{

#	op <- par(no.readonly=TRUE);

#	par(mfrow=c(1,1),pty='s',cex=cex);
	par(pty='s',cex=cex);
	rocadded = 0;
	auclist <- vector()
	sumSen <- NULL;
	blindSen <- NULL;
	ensemblePrediction <- NULL;
	ensemble.auc <- NULL;
	rout <- NULL;
	F1 <- NULL;
	
	if ((ncol(modelPredictions) == 3) && (number.of.models == 0))
	{
		number.of.models <- max(modelPredictions[,"Model"]);
	}
	number.of.runs=number.of.models;

	if ((theCVfolds>1) || (number.of.models>0))
	{
		if (number.of.runs == 0)
		{
			number.of.runs=max(modelPredictions[,"Model"]) %/% theCVfolds;
		}
		for (n in 1:number.of.runs)
		{
			mm = n;
			blindmodel <- modelPredictions[which(((modelPredictions[,"Model"]-1) %/% theCVfolds) + 1  == mm),];
			if ( (sum(blindmodel[,"Outcome"]==1) > 3) && (sum(blindmodel[,"Outcome"]==0) > 3))
			{
				auclist <- append(auclist,pROC::roc(as.vector(blindmodel[,"Outcome"]),blindmodel[,predictor],auc=TRUE,plot=TRUE,col="lightgray",lty=4,lwd=1,direction="<",quiet = TRUE)$auc)
				par(new=TRUE)
				sen <- pROC::roc(as.vector(blindmodel[,"Outcome"]),blindmodel[,predictor],auc=TRUE,plot=FALSE,ci=TRUE,progress='none',of='se',specificities=specificities,boot.n=100,smooth=FALSE,lty=3,lwd=1,direction="<",quiet = TRUE)$ci[,2]
				if (n == 1) 
				{
					blindSen <- sen;
				}
				else
				{
					blindSen <- rbind(sen,blindSen);
				}
				rocadded = rocadded +1;
			}
		}
		psta <- boxplot(modelPredictions[,predictor]~rownames(modelPredictions),plot=FALSE)
		outcomesta <- boxplot(modelPredictions[,"Outcome"]~rownames(modelPredictions),plot=FALSE)
		rout <- pROC::roc(outcomesta$stats[3,],psta$stats[3,],col="black",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3,lwd=3,direction="<",quiet = TRUE,...)
		ensemble.auc <- rout$auc
		par(new=TRUE)
		auc1 <- pROC::roc(as.vector(modelPredictions[,"Outcome"]),modelPredictions[,predictor],col="darkblue",auc=TRUE,plot=TRUE,smooth=FALSE,direction="<",quiet = TRUE,...)$auc

		ensemblePrediction <- cbind(outcomesta$stats[3,],psta$stats[3,]);
		rownames(ensemblePrediction) <- psta$names
		colnames(ensemblePrediction) <- c("Outcome",predictor);
		if (is.null(thr))
		{
			thres <- 0.5
			if (min(psta$stats[3,])<0) 
			{
				thres <- 0;
			}				
		}
		else
		{
			thres <- thr
		}
		dtable <- table(psta$stats[3,]<thres,1-outcomesta$stats[3,])
		ley.names <- c(paste("Coherence (",sprintf("%.3f",auc1),")"));
		ley.colors <- c("darkblue");
		ley.lty <- c(1);
		ley.names <- append(ley.names,paste("Ensemble (",sprintf("%.3f",ensemble.auc),")"));
		ley.colors <- append(ley.colors,"black");
		ley.lty <- append(ley.lty,3);
	}
	else
	{
		if (class(modelPredictions)[1] == "data.frame")
		{ 
			if (!is.null(modelPredictions$Model)) modelPredictions$Model <- NULL;
		}
		lastcol <- ncol(modelPredictions);
		eblindmodel <- NULL;
		if (lastcol>2)
		{
			for (n in 2:lastcol)
			{
				mm = n;
				blindmodel <- modelPredictions[,c(1,n)];
				colnames(blindmodel) <- c("Outcome",predictor);
				if ( (sum(blindmodel[,"Outcome"]==1) > 3) && (sum(blindmodel[,"Outcome"]==0) > 3))
				{
					auclist <- append(auclist,pROC::roc(blindmodel[,"Outcome"],blindmodel[,predictor],auc=TRUE,plot=TRUE,col="lightgray",lty=4,lwd=1,direction="<",quiet = TRUE)$auc)
					par(new=TRUE)
					sen <- pROC::roc(as.vector(blindmodel[,"Outcome"]),blindmodel[,predictor],auc=TRUE,plot=FALSE,ci=TRUE,progress='none',of='se',specificities=specificities,boot.n=100,smooth=FALSE,lty=3,lwd=1,direction="<",quiet = TRUE)$ci[,2]
					if (n == 1) 
					{
						blindSen <- sen;
					}
					else
					{
						blindSen <- rbind(sen,blindSen);
					}
					rocadded = rocadded +1;
				}
				eblindmodel <- rbind(eblindmodel,blindmodel);
			}
			psta <- modelPredictions[,2:lastcol];
			if (ncol(psta)>1)
			{
				psta <- rowMedians(psta,na.rm=TRUE);
			}
			outcomesta <- modelPredictions[,1];
			rout <- pROC::roc(outcomesta,psta,col="black",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3,lwd=3,direction="<",quiet = TRUE,...)
			ensemble.auc <- rout$auc
			par(new=TRUE)
			auc1 <- pROC::roc(as.vector(eblindmodel[,1]),eblindmodel[,2],col="darkblue",auc=TRUE,plot=TRUE,smooth=FALSE,direction="<",quiet = TRUE,...)$auc

			ensemblePrediction <- cbind(outcomesta,psta);
			rownames(ensemblePrediction) <- rownames(modelPredictions);
			colnames(ensemblePrediction) <- c("Outcome",predictor);
			if (is.null(thr))
			{
				thres <- 0.5
				if (min(psta$stats[3,])<0) 
				{
					thres <- 0;
				}				
			}
			else
			{
				thres <- thr
			}
			dtable <- table(psta<thres,1-outcomesta);
			ley.names <- c(paste("Coherence (",sprintf("%.3f",auc1),")"));
			ley.colors <- c("darkblue");
			ley.lty <- c(1);
			ley.names <- append(ley.names,paste("Ensemble (",sprintf("%.3f",ensemble.auc),")"));
			ley.colors <- append(ley.colors,"black");
			ley.lty <- append(ley.lty,3);
		}
		else
		{
			rout <- pROC::roc(as.vector(modelPredictions[,1]),modelPredictions[,2],direction="<",quiet = TRUE);
			specificities=seq(0, 1, .05);	
			if (nrow(modelPredictions) < 5000)
			{
				ci.sp.obj <- pROC::ci.sp(rout , sensitivities=seq(0, 1, .05), boot.n=100,progress= 'none',quiet = TRUE)
				blindSen <- pROC::ci.se(rout , specificities=seq(0, 1, .05), boot.n=100,progress= 'none',quiet = TRUE)
				pROC::plot.roc(rout,grid=c(0.1, 0.1),grid.col=c("gray", "gray"),print.auc=FALSE,quiet = TRUE,...) 
				plot(ci.sp.obj, type="s", col="gray")
				plot(blindSen, type="s", col="light gray")
			}
			else
			{
				pROC::plot.roc(rout,grid=c(0.1, 0.1),grid.col=c("gray", "gray"),print.auc=FALSE,quiet = TRUE,...) 
			}
			auclist <- rout$auc;
			if (is.null(thr))
			{
				thres <- 0.5
				if (min(modelPredictions[,2])<0) 
				{
					thres = 0;
				}				
			}
			else
			{
				thres <- thr
			}
			dtable <- table(modelPredictions[,2]<thres,1-modelPredictions[,1]);
			ley.names <- c(paste("AUC: ",sprintf("%.3f",rout$auc)));
			ley.colors <- c("black");
			ley.lty <- c(1);
		}
	}
	if (rocadded>1)
	{
		par(new=TRUE)
		boxplot(blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		par(new=TRUE)
	}

	auc = 0;
	enauc = -1;
	if (rocadded>1)
	{
#		boxplot(blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		sumSen <- colMeans(blindSen,na.rm = TRUE);
		sennames <- names(sumSen);
		sumSen <- append(0,sumSen);
		sumSen <- append(sumSen,1);
		sennames <- append("1",sennames);
		sennames <- append(sennames,"0");
		names(sumSen) <- sennames;
		spevalues <- as.numeric(names(sumSen));
		lines(spevalues,sumSen,col="red",lwd=1.0,lty=2);
		for (i in 2:length(spevalues))
		{
			auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
		}
		ley.names <- append(ley.names,c("Fold ROCs",paste("Mean Sensitivities(",sprintf("%.3f",auc),")")));
		ley.colors <- append(ley.colors,c("lightgray","red"));
		ley.lty <- append(ley.lty,c(4,2));
	}


	if (nrow(dtable)==ncol(dtable) & (ncol(dtable)>1))
	{
		colnames(dtable) <- c("O(+)","O(-)")
		rownames(dtable) <- c("T(+)","T(-)")
		Sen=dtable[1,1]/(dtable[1,1]+dtable[2,1])
		Spe=dtable[2,2]/(dtable[1,2]+dtable[2,2])
		enauc = 0.5*(Sen+Spe)	

		lines(c(1,Spe,0),c(0,Sen,1),col="green",lwd=1.0,lty=1);
	
		ley.names <- append(ley.names,paste("Class AUC (",sprintf("%.3f",enauc),")"));
		ley.colors <- append(ley.colors,"green");
		ley.lty <- append(ley.lty,1);

		legend(0.7,0.20, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n",cex=0.9*cex)
		Acc=(dtable[1,1]+dtable[2,2])/(dtable[1,1]+dtable[2,2]+dtable[1,2]+dtable[2,1])
		F1=(2*dtable[1,1])/(2*dtable[1,1]+dtable[1,2]+dtable[2,1])
		x=1.025
		y=1.025
		text(x,y,paste("(TPR=",sprintf("%.3f",Sen),",TNR=",sprintf("%.3f",Spe),",ACC=",sprintf("%.3f",Acc),",F1=",sprintf("%.3f",F1),",AUC=",sprintf("%.3f",enauc),")"),adj = c(0,1),cex=0.7*cex,col="dark green")
		par(new=TRUE,plt=c(0.6,0.8,0.37,0.57),pty='s',cex=0.8*cex)
		plot(t(dtable),main="Confusion Matrix",ylab="Test",xlab="Outcome",cex=0.8*cex)
	}
	else
	{
		legend(0.7,0.20, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n",cex=0.9*cex)
	}
#	par(op);
	
	if (nrow(dtable) == 1)
	{
		if ((rownames(dtable) == "0") || (rownames(dtable) == "FALSE"))
		{
			dtable <- rbind(dtable,c(0,0));
		}
		else
		{
			dtable <- rbind(c(0,0),dtable);
		}
		rownames(dtable) <- c("0","1")
	}
	result <- list(ROC.AUCs=auclist,
	mean.sensitivities=sumSen,
	model.sensitivities=blindSen,
	specificities=specificities,
	senAUC=auc,
	ensemblePrediction=ensemblePrediction,
	predictionTable=dtable,
	ensemble.auc=ensemble.auc,
	clasification.auc=enauc,
	roc.predictor=rout,
	F1=F1
	)
	return (result)
}
