rpredict <-	function(currentModel,DataSet,asFactor=FALSE,classLen=2,...)
{
	fclass <- class(currentModel);
	fclass <- fclass[length(fclass)];
	if (fclass == "numeric") 
	{
		pred <- rep(currentModel,nrow(DataSet));
	}
	else
	{
		pred <- try(predict(currentModel,DataSet))
	}
	if (asFactor && (classLen == 2))
	{
		if ( (fclass == "randomForest") || (fclass == "rpart")  || (fclass == "train.formula")) 
		{
			pred <- try(predict(currentModel,DataSet,type="prob"))[,"1"];
		}
		else
		{
			if (fclass == "svm")
			{
				parameters <- list(...);
				if 	(!is.null(parameters$probability))
				{
					if (parameters$probability)
					{
						pred <- try(predict(currentModel,DataSet,probability = TRUE));
						pred <- attr(pred,"probabilities")[,"1"];
					}
				}
			}
		}
	}
	if (classLen == 2)
	{
		if (fclass == "FRESA_BESS")
		{
			pred <- try(predict(currentModel,DataSet,type = "response"));
		}
		else
		if (fclass == "lm")
		{
			pred <- try(predict(currentModel,DataSet,type = "response"));
		}
	}
	if (inherits(pred, "try-error"))
	{
		pred <- numeric(nrow(DataSet));
	}
	else
	{
		if (class(pred)[1] == "list")
		{
			if (is.null(pred$posterior))
			{
				if (is.null(pred$prob))
				{
					pred <-as.numeric(as.character(pred[[1]]));
				}
				else
				{
					pred <-as.numeric(pred$prob[,2]);
				}
			}
			else
			{
				pred <- pred$posterior;
				if (classLen==2)
				{
					pred <-as.numeric(pred[,2]);
				}
			}
		}
		else
		{
#			if (class(pred)[1] == "factor")
#			{
#				pred <- as.numeric(as.character(pred));
#			}
			if (class(pred)[1] == "array")
			{
				pnames <- colnames(pred);
				pred <- pnames[apply(pred[,,1],1,which.max)];
				pred <- as.numeric(pred);
			}
			if (class(pred)[1] == "matrix") 
			{
				if (ncol(pred)>1)
				{
					pnames <- colnames(pred);
					pred <- pnames[apply(pred,1,which.max)];
					pred <- as.numeric(pred);
				}
				else
				{
					pred <- as.numeric(pred);
				}
			}
		}
	}
	if (classLen == 2)
	{
		if (!is.null(currentModel$type))
		{
			if (currentModel$type== "LOGIT")
			{
				if (class(pred)[1] == "numeric")
				{
					pred[pred == Inf] <- 36;
					pred[pred == -Inf] <- -36;
					if ((min(pred,na.rm=TRUE) < -0.01) || (max(pred,na.rm=TRUE) > 1.01 ))
					{
						pred[pred < -36] <- -36;
						pred[pred > 36] <- 36;
						pred <- 1.0/(1.0 + exp(-pred));
					}
					pred[is.na(pred)] <- 0.5;
					pred[is.nan(pred)] <- 0.5;
					names(pred) <- rownames(DataSet);
				}
			}
		}
	}
	return (pred)
}
	