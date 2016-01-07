#' comparePlp
#'
#' @description
#' Compares the performance of two or more patient level prediction models
#' @details
#' The function summarises and plots the performance of the input models for comparison
#' @param models                           A list of plp models
#' @return
#' A table summarising the performance value comparision and plots.
#'
#' @export
comparePlp <- function(models){
  
  # extract model details, training details, cv performance and validation performance
  
  performExtract <- function(x){
    featureSel <- rep('NA',3)
    if('wrapper'%in%names(x$model$modelSettings$featureSettings)){
      featureSel <- c('wrapper',
                      x$model$modelSettings$featureSettings$wrapper$method,
                      paste0('Variance:',x$model$modelSettings$featureSettings$wrapper$variance))
    }
    res <-  c(x$model$modelSettings$model, 
              paste0(names(x$model$modelSettings$modelParameters), x$model$modelSettings$modelParameters,
                     collapse=',', sep=':'),
              featureSel,
              
              ifelse(!is.null(x$model$modelLoc), x$model$modelLoc,''),
              
              x$dataSummary$trainCohort,
              x$dataSummary$trainOutcomeCount,
              x$dataSummary$testCohort,
              x$dataSummary$testOutcomeCount,
              
              ifelse(!is.null(x$model$trainAuc),x$model$trainAuc,0),
              
              x$evalType,
              as.double(x$performance$auc), # with lb/up
              x$performance$aveP,
              getTPR(x$performance$roc, FPR=0.05), # could get TPR @ 5%/10%  FPR
              getTPR(x$performance$roc, FPR=0.1)
    )
    
    names(res) <- c('model','Parameters',
                    'featureSection', 'Classifier', 'parameters',
                    'modelLocation', 
                    'trainingCount', 'trainingOutcome', 'testCount', 'testOutcome',
                    'cvAUC',
                    'Evaluation', 'testAUC', 'lower', 'upper', 'AP', 
                    'TPR@5FPR', 'TPR@10FPR'
    )
    return(res) 
  }
  
  result <- t(as.data.frame(lapply(models, performExtract)))
  result <- data.frame(modelId=paste0('model: ',1:nrow(result)), result)
  rownames(result) <- c()
  #plotting:
  
  plotData <- c()
  for(i in 1:length(models)){
  plotData <- rbind(plotData, 
                    data.frame(models[[i]]$performance$precision.recall, 
                    FPR=models[[i]]$performance$roc$FPR,                    
                    model=rep(paste0('model: ',i), 
                              nrow(models[[i]]$performance$precision.recall))))
  }
  
  plot1 <- ggplot2::ggplot(data=plotData, aes(x=FPR, y=TPR, group=model, color=model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::expand_limits(y=0) +
    ggplot2::xlab("FPR") + ggplot2::ylab("TPR") +
    ggplot2::ggtitle("ROC Plot") +
    ggplot2::scale_colour_discrete(name = "Method") +
    geom_abline(intercept = 0, slope = 1, color="grey", 
                linetype="dashed", size=1)  #+
    #ggplot2::geom_line(model)
  
  plot2 <- ggplot2::ggplot(data=plotData, aes(x=TPR, y=PPV, group=model, color=model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::expand_limits(y=0) +
    ggplot2::xlab("TPR") + ggplot2::ylab("PPV") +
    ggplot2::ggtitle("Precision-Recall Plot")
  
  require(gridExtra)

  # Create a table plot
  tbl1 <- gridExtra::tableGrob(result[,c(1:10)])
  tbl2 <- gridExtra::tableGrob(result[,c(1,11:19)])
  # Plot chart and table into one object
  gridExtra::grid.arrange(gridExtra::arrangeGrob(plot1, plot2, ncol=2), 
                          tbl1,tbl2,
                          nrow=3,
                          as.table=TRUE,
                          heights=c(2,0.5,0.5))
  
  return(result)
}


#' getTPR
#'
#' @description
#' Extracts TPR for specified FPR
#' @details
#' 
#' @param roc                           A dataframe containing TPR and FPR at range of thresholds
#' @return
#' TPR at specified FPR
#'
#' @export
getTPR <- function(roc, FPR=0.05){
  return(roc$TPR[which.min(abs(roc$FPR-FPR))])
}

#eval <- list(auc=roc$auc, aveP=aveP,
#roc=roc,
#raw = metrics[,c('TP','FP','TN','FN','FOR','accuracy')],
#precision.recall = metrics[,c('TPR','PPV')],
#F.measure = metrics[,c('Fmeasure')])
#class(eval) <- 'metric.full'