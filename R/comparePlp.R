#' function comparePlp
#'
#' @description
#' Compares the performance of two or more patient level prediction models
#' @details
#' The function summarises and plots the performance of the input models for comparison
#' @param models                           A list of plp models
#' @examples
#' modset_llr  <- list(model='lr_lasso',
#'                    param=list(variance =0.001, cohortId=c(1,2), outcomeId=2))
#' class(modset_llr) <- 'modelSettings'
#' model1 <- developModel2(plpData= plpData,
#'                         featureSettings = NULL,
#'                         modelSettings = modset_llr ,
#'                         type='year')
#'                         
#' featSet_gbm <- list(method='wrapperGA', param=list(cohortId=c(1,2), outcomeId=2, varSize=300, iter=25))
#' class(featSet_gbm) <- 'featureSettings'
#' modset_gbm <- list(model='gbm_plp',
#'                    param=list(rsampRate=0.8, ntrees=c(100,150), max_depth=c(2,4,5), cohortId=c(1,2), outcomeId=2))
#' class(modset_gbm) <- 'modelSettings'
#' model2 <- developModel2(plpData= plpData.censor,
#'                          featureSettings = featSet_gbm,
#'                          modelSettings = modset_gbm,
#'                          type='year')
#'                          
#' model3 <- developModel2(plpData= plpData.censor,
#'                          featureSettings = NULL,
#'                          modelSettings = modset_gbm,
#'                          type='year')  
#'                          
#' allModels <- list(model1[[1]], model2[[1]], model3[[1]])   
#' 
#' comparePlp(allModels)                                                                    
#' 
#' @return
#' A table summarising the performance value comparision and plots.
#' @export

comparePlp <- function(models){
  
  # extract model details, training details, cv performance and validation performance
  
  performExtract <- function(x){
    featureSel <- rep('NA',3)
    if(!is.null(x$model$metaData$featureSetting) &!'method'%in%names(x$model$metaData$featureSetting)){
      for(i in 1:min(3,length(x$model$metaData$featureSetting)) ){
          featureSel[i] <- x$model$metaData$featureSetting[[i]]$method
      }
    }
    if(!is.null(x$model$metaData$featureSetting) & 'method'%in%names(x$model$metaData$featureSetting)){
      featureSel[1] <- x$model$metaData$featureSetting$method
    }

    
    res <-  c(x$model$modelSettings$model, 
              paste(names(x$model$modelSettings$modelParameters), x$model$modelSettings$modelParameters,
                     collapse=',', sep='-'),
              featureSel,
              
              #ifelse(!is.null(x$model$modelLoc), x$model$modelLoc,''),
              format(x$time, digits=3),
              x$dataSummary$trainCohort,
              x$dataSummary$trainOutcomeCount,
              x$dataSummary$testCohort,
              x$dataSummary$testOutcomeCount,
              
              ifelse(!is.null(x$model$trainAuc),format(x$model$trainAuc, digits=3),0),
              
              x$evalType,
              format(as.double(x$performance$auc), digits=3), # with lb/up
              format(x$performance$aveP, digits=3),
              format(getTPR(x$performance$roc, FPR=0.05), digits=3), # could get TPR @ 5%/10%  FPR
              format(getTPR(x$performance$roc, FPR=0.1), digits=3),
              eval(x$model$metaData$call$cohortDatabaseSchema),
              x$validationDatabase
    )
    
    names(res) <- c('model','Parameters',
                    'featureSel1', 'featureSel2', 'featureSel3',
                    #'modelLocation', 
                    'time',
                    'train N', 'train Outcome', 'test N', 'test Outcome',
                    'cvAUC',
                    'Evaluation', 'testAUC', 'lower', 'upper', 'AP', 
                    'TPR@5FPR', 'TPR@10FPR', 'TrainDatabase', 'TestDatabase'
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
  
  plot1 <- ggplot2::ggplot(data=plotData, ggplot2::aes(x=FPR, y=TPR, group=model, color=model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::expand_limits(y=0) +
    ggplot2::xlab("FPR") + ggplot2::ylab("TPR") +
    ggplot2::ggtitle("ROC Plot") +
    ggplot2::scale_colour_discrete(name = "Method") +
    ggplot2::geom_abline(intercept = 0, slope = 1, color="grey", 
                linetype="dashed", size=1)  #+
  #ggplot2::geom_line(model)
  
  #plot2 <- models[[i]]$performance$calPlot
  #plot2 <- lapply(models, function(x) x$performance$calPlot + ggtitle("Plant growth"))
  plotCal <- list()
  length(plotCal) <- length(models)
  for(i in 1:length(models)){
    plotCal[[i]] <- models[[i]]$performance$calPlot + ggplot2::ggtitle(paste0('Model: ',i))
  }
  
  
  plot3 <- ggplot2::ggplot(data=plotData, ggplot2::aes(x=TPR, y=PPV, group=model, color=model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::expand_limits(y=0) +
    ggplot2::xlab("TPR") + ggplot2::ylab("PPV") +
    ggplot2::ggtitle("Precision-Recall Plot") +
    ggplot2::theme(legend.position="none")
  
  # Create a table plot
  tbl1 <- gridExtra::tableGrob(result[,c(1,20,21,2:7)])
  tbl2 <- gridExtra::tableGrob(result[,c(1,8:19)])
  # Plot chart and table into one object
  gridExtra::grid.arrange(gridExtra::arrangeGrob(plot1, plot3, ncol=2), 
                          tbl1,tbl2,
                          nrow=3,
                          as.table=TRUE,
                          heights=c(2,0.5,0.5))
  #plot the calibration
  do.call(gridExtra::grid.arrange, c(plotCal, ncol=ceiling(sqrt(length(plotCal)))  ))
  
  result <- list(table=result,
                 plot.roc=plot1,
                 plot.ss = plot3,
                 plot.calibration=plotCal)
  
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


#' varImportance
#'
#' @description
#' Plots the variable importance of the patient level prediction model/s and returns a table of variable importances
#' @details
#' 
#' @param models.list   A list of plp models or simple plp model
#' @return
#' A table containing the variable importances for each model
#'
#' @export
varImportance <- function(models.list){
  if(!'list'%in%class(models.list))
    models.list <- list(models.list)
  
  # get model info for label
  getModDetails <- function(mod){
  return(c(model = mod$model$modelSettings$model,
  database= eval(mod$model$metaData$call$cohortDatabaseSchema),
  validation = mod$evalType,
  parameters = paste(names(mod$model$modelSettings$modelParameters), mod$model$modelSettings$modelParameters,
                     collapse=',', sep='-')))}
  modDetails <- data.frame(modelId=paste0('Model: ', 1:length(models.list)),t(sapply(models.list,getModDetails)))

  
  varImport <- function(mod){
    if(mod$model$modelSettings$model=='lr_lasso'){
      varImp <- data.frame(covariate=names(mod$model$model$coefficients), lr_lasso=mod$model$model$coefficients)
      colnames(varImp)[2] <- paste0('Model: ', mod$id)
      varImp <- varImp[varImp$covariate!='(Intercept)',]
      varImp$covariate <- as.double(as.character(varImp$covariate))
    }
    if(mod$model$modelSettings$model %in% c('lr_enet_plp','gbm_plp','randomForest_plp')){
      varImp <- as.data.frame(mod$model$model@model$variable_importances)[,c('variable','scaled_importance')]
      colnames(varImp) <- c('covariate', paste0('Model: ', mod$id) )
    }
    if(mod$model$modelSettings$model %in% c('nnet')){
      varImp <- caret::varImp(mod$model$model)
      varImp <- data.frame('covariate'=gsub('X', '', rownames(varImp$importance)), 'nnet_plp'=varImp$importance/100)
      colnames(varImp)[2] <- paste0('Model: ', mod$id)
      rownames(varImp) <- NULL
    }
    if(mod$model$modelSettings$model%in%c('knn_plp','svmRadial_plp')){
      warning('Classifier does not currently support variable importance')
      return(NULL)
    }
    
    varImp <- varImp[varImp[,2]!=0, ]
    varImp <- varImp[order(-abs(varImp[,2])),]
    
    return(varImp)
  }
  
  for(i in 1:length(models.list)){models.list[[i]]$id <- i}
  varImps <- lapply(models.list, varImport)
  
  if (sum(sapply(varImps, is.null))>0)
    varImps[-(which(sapply(varImps,is.null),arr.ind=TRUE))]
  
  allImps <- merge(varImps[[1]], ff::as.ram(models.list[[1]]$model$covariateRef)[,c('covariateId','covariateName')], 
                   by.x='covariate', by.y='covariateId', all.y=T)
  if(length(varImps)>1){
    for (i in 2:length(varImps)){
      allImps <- merge(allImps, varImps[[i]], by.x='covariate', by.y='covariate', all.x=T)
    }
  }
  
  # order by overall rank and plot top 100
  allImps[,!colnames(allImps)%in%c('covariate','covariateName')][is.na(allImps[,!colnames(allImps)%in%c('covariate','covariateName')])] <- 0
  if(sum(!colnames(allImps)%in%c('covariate','covariateName'))==1)
    allImps <- allImps[order(-abs(allImps[,!colnames(allImps)%in%c('covariate','covariateName')])),]
  if(sum(!colnames(allImps)%in%c('covariate','covariateName'))>1){
    ranks <- apply(allImps[,!colnames(allImps)%in%c('covariate','covariateName')],
                   2, function(x) rank(-abs(x)))
    ranks.mean <- apply(ranks, 1, mean)
    allImps <- allImps[order(ranks.mean),]
  }
  
  extractName <- function(x){
    if(length(strsplit(x,':')[[1]])>1 && nchar(x)>30)
      x <- paste(strsplit(x,':')[[1]][-1], sep=' ')
    if(nchar(x)>40)
      x <- substring(x,1,40)
    return(x)
  }
  
  allImps$covariateName <- sapply(as.character(allImps$covariateName), extractName)
  # do this to get order by value in ggplot:
  allImps$covariateName=factor(allImps$covariateName,levels=rev(allImps$covariateName))
  
  # plot top 20 overal variables
  melted <- reshape2::melt(allImps[1:20,], id.vars=c("covariate", "covariateName"))
  
  p1<-ggplot2::ggplot(melted,ggplot2::aes(x=factor(covariateName),y=as.double(value),fill=variable))+
    ggplot2::geom_bar(stat="identity", position="dodge") + 
    ggplot2::coord_flip() +
    ggplot2::ylab("Variable Importance") +
    ggplot2::xlab("Variable") +
    ggplot2::ggtitle("Model concensus top 20 variables") +
    ggplot2::scale_fill_discrete(name="ModelId")
  
  #+
  #ggplot2::facet_wrap(variable, nrow = 1, scales = "fixed", shrink = TRUE,drop = TRUE) 
  
  # Create a table plot
  tbl1 <- gridExtra::tableGrob(modDetails)
  # Plot chart and table into one object
  gridExtra::grid.arrange(p1, 
                          tbl1,
                          nrow=2,
                          as.table=TRUE,
                          heights=c(4,1))
  
  result <- list(table=allImps,
                 plot=p1)
  return(result)
}



