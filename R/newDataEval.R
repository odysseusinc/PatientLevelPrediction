#' Make predictions using model on new data
#'
#' @details
#' Computes the performance of the model on new data
#'
#' @param newData               An object of type \code{plpData}.
#' @param plpModel              A model returned by calling developModel2()

#' @export
newDataEval <- function(newData, plpModel){
  if(!'overall.plpModel' %in% class(plpmodel)){
    stop('You did not enter a correct model: try adding [[1]] to your input, e.g., mod[[1]] instead of mod')
  }
  start.all <-  Sys.time()
  # clone newData
  data <- list(outcomes =ff::clone(newData$outcomes),
                  cohorts =ff::clone(newData$cohorts),
                  covariates =ff::clone(newData$covariates),
                  exclude =ff::clone(newData$exclude),
                  covariateRef=ff::clone(newData$covariateRef),
                  metaData=newData$metaData
  )
  
  # apply feature selection
  # apply the feature transformations
  if('list'%in%class(plpModel$model$metaData$featureSettings) && !'transform'%in%names(plpModel$model$metaData$featureSettings)){
    for(x in plpModel$model$metaData$featureSettings){
      data <- do.call(x$transform, list(plpData2=data))
    }
  }
  if('transform'%in%names(plpModel$model$metaData$featureSettings) ){
    data <- do.call(plpModel$model$metaData$featureSettings$transform, list(plpData2=data))
  }
  
  # do the predction on the new data
  if(class(plpModel$model)=='plpModel'){
    # extract the classifier type
    pred <- paste0('predict.',attr(plpModel$model, 'type'))
    prediction <-  do.call(pred, list(plpModel=plpModel$model, plpData=data))
  }
  
  attr(prediction, "modelType") <- "logistic"
  attr(prediction, "outcomeId") <- plpModel$model$modelSettings$outcomeId
  
  
  # calculate metrics
  performance <- evaluatePlp(prediction, data)
  writeLines('3) Performance calculated')
  
  comp <- Sys.time() - start.all
 results <- list(prediction=prediction,
                 performance=performance,
                 time=plpModel$time,
                 predictionTime=comp,
                 trainDatabase= eval(plpModel$model$metaData$call$cohortDatabaseSchema),
                 validationDatabase= eval(newData$metaData$call$cdmDatabaseSchema),
                 model= list(metaData=plpModel$model$metaData,
                             trainAUC= plpModel$model$trainAUC,
                             dataSummary=plpModel$dataSummary, # update test with newdata counts!
                             modelSettings=plpModel$model$modelSettings)
 )
  
}