# default patient level prediction prediction  
predict.plp <- function(plpModel, plpData){
  prediction <- predictProbabilities(plpModel$model, plpData)
  
  return(prediction)
}

# caret model prediction 
predict.caret <- function(plpModel, plpData){
  writeLines('Predicting on validation set...')
  # convert plpData to matrix:
  covariates <- ff::clone(plpData$covariates)
  # if plpModel$features in not null filter non features from covariates

  # now convert into matrix using reshape2::dcast
  writeLines('Converting data into matrix...')
  plpData.matrix <- reshape2::dcast(ff::as.ram(covariates), rowId~covariateId, value.var='covariateValue', fill=0)

  # add people with no covarites:
  ppl <- ff::as.ram(plpData$cohorts$rowId)
  miss.ppl <- ppl[!ppl%in%plpData.matrix$rowId]
  if(length(miss.ppl)>0){
    cov.add <- matrix(rep(0, (ncol(plpData.matrix)-1)*length(miss.ppl)   ), ncol=(ncol(plpData.matrix)-1))
    cov.add <- data.frame(miss.ppl,cov.add)
    colnames(cov.add) <- colnames(plpData.matrix)
    plpData.matrix <-  rbind(plpData.matrix, cov.add)
  }
  
  # convert the columns to have X
  colnames(plpData.matrix)[!colnames(plpData.matrix)%in%c('rowId')] <- paste0('X',colnames(plpData.matrix)[!colnames(plpData.matrix)%in%c('rowId')])
  
  # remove features not in model:
  if(!"ksvm"%in%class(plpModel$model$finalModel))
    modelCoef <-plpModel$model$finalModel$xNames
  if("ksvm"%in%class(plpModel$model$finalModel))
    modelCoef <- colnames(plpData.matrix) #ifelse(length(plpModel$model$coefnames)==0, colnames(plpData.matrix), plpModel$model$coefnames)
  
  writeLines(paste0('model: ', length(modelCoef), ' -- data: ', ncol(plpData.matrix)))
  plpData.matrix <- plpData.matrix[,colnames(plpData.matrix)%in%c(modelCoef,'rowId')]
  
  missCoef <- modelCoef[!modelCoef%in%colnames(plpData.matrix)]
  if(length( missCoef )>0){
    writeLines(paste0('missing these covariates: ', paste(missCoef, collapse=',', sep=',')))
    addval <- matrix(rep(0, length( missCoef )*nrow(plpData.matrix)), ncol=length( missCoef))
    colnames(addval) <- missCoef
    plpData.matrix <- cbind(plpData.matrix, addval)
  }
  writeLines(paste0('model: ', length(modelCoef), ' -- data: ', ncol(plpData.matrix[,-1])))
  
  
    #writeLines(paste(plpModel$model$finalModel$xNames, collapse=',',sep=','))
  
  #writeLines(paste(colnames((plpData.matrix[,-1])), collapse=',',sep=','))
  
  resp <- 'raw'
  if(class(plpModel$model$finalModel)=="ksvm")
    resp <- 'probabilities'
  prediction <- predict(plpModel$model$finalModel, plpData.matrix[,-1], type=resp)
  ##prediction <- ffdf(rowId = as.ff(plpData.matrix$rowId), value = as.ff(prediction))
  prediction <- data.frame(rowId = plpData.matrix$rowId, value = prediction[,1])
  ##writeLines(paste(colnames(prediction), sep='', collapse='_'))
  writeLines('Prediction complete...')
  # check whether I need to add rowId name and column.
  prediction <- merge(ff::as.ram(plpData$cohorts), prediction, by='rowId', all.x=T)
  prediction$value[is.na(prediction$value)] <- 0
   
  return(prediction)
}





# default h2o prediction
predict.h2o <- function(plpModel, plpData){
  covariates <- ff::clone(plpData$covariates)

  # now convert into matrix using reshape2::dcast
  cov <- reshape2::dcast(ff::as.ram(covariates), rowId~covariateId, value.var='covariateValue', fill=0)
  
  # add people with no covarites:
  ppl <- ff::as.ram(plpData$cohorts$rowId)
  miss.ppl <- ppl[!ppl%in%cov$rowId]
  if(length(miss.ppl)>0){
    cov.add <- matrix(rep(0, (ncol(cov)-1)*length(miss.ppl)   ), ncol=(ncol(cov)-1))
    cov.add <- data.frame(miss.ppl,cov.add)
    colnames(cov.add) <- colnames(cov)
    cov <-  rbind(cov, cov.add)
  }
  
  cov.h2o <-h2o::as.h2o(cov)
  value <- h2o::h2o.predict(plpModel$model, cov.h2o)
  #writeLines(paste(colnames(value), sep='-',collapse='-'))
  #writeLines(paste(as.data.frame(value)[1,], sep='-',collapse='-'))
  pred <- data.frame(rowId=cov$rowId, value=as.data.frame(value)[,3])
  #writeLines(paste(pred[1,], sep='-',collapse='-'))
  #writeLines(paste(colnames(ff::as.ram(plpData$cohorts)), sep='-',collapse='-'))
  
  prediction <- merge(ff::as.ram(plpData$cohorts), pred, by='rowId', all.x=T)
  
  #writeLines(paste(prediction[1,], sep='-',collapse='-'))
  return(prediction)
}



predict.knn <- function(plpData, plpModel){
  
  prediction <- BigKnn::predictKnn(covariates = plpData$covariates,
                                   indexFolder = plpModel$model,
                                   k = plpModel$modelSettings$modelParameters$k,
                                   weighted = TRUE)
  
  # return the cohorts as a data frame with the prediction added as 
  # a new column with the column name 'value'
  prediction <- merge(ff::as.ram(plpData$cohorts), prediction, by='rowId', 
                      all.x=T, fill=0)
  prediction <- prediction[!is.na(prediction$value),]
  
  return(prediction)
}
