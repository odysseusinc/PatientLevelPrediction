# @file developModel.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
#
# This file is part of PatientLevelPrediction
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' developModel - Train and elavuate the model
#'
#' @description
#' Trains various machine learning models and elavulates them
#' @details
#' Users can define a risk period of interest for the prediction of the outcome relative to index or use
#' the cohprt dates.  The user can then specify whether they wish to exclude patients who are not observed
#' during the whole risk period, cohort period or experienced the outcome prior to the risk period.
#'
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param modelSettings                    A list defining the model to train, the cohortId, the outcomeId,
#'                                         preprocessing options and model parameters.
#'
#'                                         The potential models are:
#'                                         Logistic regression with lasso regularisation: 'lr-lasso'
#'                                         Random forest: 'rf'
#'                                         Gradient boosting machineL 'gbm'
#'                                         Logistic regression with elastic new regularisation: 'lr-enet'
#'
#'                                         The potential preprocessing settings are:
#'                                         Wrapper feature selection using lasso logistic regression: 'lr-lasso'
#'                                         Selecting only the condition/drug era features: 'allEra'
#'
#'                                         The parameter setting depend on the model-
#'                                         model='lr-lasso' has the parameter val specificying the initial variance
#'                                         model='rf' has the parameters ntree, max_depth and mtry
#'                                         model='gbm' has the parameters ntree, bal (class balance), nrowSample (fraction of training data people to use per tree),
#'                                                ncolSample (fraction of training data features to include into each tree)
#'                                         model='lr-enet' has the parameters alpha, ...
#'
#'                                         example: list(model='gbm', cohortId=NULL, outcomeId=c(1,2),
#'                                                       preprocess='lr-lasso', param=list(ntree=50, bal=F, nrowsample=0.6))
#'                                         Would train a gradient boosting machine to predict the outcomes 1 and 2
#'                                         using the features selected by logistic regression with lasso regularisation
#'                                         with the model settings ntree 50, no class label balance and 0.6 training data rows
#'                                         used per tree.
#'
#' @param validationFraction               The fraction of the data to be used as the validation set in the patient
#'                                         split evaluation.
#' @param fileLoc                          The path to the directory where the models will be saved
#' @param type                             A subset of c('year','both','patient') specifying the type of evaluation used.
#'                                         'year' splits the date prior to 2013 into the training set and post 2013 into the test set
#'                                         'both' splits the data by the year 2013 but removes any patient in the test set from the training set
#'                                         'patient' splits the data into test (1-validationFraction of the data) and
#'                                         train (validationFraction of the data) sets.  The split is stratified by the class label.

#'
#' @return
#' An object containing the model or location where the model is save, the data selection settings, the preprocessing
#' and training settings as well as various performance measures obtained by the model.
#' \describe{
#' \item{model}{The trained prediction model}
#' \item{dataSummary}{A list detailing the size of the train/test sets and outcome prevalence}
#' \item{evalType}{The type of evaluation that was performed}
#' \item{prediction}{The model prediction on the test set }
#' \item{performance}{A list detailing the performance of the model}}
#'
#'
#' @export
developModel <- function(plpData,
                         modelSettings=list(model='lr-lasso', param=NULL,
                                            cohortId=NULL, outcomeId=2),
                         featureSettings = list(analysisSelector=NULL,
                                                covariateSelector=NULL,
                                                wrapper=list(method='lr-lasso', variance=0.01),
                                                matrixFactor=F),
                         validationFraction=0.2, fileLoc=file.path(getwd(),'models'),
                         type = c('year','both','patient')
){


  # extract censoring details
  # train 3 models 1) will test/train on old years and validate on last year of data (check using describePlpData)
  #                2) will also split on people (remove any in last year from training)
  #                3) will test/train on n% of people and validate on (100-n)%
  if(!dir.exists(fileLoc)){dir.create(fileLoc)}
  settings <- list()
  length(settings) <- length(type)
  results <- list()
  length(results) <- length(type)



  for (i in 1:length(type)){
    # split data
    data <- spliter(plpData, type[i], validationFraction)
    settings[[i]] <- list(model = modelSettings$model,
                          param = modelSettings$param,
                          featureSettings = featureSettings,
                          cohortId = modelSettings$cohortId,
                          outcomeId = modelSettings$outcomeId,
                          data=data[[1]], loc=fileLoc)

    dataSummary= list(trainCohort =nrow(data[[1]]$cohorts),
                      trainOutcomeCount =nrow(data[[1]]$outcomes),
                      testCohort =nrow(data[[2]]$cohorts),
                      testOutcomeCount =nrow(data[[2]]$outcomes),
                      covariateCount = length(unique(data[[1]]$covariates$covariateId))
    )

    # train model
    model <- do.call(fitPlp, settings[[i]])
    writeLines('1: Model Trained')

    prediction <- NULL
    performance <- NULL
    if(model!='no cov'){
      # do prediction
      prediction <- predictPlp(model, data[[2]])
      writeLines('2) Prediction Calculated')
      # calculate metrics
      performance <- evaluatePlp(prediction, data[[2]])
      writeLines('3) Performance calculated')
    }

    results[[i]] <- list(model=model, dataSummary=dataSummary,
                         evalType=type[i],
                         prediction=prediction,
                         performance=performance)
  }

  #return: trainingSummary, valSummary, settings

  # return table of form:
  #                   train: auc, pAuc, p10, map, ... val: auc, pAuc, p10, map, ...
  # time-split:
  # person-split
  # timeperson-split:
  return(results)

}






#' fitModel
#'
#' @description
#' Train various models using a default parameter gird search or user specified parameters
#'
#' @details
#' The user can define the machine learning model to train (regularised logistic regression, random forest,
#' gradient boosting machine, neural network and )
#' @param model                            A character string specifiying the machine learning model to train.
#'
#'                                         The potential models are:
#'                                         \description{
#'                                         \item{'lr-lasso'}{Logistic regression with lasso regularisation}
#'                                         \item{'rf'}{Random forest}
#'                                         \item{'gbm'}{Gradient boosting machines}
#'                                         \item{'lr-enet'}{Logistic regression with elastic new regularisation}
#'                                         }
#' @param param                            The parameter setting depend on the model-
#'                                         model='lr-lasso' has the parameter val specificying the initial variance
#'                                         model='rf' has the parameters ntree, max_depth and mtry
#'                                         model='gbm' has the parameters ntree, bal (class balance), nrowSample (fraction of training data people to use per tree),
#'                                                ncolSample (fraction of training data features to include into each tree)
#'                                         model='lr-enet' has the parameters alpha, ...
#'
#'
#' @param featureSettings                  A list containing any parameters for the feature selection/engineering.
#'                                         \description{
#'                                         \item{...}{The initial variance setting for the lasso logistic regression}
#'                                         }
#' @param outcomeId                        The outcomeId the user is aiming to predict
#' @param cohortId                         The id of the cohort being used.  Default is NULL.
#' @param data                             An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param loc                              The path to the directory where the model will be saved

#'
#' @return
#' An object of class plpModel
#' \describe{
#' \item{model}{The trained prediction model}
#' \item{modelLoc}{The path to where the model is saved}
#' \item{trainAuc}{The AUC obtained on the training set}
#' \item{trainCalibration}{The calibration obtained on the training set}
#' \item{modelSettings}{A list specifiying the model, preprocessing, outcomeId and cohortId}
#' \item{features}{The features used by the model}}
#'
#'

#' @export
fitPlp <- function(model,param, featureSettings, outcomeId, cohortId, data, loc){
  plpData <- list(outcomes =ff::clone(data$outcomes),
                  cohorts =ff::clone(data$cohorts),
                  covariates =ff::clone(data$covariates),
                  exclude =ff::clone(data$exclude),
                  covariateRef=ff::clone(data$covariateRef),
                  metaData=plpData$metaData)

  # set saving locations
  saveLoc <- file.path(loc, paste0(model,'_',makeRandomString(),'.rds'))
  modelLoc <- file.path(loc, paste0(model,'_',makeRandomString()))

  # ================== STANDARD MODEL ================================
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(model=='lr-lasso'){
    val <- 0.003
    if(!is.null(param$val )) val <- param$val
    writeLines(paste0('Training ',model, ' model'))
    start <- Sys.time()
    modelTrained <- fitPredictiveModel(plpData = plpData,
                                       modelType = "logistic",
                                       removeDropoutsForLr = F,
                                       cohortId = cohortId,
                                       outcomeId = outcomeId,
                                       prior = createPrior("laplace",
                                                           exclude = c(0),
                                                           variance = val))
    comp <- Sys.time() - start
    writeLines(paste0('Model ',model, ' trained - took:',  format(comp, digits=3)))
    writeLines('Saving model...')
    result <- list(model = modelTrained,
                   modelLoc = NULL,    # did I actually save this!?
                   modelType = "logistic",
                   trainAuc = NULL,
                   trainCalibration=NULL,
                   modelSettings = list(model=model, param= param,
                                        featureSettings= NULL,
                                        outcomeId=outcomeId, cohortId=cohortId),
                   metadata = plpData$metaData
    )
    class(result) <- 'plpModel'
    saveRDS(result, saveLoc)
  }

  # ================== EXTRA MODELS ====================================
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(model%in%c('randomForest','gbm','lr-enet', 'nnet','naiveBayes')  ){
    writeLines(paste0('Training ',model, ' model'))
    start.main <- Sys.time()

    # ================== FEATURE REDUCTION =============================
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #perform the feature reducion:
    ##plpData <- featureReducer(plpData, analysisSelector=featureSettings$analysisSelector,
    ##                          covariateSelector=featureSettings$covariateSelector,
    ##                          wrapper=featureSettings$wrapper,
    ##                          cohortId = featureSettings$cohortId, outcomeId = featureSettings$outcomeId,
    ##                          matrixFactor=featureSettings$matrixFactor)
    featureSettings$plpData <- plpData
    plpData <- do.call(featureReducer,featureSettings )
    if(is.null(plpData$covariates))
      return('no cov')
    # now convert into h2o matrix
    writeLines('Converting sparse data into matrix...')
    start <- Sys.time()
    cov <- reshape2::dcast(ff::as.ram(plpData$covariates), rowId~covariateId, value.var='covariateValue', fill=0)
    allData <- merge(cov, ff::as.ram(plpData$outcomes[,c('rowId','outcomeCount')]), by='rowId', all.x=T)
    allData$outcomeCount[is.na(allData$outcomeCount)] <- 0
    allData$outcomeCount <- as.factor(allData$outcomeCount)
    writeLines(paste0('Conversion tooK:', format(Sys.time()-start, digits=3)))

    saveLoc2 <- file.path(loc, paste0(model,'_',makeRandomString(),'.csv'))
    writeLines(paste0('Saving matrix to ', saveLoc2))
    #writeLines(paste0(allData[1,1]))
    start <- Sys.time()
    write.csv(allData, saveLoc2, row.names=F)
    rm(allData)
    writeLines(paste0('Saving took: ', format(Sys.time()-start, digits=3)))

    # ================== CARET MODELS ====================================
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(model%in%c('nnet','naiveBayes')){
      require(caret)
      allData <- read.csv( saveLoc2)
      newLab <- rep('no', nrow(allData))
      newLab[allData$outcomeCount==1] <-  'yes'
      allData$outcomeCount <- newLab
      if(model=='nnet'){

        size <- c(2,20,50)
        if(!is.null(param$size))
          size <- param$size
        decay <- c(0,0.1, 0.05)
        if(!is.null(param$size))
          decay <- param$decay
        weights <- rep(1, nrow(allData))
        weights[allData$outcomeCount=='yes'] <- sum(allData$outcomeCount=='no')/sum(allData$outcomeCount=='yes')

        tuneGrid <- expand.grid(size=size, decay=decay)
        fitControl <- caret::trainControl(method = "repeatedcv", number = 3,repeats = 1,
                                   verboseIter = FALSE,classProbs = TRUE,
                                   summaryFunction=twoClassSummary)

        model <- caret::train(x=allData[,!colnames(allData)%in%c('outcomeCount','rowId')],
                              y=as.factor(allData$outcomeCount),
                              method = "nnet",
                              #preProcess = NULL,
                              weights = weights,
                              metric = 'ROC',
                              maximize = TRUE,
                              trControl = fitControl,
                              tuneGrid = tuneGrid,
                              maxit=500,MaxNWts=20000)

        param.string <- paste(paste0(c('size','decay'),':',model$results[which.max(model$results$ROC),c('size','decay')]), collapse=',')
        writeLines(paste0('Neural Network with parameters ',param.string,' obtained AUC: ', model$results$ROC[which.max(model$results$ROC)]))

        param.best <- model$results[which.max(model$results$ROC),]

        result <- list(model = model,
                       modelLoc = NULL,
                       trainAuc = model$results$ROC[which.max(model$results$ROC)],
                       trainCalibration= NULL,
                       modelSettings = list(model=model,modelParameters=param.best,
                                            featureSettings=featureSettings,
                                            outcomeId=outcomeId, cohortId=cohortId),
                       metaData = plpData$metaData
        )
        class(result) <- 'plpModel'
        saveRDS(result, saveLoc)
        writeLines(paste0('Saved model to : ',saveLoc))

      }

      if(model=='naiveBayes'){

        fitControl <- caret::trainControl(method = "repeatedcv", number = 3,repeats = 1,
                                   verboseIter = FALSE,returnResamp = "all",classProbs = TRUE,
                                   summaryFunction=twoClassSummary)

        model <- caret::train(x=allData[,!colnames(allData)%in%c('outcomeCount','rowId')],
                              y=as.factor(allData$outcomeCount),
                              method = "nb",
                              preProcess = NULL,
                              weights = NULL,
                              metric = 'ROC',
                              maximize = TRUE,
                              trControl = fitControl)

        writeLines(paste0('Naive Bayes obtained AUC: ', model$results$ROC))

        result <- list(model = model,
                       modelLoc = NULL,
                       trainAuc = model$results$ROC,
                       trainCalibration= NULL,
                       modelSettings = list(model=model,modelParameters=NULL,
                                            featureSettings=featureSettings,
                                            outcomeId=outcomeId, cohortId=cohortId),
                       metaData = plpData$metaData
        )
        class(result) <- 'plpModel'
        saveRDS(result, saveLoc)
        writeLines(paste0('Saved model to : ',saveLoc))

      }


    }

    # ================== H2O MODELS ====================================
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(model%in%c('randomForest','gbm','lr-enet')){
      require(h2o)
      # now load into h2o:
      writeLines(paste0('Loading into h2o data frame' ))
      start <- Sys.time()
      h2oData <- h2o::h2o.importFile(path = saveLoc2, header = T)
      h2oData$outcomeCount <- h2o::as.factor(1*(h2oData$outcomeCount>0))
      writeLines(paste0('Loading took: ', format(Sys.time()-start, digits=3)))


      # now apply cross val to fit the model:
      writeLines(paste0('Training model...' ))
      start <- Sys.time()
      if(model =='randomForest'){

        #bal=bal, sample_rate=rsampRate, ntrees=ntrees, mtries=mtries,nbins=nbins, max_depth=max_depth, min_rows=min_rows
        rfTrainer <- function(sample_rate=0.5,mtries=-1, ntrees=50, bal=F,
                              nbins=20, max_depth=4, min_rows=20){
          modelTrained <- h2o::h2o.randomForest(x=2:(ncol(h2oData)-1) , y=ncol(h2oData),
                                           training_frame = h2oData, sample_rate=sample_rate,
                                           mtries=mtries, nbins=nbins,
                                           ntrees = ntrees, max_depth=max_depth,
                                           balance_classes = bal, nfolds = 3
          )
          param.string <- paste(paste0(names(as.list(match.call()) ),':',as.list(match.call()))[-1], collapse=',')
          writeLines(paste0('Random forest model with params: ',param.string,' obtained AUC: ',modelTrained@model$cross_validation_metrics@metrics$AUC))
          auc <- modelTrained@model$cross_validation_metrics@metrics$AUC
          model <- modelTrained
          return(list(auc=auc, model=model))
        }

        # default grid search:
        if(!is.null(param))
          param <- do.call(paramSettings, param)
        if(is.null(param))
          param <- split(expand.grid(bal=c(T,F), mtries=c(-1), ntrees=c(20,50,100)), 1:6)
        
        res <- lapply(param, function(x) do.call(rfTrainer, x ))
        modelTrained <- res[[which.max(unlist(lapply(res, function(x) x$auc)))]]$model
        param.best <- param[[which.max(unlist(lapply(res, function(x) x$auc)))]]


      }
      if(model =='lr-enet'){
        glmTrainer <- function(alpha=0.5, lambda=0.000001, lambda_search=T,
                               lambda_min_ratio = 1/1000,nlambdas = 100){
          modelTrained <- h2o::h2o.glm(x=2:(ncol(h2oData)-1) , y=ncol(h2oData),
                                  training_frame = h2oData, family= "binomial",
                                  alpha=alpha,
                                  lambda = lambda, lambda_search = lambda_search,
                                  lambda_min_ratio = lambda_min_ratio, nlambdas = nlambdas)
          param.string <- paste(paste0(names(as.list(match.call()) ),':',as.list(match.call()))[-1], collapse=',')
          writeLines(paste0('Elastic net logistic regression model with params: ',param.string,' obtained AUC: ',modelTrained@model$training_metrics@metrics$AUC))
          auc <- modelTrained@model$training_metrics@metrics$AUC
          model <- modelTrained
          return(list(auc=auc, model=model))
        }

        # default grid search:
        if(!is.null(param))
          param <- do.call(paramSettings, param)
        if(is.null(param))
          param <- split(expand.grid(lamba=c(0.000001), alpha=c(0,0.2,0.5)), 1:3)


        res <- lapply(param, function(x) do.call(glmTrainer, x ))
        modelTrained <- res[[which.max(unlist(lapply(res, function(x) x$auc)))]]$model
        param.best <- param[[which.max(unlist(lapply(res, function(x) x$auc)))]]

      }
      if(model =='gbm'){

        gbmTrainer <- function(rsampRate=0.5,csampRate=1, ntrees=50, bal=F,
                               nbins=20, max_depth=4, min_rows=20, learn_rate=0.1){
          modelTrained <- h2o::h2o.gbm(x=2:(ncol(h2oData)-1) , y=ncol(h2oData),
                                  training_frame = h2oData,distribution = "bernoulli",
                                  sample_rate = rsampRate, col_sample_rate=csampRate,
                                  balance_classes = bal, ntrees = ntrees,
                                  max_depth = max_depth, min_rows = min_rows,learn_rate = learn_rate,
                                  nbins=nbins,nfolds = 3)
          param.string <- paste(paste0(names(as.list(match.call()) ),':',as.list(match.call()))[-1], collapse=',')
          writeLines(paste0('GBM model with params: ',param.string,' obtained AUC: ',modelTrained@model$cross_validation_metrics@metrics$AUC))
          auc <- modelTrained@model$cross_validation_metrics@metrics$AUC
          model <- modelTrained
          return(list(auc=auc, model=model))
        }

        # default grid search:
        if(!is.null(param))
          param <- do.call(paramSettings, param)
        if(is.null(param))
          param <- split(expand.grid(bal=c(T,F), rsampRate=c(0.7,0.9,1), ntrees=c(20,50,100)), 1:18)
        
        res <- lapply(param, function(x) do.call(gbmTrainer, x ))
        modelTrained <- res[[which.max(unlist(lapply(res, function(x) x$auc)))]]$model
        param.best <- param[[which.max(unlist(lapply(res, function(x) x$auc)))]]

      }
      writeLines(paste0('Training took: ', format(Sys.time()-start, digits=3)))

      modelLoc <- h2o::h2o.saveModel(modelTrained, path=paste0('file:///',modelLoc), force=T)

      comp <- Sys.time() - start.main
      writeLines(paste0('Training of Model ',model, ' including all formating took:',  format(comp, digits=3)))
      writeLines('Finally saving model details...')

      result <- list(model = NULL,
                     modelLoc = modelLoc,
                     trainAuc = ifelse(is.null(modelTrained@model$cross_validation_metrics@metrics$AUC),
                                       modelTrained@model$training_metrics@metrics$AUC,
                                       modelTrained@model$cross_validation_metrics@metrics$AUC),
                     trainCalibration= NULL,
                     modelSettings = list(model=model,modelParameters=param.best,
                                          featureSettings=featureSettings,
                                          outcomeId=outcomeId, cohortId=cohortId),
                     metaData = plpData$metaData
      )
      class(result) <- 'plpModel'
      saveRDS(result, saveLoc)

    }
  }

  return(result)

}


#' predictPlp
#'
#' @description
#' Predict the risk of the outcome using the input plpModel for the input plpData
#' @details
#' The function applied the trained model on the plpData to make predictions
#' @param plpModel                         An object of type \code{plpModel} - a patient level prediction model
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#'
#' @return
#' An ffdf object containing the prediction for each person in the cohort
#'

#' @export
predictPlp <- function(plpModel, plpData){

  if((!is.null(plpModel$model)) && (plpModel$modelSettings$model=='lr-lasso')){
    prediction <- predictProbabilities(plpModel$model, plpData)
  }

  if( (!is.null(plpModel$model)) && plpModel$modelSettings$model%in%c('nnet','naiveBayes')){
    writeLines('Predicting on validation set...')
    # convert plpData to matrix:
    covariates <- ff::clone(plpData$covariates)
    # if plpModel$features in not null filter non features from covariates
    if(!is.null(plpModel$metaData$usedCovariateIds)){
      writeLines('Extracting covariates used in model...')
      t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(plpModel$metaData$usedCovariateIds$covariateId))
      covariates <- covariates[ffbase::ffwhich(t, !is.na(t)),]
    }

    # now convert into matrix using reshape2::dcast
    writeLines('Converting data into matrix...')
    plpData.matrix <- reshape2::dcast(ff::as.ram(covariates), rowId~covariateId, value.var='covariateValue', fill=0)

    # convert the columns to have X
    colnames(plpData.matrix)[!colnames(plpData.matrix)%in%c('rowId')] <- paste0('X',colnames(plpData.matrix)[!colnames(plpData.matrix)%in%c('rowId')])


    prediction <- predict(plpModel$model$finalModel, plpData.matrix[,-1], type='raw')
    ##prediction <- ffdf(rowId = as.ff(plpData.matrix$rowId), value = as.ff(prediction))
    prediction <- data.frame(rowId = plpData.matrix$rowId, value = prediction[,1])
    ##writeLines(paste(colnames(prediction), sep='', collapse='_'))
    writeLines('Prediction complete...')
    # check whether I need to add rowId name and column.
    prediction <- merge(ff::as.ram(plpData$cohorts), prediction, by='rowId', all.x=T)
    prediction$value[is.na(prediction$value)] <- 0
  }

  # join cohort with prediciton named value
  if(is.null(plpModel$model)  ){

    writeLines('loading model...')
    plpmod <- h2o::h2o.loadModel(plpModel$modelLoc)

    covariates <- ff::clone(plpData$covariates)
    # if plpModel$features in not null filter non features from covariates
    if(!is.null(plpModel$metaData$usedCovariateIds)){
      t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(plpModel$metaData$usedCovariateIds$covariateId))
      covariates <- covariates[ffbase::ffwhich(t, !is.na(t)),]
    }

    # now convert into matrix using reshape2::dcast
    cov <- reshape2::dcast(ff::as.ram(covariates), rowId~covariateId, value.var='covariateValue', fill=0)
    cov.h2o <-h2o::as.h2o(cov)
    value <- h2o::h2o.predict(plpmod, cov.h2o)
    pred <- data.frame(rowId=cov$rowId, value=as.data.frame(value)[,3])
    prediction <- merge(ff::as.ram(plpData$cohorts), pred, by='rowId', all.x=T)
  }
  attr(prediction, "modelType") <- "logistic"
  attr(prediction, "outcomeId") <- plpModel$modelSettings$outcomeId
  return(prediction)
}





#' spliter
#'
#' @description
#' Various train/test splitting techniques
#' @details
#' The function applied the trained model on the plpData to make predictions
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param type                             The type of train/test split:
#'                                         \describe{
#'                                         \item{'time'}{Split the data based on cohort start date -
#'                                         test:pre 2013/train:post2013 }
#'                                         \item{'both'}{Split the data on year 2013 but exclude any train
#'                                         people who are in the test set.}
#'                                         \item{'both'}{Split the data into (1-frac) train set and frac test set.
#'                                         This is stratified by the outcome}
#'                                         }
#' @param frac                              The fraction of people who will go into the test set
#'
#'
#' @return
#' A list with the training set as the first element and the testing set as the second element
#'

#' @export
spliter <- function(plpData, type, frac){
  if(type=='year'){
    split <- list()
    length(split) <- 2

    split[[1]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c('2000-01-01','2012-12-31'),
                               minPriorObservation= NULL #washoutWindow
                               , excludeOutcomeOccurrence=list('1'=c('inf',0)),
                               classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                         insufficientPredictionPeriod = c('include','include'),
                                                         minPostObservation=NULL,
                                                         insufficientPostObservation = c('include','include'),
                                                         survivalCensor=list()
                               ))
    split[[2]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c('2013-01-01','2014-12-31'),
                               minPriorObservation= NULL #washoutWindow
                               , excludeOutcomeOccurrence=list('1'=c('inf',0)),
                               classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                         insufficientPredictionPeriod = c('include','include'),
                                                         minPostObservation=NULL,
                                                         insufficientPostObservation = c('include','include'),
                                                         survivalCensor=list()
                               ))
    return(split)

  }

  if(type=='both'){
    # same as above by remove any rowIds in split[[2]] from split[[1]]
    split <- list()
    length(split) <- 2

    split[[1]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c('2000-01-01','2012-12-31'),
                               minPriorObservation= NULL #washoutWindow
                               , excludeOutcomeOccurrence=list('1'=c('inf',0)),
                               classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                         insufficientPredictionPeriod = c('include','include'),
                                                         minPostObservation=NULL,
                                                         insufficientPostObservation = c('include','include'),
                                                         survivalCensor=list()
                               ))
    split[[2]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c('2013-01-01','2014-12-31'),
                               minPriorObservation= NULL #washoutWindow
                               , excludeOutcomeOccurrence=list('1'=c('inf',0)),
                               classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                         insufficientPredictionPeriod = c('include','include'),
                                                         minPostObservation=NULL,
                                                         insufficientPostObservation = c('include','include'),
                                                         survivalCensor=list()
                               ))

    idExclude <- unique(split[[2]]$covariates$rowId)
    t <- ffbase::ffmatch(split[[1]]$covariates$rowId, table=idExclude)
    split[[1]]$covariates <- split[[1]]$covariates[ffbase::ffwhich(t, is.na(t)),]
    t <- ffbase::ffmatch(split[[1]]$cohorts$rowId, table=idExclude)
    split[[1]]$cohorts <- split[[1]]$cohorts[ffbase::ffwhich(t, is.na(t)),]
    t <- ffbase::ffmatch(split[[1]]$outcomes$rowId, table=idExclude)
    split[[1]]$outcomes <- split[[1]]$outcomes[ffbase::ffwhich(t, is.na(t)),]

    return(split)
  }



  if(type=='patient'){
    # stratified splitting of patients
    covariates <- ff::clone(plpData$covariates)
    cohorts <- ff::clone(plpData$cohorts)
    outcomes <- ff::clone(plpData$outcomes)

    outs <- unique(outcomes$rowId)
    t <-  ffbase::ffmatch(cohorts$rowId, table=outs)
    noOut <- unique(cohorts$rowId[ffbase::ffwhich(t, is.na(t))])

    set.seed =1
    out.train <- sample(ff::as.ram(outs), frac*length(outs))
    set.seed =2
    noout.train <- sample(ff::as.ram(noOut), frac*length(noOut))

    val <- ff::as.ff(c(out.train, noout.train))

    t <- ffbase::ffmatch(covariates$rowId, table=val)
    covariates.val <- covariates[ffbase::ffwhich(t,!is.na(t)),]
    covariates.test <- covariates[ffbase::ffwhich(t,is.na(t)),]
    t <- ffbase::ffmatch(cohorts$rowId, table=val)
    cohorts.val <- cohorts[ffbase::ffwhich(t,!is.na(t)),]
    cohorts.test <- cohorts[ffbase::ffwhich(t,is.na(t)),]
    t <- ffbase::ffmatch(outcomes$rowId, table=val)
    outcomes.val <- outcomes[ffbase::ffwhich(t,!is.na(t)),]
    outcomes.test <- outcomes[ffbase::ffwhich(t,is.na(t)),]

    split <- list()
    length(split) <- 2
    split[[1]] <- list(outcomes = outcomes.test,
                       cohorts = cohorts.test,
                       covariates = covariates.test,
                       covariateRef = ff::clone(plpData$covariateRef),
                       metaData=plpData$metaData
    )

    split[[2]] <- list(outcomes = outcomes.val,
                       cohorts = cohorts.val,
                       covariates = covariates.val,
                       covariateRef = ff::clone(plpData$covariateRef),
                       metaData=plpData$metaData
    )

    return(split)


  }
}



#' makeRandomString
#'
#' @description
#' A function for making a random string
#' @details
#' The function creates n random strings of size length
#' @param n                                An integer - the number of random string to generate
#' @param length                           An integer - the number of characters for each string
#'
#' @return
#' A list containing n random strings with the number of characters specified by the use input length
#'
makeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}


#' paramSettings
#'
#' @description
#' A function for specifying the hyperparameters to create a grid search when training the GBM model
#' @details
#' The function takes a list of the model's hyperparameters and values to investigate while training the model
#' @param rsampRate                        A vector of values between 0 and 1 specifying the fraction of rows to use for each tree
#' @param csampRate                        A vector of values between 0 and 1 specifying the fraction of features to use for each tree
#' @param bal                              A vector of boolean values - specifying whether to balance the class labels during training
#' @param ntrees                           A vector of integers -specifying the number of trees to train
#'
#'
#' @return
#' A list the parameters expanded out like a grid ready to be investigated during model training
#'

#' @export
###paramSettings <- function(rsampRate=1, csampRate=1, bal=F, ntrees=50){
###  return(split(expand.grid(bal=bal, rsampRate=rsampRate, ntrees=ntrees, csampRate=csampRate),
###        1:(length(rsampRate)*length(csampRate)*length(bal)*length(ntrees))))
###}
paramSettings <-function(model='gbm',bal=F, ntrees=100, rsampRate=1, csampRate=1,
                         learn_rate=0.01,
                         nbins=10, max_depth=4, min_rows=20, mtries =-1,
                         alpha=0.5, lambda=1e-06, lambda_search=T,
                         nlambdas=50, lambda_min_ratio=1e-03){

  #sample_rate = rsampRate

  if(model=='gbm'){
    return(split(expand.grid(bal=bal, rsampRate=rsampRate, ntrees=ntrees, csampRate=csampRate,
                             nbins=nbins, max_depth=max_depth, min_rows=min_rows, learn_rate=learn_rate),
                 1:(length(rsampRate)*length(csampRate)*length(bal)*length(ntrees)*length(nbins)*length(max_depth)*length(min_rows)*length(learn_rate)  )))
  }
  if(model=='randomForest'){
    return(split(expand.grid(bal=bal, sample_rate=rsampRate, ntrees=ntrees, mtries=mtries,
                             nbins=nbins, max_depth=max_depth, min_rows=min_rows),
                 1:(length(rsampRate)*length(mtries)*length(bal)*length(ntrees)*length(nbins)*length(max_depth)*length(min_rows)  )))
  }
  if(model=='lr-enet'){
    return(split(expand.grid(alpha=alpha, lambda=lambda, lambda_search=lambda_search,
                             nlambdas=nlambdas, lambda_min_ratio=lambda_min_ratio),
                 1:(length(alpha)*length(lambda)*length(lambda_search)*length(nlambdas)*length(lambda_min_ratio)  )))
  }
}


