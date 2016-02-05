# @file developModelFramework.R
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

#' developModel - Train and evaluate the model
#'
#' @description
#' This provides a general framework for training patient level prediction models.  The user can select 
#' various default feature selection methods or incorporate their own,  The user can also select from
#' a range of default classifiers or incorporate their own.  There are three types of evaluations for the model
#' patient (randomly splits people into train/validation sets) or year (randomly splits data into train/validation sets
#' based on index year - older in training, newer in validation) or both (same as year spliting but checks there are
#' no overlaps in patients within training set and validaiton set - any overlaps are removed from validation set)
#' 
#' @details
#' Users can define a risk period of interest for the prediction of the outcome relative to index or use
#' the cohprt dates.  The user can then specify whether they wish to exclude patients who are not observed
#' during the whole risk period, cohort period or experienced the outcome prior to the risk period.
#'
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param modelSettings                    A list of class \code{modelSettings} containing:
#'                                         \itemize{
#'                                         \item{model -}{ a string specifying the name of classifier function (e.g. 'lr-lasso')}
#'                                         \item{param -}{ a list containing the model parameters, cohortIds, the outcomeIds.}
#'                                         }
#'                                         The default models are:
#'                                         \itemize{
#'                                         \item{lr-lasso -}{ Logistic regression with lasso regularisation - parameters: variance}
#'                                         \item{nnet_plp -}{ Neural network from caret package- parameters: size/decay}
#'                                         \item{svmRadial_plp -}{ SVM with radial kernal from caret package - parameters: C, ...}
#'                                         \item{randomForest_plp -}{ Random forest from h2o package}
#'                                         \item{gbm_plp -}{ Gradient boosting machine from h2o package}
#'                                         \item{lr_enet_plp -}{ Logistic regression with elastic new regularisation from h2o package}
#'                                         }
#' @param featureSettings                  A list of class \code{featureSettings} containing:
#'                                         \itemize{
#'                                         \item{method -}{ a string specifying the name of feature modifying function (e.g. 'wrapperGA')}
#'                                         \item{param -}{ a list containing the method parameters, cohortIds, the outcomeIds.}
#'                                         }
#'
#'                                         The default preprocessing methods are:
#'                                         \itemize{
#'                                         \item{lassolr -}{ Feature selection using lasso logistic regression}
#'                                         \item{wrapperGA -}{ Genetic algorithm wrapper}
#'                                         \item{glrm -}{ (IN PROGRESS)Generaised low rank models}
#'                                         \item{varImp -}{ Variable importance}
#'                                         \item{filterCovariates -}{ Filtering covariates}
#'                                          }
#' @param type                             A subset of c('year','both','patient') specifying the type of evaluation used.
#'                                         'year' find the date where validationFraction of patients had an index after the date and assigns patients with an index prior to this date into the training set and post the date into the test set
#'                                         'both' splits the data by the year but removes any patient in the test set from the training set
#'                                         'patient' splits the data into test (1-validationFraction of the data) and
#'                                         train (validationFraction of the data) sets.  The split is stratified by the class label.
#' @param validationFraction               The fraction of the data to be used as the validation set in the patient
#'                                         split evaluation.
#' @param fileLoc                          The path to the directory where the models will be saved

#'
#' @return
#' An object containing the model or location where the model is save, the data selection settings, the preprocessing
#' and training settings as well as various performance measures obtained by the model.
#'
#' \item{model}{A list of class \code{plpModel} containing the model, training metrics and model metadata}
#' \item{dataSummary}{A list detailing the size of the train/test sets and outcome prevalence}
#' \item{evalType}{The type of evaluation that was performed}
#' \item{prediction}{An ffdf object containing the prediction for each person in the validation set }
#' \item{performance}{A list detailing the performance of the model}
#' \item{time}{The complete time taken to do the model framework}
#'
#'
#' @export
#' @examples
#' #******** EXAMPLE 1 ********* 
#' #lasso logistic regression oredicting outcome 2 in cohorts 1 and 3 
#' #using no feature selection with a year split evaluation:
#' modset_llr  <- list(model='lr_lasso',
#'                    param=list(variance =0.001, cohortId=c(1,2), outcomeId=2))
#' class(modset_llr) <- 'modelSettings'
#' mod_llr <- developModel2(plpData= plpData,
#'                         featureSettings = NULL,
#'                         modelSettings = modset_llr ,
#'                         type='year')
#'  
#' #******** EXAMPLE 2 *********                                               
#' # Gradient boosting machine using a genetic algorimth wrapper to 
#' # select the feature subset and a grid search to select hyper parameters                         
#' featSet_gbm <- list(method='wrapperGA', param=list(cohortId=c(1,2), outcomeId=2, varSize=300, iter=25))
#' class(featSet_gbm) <- 'featureSettings'
#' modset_gbm <- list(model='gbm_plp',
#'                    param=list(rsampRate=0.8, ntrees=c(100,150), max_depth=c(2,4,5), cohortId=c(1,2), outcomeId=2))
#' class(modset_gbm) <- 'modelSettings'
#' mod_gbm <- developModel2(plpData= plpData.censor,
#'                          featureSettings = featSet_gbm,
#'                          modelSettings = modset_gbm,
#'                          type='year')
#' 
developModel2 <- function(plpData,
                         modelSettings,
                         featureSettings,
                         type = c('year','both','patient'), validationFraction=0.2, 
                         fileLoc=file.path(getwd(),'models')
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
    start.all <- Sys.time()
    data <- spliter(plpData, type[i], validationFraction)
    settings[[i]] <- list(data=data[[1]],
                          modelSettings = modelSettings,
                          featureSettings = featureSettings,
                          cohortId = modelSettings$cohortId,
                          outcomeId = modelSettings$outcomeId,
                          loc=fileLoc)
    
    dataSummary= list(trainCohort =nrow(data[[1]]$cohorts),
                      trainOutcomeCount =nrow(data[[1]]$outcomes),
                      testCohort =nrow(data[[2]]$cohorts),
                      testOutcomeCount =nrow(data[[2]]$outcomes),
                      covariateCount = length(unique(data[[1]]$covariates$covariateId))
    )
    
    # train model
    model <- do.call(fitPlp2, settings[[i]])
    writeLines('1: Model Trained')
    
    prediction <- NULL
    performance <- NULL
    if(model!='no cov'){
      # do prediction
      prediction <- predictPlp2(plpModel=model, plpData=data[[2]])
      writeLines('2) Prediction Calculated')
      # calculate metrics
      performance <- evaluatePlp(prediction, data[[2]])
      writeLines('3) Performance calculated')
    }
    comp <- Sys.time() - start.all
    results[[i]] <- list(model=model,
                         dataSummary=dataSummary,
                         evalType=type[i],
                         prediction=prediction,
                         performance=performance,
                         time=comp)
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
#' 
#' @param data                             An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param modelSettings                    A list of class \code{modelSettings} containing:
#'                                         \itemize{
#'                                         \item{model -}{ a string specifying the name of classifier function (e.g. 'lr-lasso')}
#'                                         \item{param -}{ a list containing the model parameters, cohortIds, the outcomeIds.}
#'                                         }
#'                                         The default models are:
#'                                         \itemize{
#'                                         \item{lr-lasso -}{ Logistic regression with lasso regularisation - parameters: variance}
#'                                         \item{nnet_plp -}{ Neural network from caret package- parameters: size/decay}
#'                                         \item{svmRadial_plp -}{ SVM with radial kernal from caret package - parameters: C, ...}
#'                                         \item{randomForest_plp -}{ Random forest from h2o package}
#'                                         \item{gbm_plp -}{ Gradient boosting machine from h2o package}
#'                                         \item{lr_enet_plp -}{ Logistic regression with elastic new regularisation from h2o package}
#'                                         }
#' @param featureSettings                  A list of class \code{featureSettings} containing:
#'                                         \itemize{
#'                                         \item{method -}{ a string specifying the name of feature modifying function (e.g. 'wrapperGA')}
#'                                         \item{param -}{ a list containing the method parameters, cohortIds, the outcomeIds.}
#'                                         }
#'
#'                                         The default preprocessing methods are:
#'                                         \itemize{
#'                                         \item{lassolr -}{ Feature selection using lasso logistic regression}
#'                                         \item{wrapperGA -}{ Genetic algorithm wrapper}
#'                                         \item{glrm -}{ (IN PROGRESS)Generaised low rank models}
#'                                         \item{varImp -}{ Variable importance}
#'                                         \item{filterCovariates -}{ Filtering covariates}
#'                                          }
#' @param outcomeId                        The outcomeId the user is aiming to predict
#' @param cohortId                         The id of the cohort being used.  Default is NULL.
#' @param loc                              The path to the directory where the model will be saved

#'
#' @return
#' An object of class \code{plpModel} containing:
#' 
#' \item{model}{The trained prediction model}
#' \item{modelLoc}{The path to where the model is saved (if saved)}
#' \item{trainAuc}{The AUC obtained on the training set}
#' \item{trainCalibration}{The calibration obtained on the training set}
#' \item{modelSettings}{A list specifiying the model, preprocessing, outcomeId and cohortId}
#' \item{metaData}{The model meta data}
#' \item{trainingTime}{The time taken to train the classifier}
#'
#'

#' @export
fitPlp2 <- function(data, modelSettings, featureSettings, outcomeId, cohortId, loc){
  plpData <- list(outcomes =ff::clone(data$outcomes),
                  cohorts =ff::clone(data$cohorts),
                  covariates =ff::clone(data$covariates),
                  exclude =ff::clone(data$exclude),
                  covariateRef=ff::clone(data$covariateRef),
                  #metaData=plpData$metaData
                  metaData=data$metaData
  )
  
  # set saving locations
  #saveLoc <- file.path(loc, paste0(model,'_',makeRandomString(),'.rds'))

  # apply each featureSetting option in order of the list entry to do feature engineering/selection
  if (class(featureSettings) == "featureSettings") {
    #fun <- attr(featureSettings, "fun")
    fun <- featureSettings$method
    args <- c(list(plpData =plpData),featureSettings$param)
    plpData <- do.call(fun, args)
    
    if (nrow(plpData$covariates) == 0) {
      warning("No features remaining")
    } 
  } else if (is.list(featureSettings)) {
    for (i in 1:length(featureSettings)) {
      #fun <- attr(featureSettings[[i]], "fun")
      fun <- featureSettings[[i]]$method
      args <- c(list(plpData=plpData),featureSettings[[i]]$param)
      plpData <- do.call(fun, args)
      
      if (nrow(plpData$covariates) == 0) {
        warning("No features remaining")
      } 
    }
  }
  
  # Now apply the classifier:
    #fun <- attr(modelSettings, "fun")
    fun <- modelSettings$model
    args <- list(plpData =plpData,param =modelSettings$param)
    plpModel <- do.call(fun, args)
    
    return(plpModel)

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
predictPlp2 <- function(plpModel, plpData){
  # in the model part add an attribute type - plp, h2o, caret, ... then apply prediciton for that type or allow custom
  
  # apply the feature transformations
  if(class(plpModel$metaData$featureSettings)=='list' & !'transform'%in%names(plpModel$metaData$featureSettings)){
    for(x in plpModel$metaData$featureSettings){
      plpData <- do.call(x$transform, list(plpData2=plpData))
    }
  }
  if('transform'%in%names(plpModel$metaData$featureSettings) ){
    plpData <- do.call(plpModel$metaData$featureSettings$transform, list(plpData2=plpData))
  }
    
  # do the predction on the new data
  if(class(plpModel)=='plpModel'){
    # extract the classifier type
    pred <- paste0('predict.',attr(plpModel, 'type'))
  prediction <-  do.call(pred, list(plpModel=plpModel, plpData=plpData))
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
#' The function splits the plpData into test/train sets
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param type                             The type of train/test split:
#'                                         \itemize{
#'                                         \item{'year'}{Split the data based on cohort start date -
#'                                         test:pre 2013/train:post2013 }
#'                                         \item{'both'}{Split the data on year 2013 but exclude any train
#'                                         people who are in the test set.}
#'                                         \item{'patient'}{Split the data into (1-frac) train set and frac test set.
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
    writeLines(paste0("Splitting data by date based on approximately ",frac*100,"% being in training set"))
    split <- list()
    length(split) <- 2
    
    cohorts <- ff::clone(plpData$cohorts)
    cohortDates <- as.Date(ff::as.ram(cohorts$cohortStartDate), format = "%Y-%m-%d")
    cohortDates <- cohortDates[order(cohortDates)]
    ind <- floor(length(cohortDates)*(1-frac))
    split.d1 <- cohortDates[ind]
    
    writeLines(paste0('Splitting on date: ', split.d1))
    t <- as.Date(ff::as.ram(cohorts$cohortStartDate), format = "%Y-%m-%d") <= split.d1     
    t <- ff::as.ff(t)
    train.cohorts <- cohorts[ffbase::ffwhich(t, t==T),]
    t <- as.Date(ff::as.ram(cohorts$cohortStartDate), format = "%Y-%m-%d") > split.d1     
    t <- ff::as.ff(t)
    test.cohorts <- cohorts[ffbase::ffwhich(t, t==T),]
    writeLines(paste0('Train size: ', nrow(train.cohorts), ' -- Test size:',nrow(test.cohorts)))
    
    covariates <- ff::clone(plpData$covariates)
    t <- ffbase::ffmatch(covariates$rowId, table=ff::as.ff(unique(test.cohorts$rowId)))
    covariates.test <- covariates[ffbase::ffwhich(t, !is.na(t)),]
    t <- ffbase::ffmatch(covariates$rowId, table=ff::as.ff(unique(train.cohorts$rowId)))
    covariates.train <- covariates[ffbase::ffwhich(t, !is.na(t)),]
    
    outcomes <- ff::clone(plpData$outcomes )
    t <- ffbase::ffmatch(outcomes$rowId, table=ff::as.ff(unique(test.cohorts$rowId)))
    outcomes.test <- outcomes[ffbase::ffwhich(t, !is.na(t)),]
    t <- ffbase::ffmatch(outcomes$rowId, table=ff::as.ff(unique(train.cohorts$rowId)))
    outcomes.train <- outcomes[ffbase::ffwhich(t, !is.na(t)),]
    
    writeLines(paste0(ncol(outcomes.test), '-', ncol(covariates.test)))
    writeLines(paste0(ncol(outcomes.train), '-', ncol(covariates.train)))
    
    split[[1]] <- list(cohorts=train.cohorts,
                       covariates=covariates.train,
                       outcomes=outcomes.train,
                       covariateRef = ff::clone(plpData$covariateRef),
                       metaData = plpData$metaData)
    split[[2]] <- list(cohorts=test.cohorts,
                       covariates=covariates.test,
                       outcomes=outcomes.test,
                       covariateRef = ff::clone(plpData$covariateRef),
                       metaData = plpData$metaData)
    
    return(split)
    
  }
  
  if(type=='both'){
    # same as above by remove any rowIds in split[[2]] from split[[1]]
    split <- list()
    length(split) <- 2
    
    cohortDates <- as.ram(plpData$cohorts$cohortStartDate)
    cohortDates <- cohortDates[order(cohortDates)]
    ind <- floor(length(cohortDates)*(1-frac))
    split.d1 <- cohortDates[ind]
    split.d2 <- cohortDates[ind]+1
    
    split[[1]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c(min(cohortdates),split.d1),
                               minPriorObservation= NULL #washoutWindow
                               , excludeOutcomeOccurrence=list('1'=c('inf',0)),
                               classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                         insufficientPredictionPeriod = c('include','include'),
                                                         minPostObservation=NULL,
                                                         insufficientPostObservation = c('include','include'),
                                                         survivalCensor=list()
                               ))
    split[[2]]<- censorPlpData(plpData, predictionPeriod =NULL,  dateInterval=c(split.d2,max(cohortDates)),
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
                         nlambdas=50, lambda_min_ratio=1e-03,...){
  
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


