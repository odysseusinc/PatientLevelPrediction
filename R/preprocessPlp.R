# @file featureReduction.R
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

#' featureReduction
#'
#' @description
#' Applies a variety of feature selection/engineering techniques
#' @details
#' Users can pick from selecting a type of feature, applying wrapper feature selection using logistic regression
#' with the lasso regularisation or using matrix factorisation.
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

#'
#' @return
#' A list of class plpData containing the original plpData but with a reduced number of covariates and
#' details of the feature selection/engineering stored in the metaData
#' \describe{ \item{cohorts}{An ffdf object listing all persons and their prediction periods. This
#' object will have these fields: row_id (a unique ID per period), person_id, cohort_start_date,
#' cohort_id, time (number of days in the window).} \item{outcomes}{An ffdf object listing all
#' outcomes per period. This object will have these fields: row_id, outcome_id, outcome_count,
#' time_to_event.} \item{exclude}{Either NULL or an ffdf object listing per outcome ID which windows
#' had the outcome prior to the window. This object will have these fields: rowId, outcomeId.}
#' \item{covariates}{An ffdf object listing the baseline covariates per person in the cohorts. This is
#' done using a sparse representation: covariates with a value of 0 are omitted to save space. The
#' covariates object will have three columns: rowId, covariateId, and covariateValue. }
#' \item{covariateRef}{An ffdf object describing the covariates that have been extracted.}
#' \item{metaData}{A list of objects with information on how the plpData object was constructed
#' and feature reduction details.  The list member named 'featureReducer' contains the parameters
#' input into the featureReduction and the member named 'usedCovariateIds' contains the set of covariteIds included.} }
#'
#' @export


featureReducer <- function(plpData, analysisSelector=c(-(1:16),4,201,505),
                           covariateSelector=NULL, wrapper=NULL,
                           matrixFactor=F, outcomeId=1, cohortId=NULL){

  riskfact <- NULL

  # feature selection by simply picking the types of covariates
  if(!is.null(analysisSelector) | !is.null(covariateSelector)){
    if(!is.null(analysisSelector)){
      writeLines(paste0('Extracting all covariates with analysisId in:', paste(analysisSelector, collapse=', ')))
      covsel <- plpData$covariateRef$covariateId[ff::as.ram(plpData$covariateRef$analysisId)%in%analysisSelector]
    }
  if(!is.null(covariateSelector)){
    writeLines(paste0('Extracting all covariates with covariateId in:', paste(covariateSelector, collapse=', ')))

    if(!is.null(analysisSelector))
      covsel <- c(covsel, ff::as.ff(covariateSelector))
    if(is.null(analysisSelector))
      covsel <- ff::as.ff(covariateSelector)
  }

    t <- ffbase::ffmatch(plpData$covariates$covariateId, table=ff::as.ff(covsel))
    plpData$covariates <-  plpData$covariates[ffbase::ffwhich(t, !is.na(t)),]

    writeLines(paste0('Found ',length(unique(plpData$covariates$covariateId)),' with selected analysisIds/covariateIds'))
    riskfact <- data.frame(covariateId=ff::as.ram(unique(plpData$covariates$covariateId)))
  }


  # feature selection using matric factorisation
  if(matrixFactor==T){
    # convert to svm light
    rowIds <- unique(plpData$cohorts$rowId)
    # cast by 10,000 patients at a time - TO SLOW for 30k features!
    for (i in 1:ceiling(length(rowIds)/10000)){
      t <- ffbase::ffmatch(plpData$covariates$rowId, table=ff::as.ff(rowIds[( (i-1)*10000+1):min(10000*i, length(rowIds))]))
      cov <- reshape2::dcast(ff::as.ram(plpData$covariates[ffbase::ffwhich(t,!is.na(t)),]),
                   rowId~covariateId, value.var='covariateValue', fill=0)
      allData <- merge(cov, ff::as.ram(plpData$outcomes[,c('rowId','outcomeCount')]), by='rowId', all.x=T)
      allData$outcomeCount[is.na(allData$outcomeCount)] <- 0
      allData$outcomeCount <- as.factor(allData$outcomeCount)
      if(i==1)
        allData.h2o <- h2o::as.h2o(allData)
      if(i>1)
        allData.h2o <- h2o::h2o.rbind(allData.h2o, h2o::as.h2o(allData) )
    }
  }


  # perfprming lasso lr wrapper feature selection - this is forced when number of features > 500
  if(wrapper$method=='lr-lasso' | (is.null(wrapper$method)&length(unique(plpData$covariates$covariateId))>500) ){
    writeLines(paste0('Performing feature selection using lasso logistic regression - this may be due to more than 500 covariates'))
    start <- Sys.time()
    val <- 0.002
    if(!is.null(wrapper$variance)) val <- wrapper$variance
    feature <- fitPredictiveModel(plpData = plpData,
                                  modelType = "logistic",
                                  removeDropoutsForLr = F,
                                  cohortId = cohortId,
                                  outcomeId = outcomeId,
                                  prior = createPrior("laplace",
                                                      exclude = c(0),
                                                      variance = val))

    riskfact <- data.frame(covariateId=names(feature$coefficients)[feature$coefficients!=0],
                           value=feature$coefficients[feature$coefficients!=0])

    t <- ffbase::ffmatch(plpData$covariates$covariateId, table=ff::as.ff(riskfact$covariateId))
    if(sum(!is.na(t))>0)
      plpData$covariates <- plpData$covariates[ffbase::ffwhich(t, !is.na(t)),]
    if(sum(!is.na(t))==0)
      plpData$covariates <- NULL

    writeLines(paste0('Feature select using lasso logistic regression took: ', format(Sys.time()-start, digits=3)))
    writeLines(paste0('Found ', nrow(riskfact), ' risk factors'))

  }

  # save the details into the plpData
  plpData$metaData$featureReducer <- match.call()
  if(!is.null(riskfact))
    plpData$metaData$usedCovariateIds <- riskfact

  return(plpData)

}


