.testcode <- function(){
  library(ff)
  library(ffbase)
  library(reshape2)
  library(h2o)
  library(caret)
  library(PatientLevelPrediction)
  options(fftempdir = "s:/FFtemp")
  #h2o.init(nthreads = -1, max_mem_size = '50G')

  setwd('C:/Users/jreps/Documents')

  cdmDatabaseSchemaList <- c("cdm_jmdc_v5.dbo","cdm_optum_v5.dbo","cdm_cprd_v5.dbo","cdm_truven_mdcd_v5.dbo","cdm_truven_ccae_v5.dbo")
  cdmDatabaseSchema <- cdmDatabaseSchemaList[1]
  database <- strsplit(cdmDatabaseSchema, '\\.')[[1]][1]

  plpData <- loadPlpData(file.path("s:/gestDia/data_new", database))
  dim(plpData$cohorts)
  dim(plpData$outcomes)
  dim(plpData$covariates)
  dim(plpData$exclude)

  classificationCensor=list(insufficientCohortObservation = c('exclude','exclude'),
                            insufficientPredictionPeriod = c('include','include'),
                            minPostObservation=365,
                            insufficientPostObservation = c('include','exclude')
  )

  plpData.censored <- censorPlpData(plpData, outcomeTime=NULL, newOutcome=NULL,
                                    minCohortTime=150,
                                    predictionPeriod =c(147,200),  dateInterval=c('2008-01-01','2013-12-31'),
                                    minPriorObservation= 2*365 #washoutWindow
                                    ,excludeOutcomeOccurrence=list('3'=c(2*365,30),'1'=c('inf',147),
                                                                   '4'=c(2*365,30)),
                                    classificationCensor,
                                    survivalCensor=list())

  #unique(plpData.censored$exclude$outcomeId)

  dim(plpData.censored$cohorts)
  dim(plpData.censored$outcomes)
  t <- plpData.censored$outcomes$outcomeId==2
  dim(plpData.censored$outcomes[ffwhich(t, t==T),])
  dim(plpData.censored$covariates)
  dim(plpData.censored$exclude)

  # get a description of the data:
  dataDesc <- describePlpData(plpData.censored, covariateVals=c('riskFactors','age'), cdmDatabase=database,
                              agehist=TRUE,
                              plot=T, plotFile=NULL, saveTable=F, tableFile=NULL,
                              outcomeName='Pregnancy with GDM',
                              cohortName='Total pregnancy',
                              perYear=T)

  # get lasso lr benchmark:
  ##plpData.censored$outcomes$outcomeId <- ff(rep(2, length(plpData.censored$outcomes$outcomeId)))
  lr <- developModel(plpData.censored,modelSettings=list(model='lr-lasso', cohortId=NULL, outcomeId=2,
                                                         param=list(val=0.5)),
                     validationFraction=0.1, fileLoc=file.path(getwd(),'test','models'),
                     type='both'
  )
  lr[[1]]$performance$auc
  var.imp.lr <- data.frame(covariateId = names(lr[[1]]$model$model$coefficients[lr[[1]]$model$model$coefficients!=0]),
                           imp =lr[[1]]$model$model$coefficients[lr[[1]]$model$model$coefficients!=0])
  var.imp.lr <- merge(as.ram(plpData.censored$covariateRef), var.imp.lr, by='covariateId')
  var.imp.lr <- var.imp.lr[order(-abs(var.imp.lr$imp)),]

  trainModel <- developModel(plpData.censored,modelSettings=list(model='nnet', cohortId=NULL, outcomeId=2,
                                                                 #preprocess=c('allEra','lr-lasso')  ),
                                                                 #param=list(model='gbm',rsampRate=c(0.8), ntrees=c(100),
                                                                #            learn_rate=c(0.01,0.1,0.05))),
                                                                #param=list(model='randomForest',rsampRate=c(0.8), ntrees=c(100),
                                                                #           mtries=c(-1, 10,50))),
                                                                #param=list(model='lr-enet',lambda=c(0.000001), alpha=c(0,0.2,0.5,0.8,1)
                                                                param=list(model='nnet',size=c(2), decay=c(0.2)

                                                                )),
                             featureSettings = list(analysisSelector=c(-(1:16),4,201,505),
                                                    covariateSelector=NULL,
                                                    wrapper=list(method='lr-lasso', variance=0.1),
                                                    matrixFactor=F,outcomeId=2),
                             validationFraction=0.1, fileLoc=file.path(getwd(),'test','models'),
                             type='both'
                             )

  # lasso benchmark: 0.7154
  trainModel[[1]]$performance$auc
  mod <- h2o.loadModel(trainModel[[1]]$model$modelLoc)
  mod@allparameters
  as.data.frame(mod@model$variable_importances)[1:100,]

  varImp <- as.data.frame(mod@model$variable_importances)[1:100,]
  colnames(varImp)[1] <- 'covariateId'
  varImp <- merge(varImp, plpData.censored$covariateRef, by='covariateId')
  varImp.gbm.optum <- varImp[order(-varImp$relative_importance),]
  trainModel[[2]]$performance$auc
  trainModel[[3]]$performance$auc
  #
}

