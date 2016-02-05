filterCovariates <- function(plpData, include, conceptIds=NULL, 
                             covariateIds=NULL, analysisIds=NULL, quiet=F,...){
  covariates <- clone(plpData$covariates)
  allCovs <- c()
  #covRef: covariateId  covariateName analysisId conceptId
  
  # include/exclude covariates as specified
  if(!is.null(conceptIds)){
    if(!quiet)
      writeLines('Extracting concepts used in model...')
    t <- ffbase::ffmatch(plpData$covariateRef$conceptId, table=ff::as.ff(conceptIds))
    covarIds <- unique(plpData$covariateRef$covariateId[ffbase::ffwhich(t, !is.na(t)),])
    allCovs <- c(allCovs, ff::as.ram(covarIds)) 
    t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(covarIds))
    if(include)
      covariates<- covariates[ffbase::ffwhich(t, !is.na(t)),]
    if(!include)
      covariates<- covariates[ffbase::ffwhich(t, is.na(t)),]
  }
  
  if(!is.null(covariateIds)){
    if(!quiet)
      writeLines('Extracting covarites used in model...')
    allCovs <- c(allCovs, covariateIds) 
    t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(covariateIds))
    if(include)
      covariates<- covariates[ffbase::ffwhich(t, !is.na(t)),]
    if(!include)
      covariates<- covariates[ffbase::ffwhich(t, is.na(t)),]
  }
  
  if(!is.null(analysisIds)){
    if(!quiet)
      writeLines('Extracting analysis concepts used in model...')
    t <- ffbase::ffmatch(plpData$covariateRef$analysisId, table=ff::as.ff(analysisIds))
    if(sum(!is.na(t))==0){return(warning('No covariates'))}
    covarIds <- unique(plpData$covariateRef$covariateId[ffbase::ffwhich(t, !is.na(t)),])
    allCovs <- c(allCovs, ff::as.ram(covarIds)) 
    t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(covarIds))
    if(include)
      covariates<- covariates[ffbase::ffwhich(t, !is.na(t)),]
    if(!include)
      covariates<- covariates[ffbase::ffwhich(t, is.na(t)),]
  }
  
  # create a transform function to replicate this for new data:
  transformData <- function(plpData2){
    writeLines('transform function running...')
    newcovs <- clone(plpData2$covariates)
    t <- ffbase::ffmatch(newcovs$covariateId, table=ff::as.ff(allCovs))
    if(include)
      newcovs<- newcovs[ffbase::ffwhich(t, !is.na(t)),]
    if(!include)
      newcovs<- newcovs[ffbase::ffwhich(t, is.na(t)),]
    
    metaData <- plpData2$metaData
    if(!is.null(metaData$featureInclude) & include)
      metaData$featureInclude <- c(metaData$featureInclude,allCovs)
    if(!is.null(metaData$featureExclude) & !include)
      metaData$featureExclude <- c(metaData$featureExclude,allCovs)
    if(is.null(metaData$featureInclude) & include)
      metaData$featureInclude <- allCovs
    if(is.null(metaData$featureExclude) & !include)
      metaData$featureExclude <- allCovs

    newData <- list(covariates = newcovs,
                    cohorts = ff::clone(plpData2$cohorts),
                    outcomes = ff::clone(plpData2$outcomes),
                    covariateRef = ff::clone(plpData2$covariateRef),
                    metaData = metaData
                    )
   return(newData) 
  }
  
# update metaData with transform functiuon and included/excluded covs:
  metaDataUpdate <- plpData$metaData
  if(!is.null(metaDataUpdate$featureInclude) & include)
    metaDataUpdate$featureInclude <- c(metaDataUpdate$featureInclude,allCovs)
  if(!is.null(metaDataUpdate$featureExclude) & !include)
    metaDataUpdate$featureExclude <- c(metaDataUpdate$featureExclude,allCovs)
  if(is.null(metaDataUpdate$featureInclude) & include)
    metaDataUpdate$featureInclude <- allCovs
  if(is.null(metaDataUpdate$featureExclude) & !include)
    metaDataUpdate$featureExclude <- allCovs
    
  featureSet <- list(method='filterCovariates', 
                        covariateIds=covariateIds,
                        conceptIds = conceptIds,
                        analysisIds = analysisIds,
                        include=include,
                        all = allCovs,
                        transform=transformData)
  if(!is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- list(metaDataUpdate$featureSettings, featureSet )
  if(is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- featureSet 
  
  result <- list(covariates = covariates,
                  cohorts = ff::clone(plpData$cohorts),
                  outcomes = ff::clone(plpData$outcomes),
                  covariateRef = ff::clone(plpData$covariateRef),
                  metaData = metaDataUpdate
                 )
  return(result)
}

lassolr <- function(plpData, cohortId, outcomeId, variance=NULL,quiet=T, ...){
  # apply lasso lr to select subset of features:
  if(is.null(variance)){
    variance <- 0.01
  }
  modelTrained <- fitPredictiveModel(plpData = plpData,
                                     modelType = "logistic",
                                     removeDropoutsForLr = F,
                                     cohortId = cohortId,
                                     outcomeId = outcomeId,
                                     prior = createPrior("laplace",
                                                         exclude = c(0),
                                                         variance = variance))
  
  # find the selected covariates:
  covs <- clone(plpData$covariates)
  lassocovs <- as.double(names(modelTrained$coefficients[-1])[!modelTrained$coefficients[-1]==0])
  
  if(!quiet)
    writeLines('Extracting covariates found by lasso logistic regression...')
  t <- ffbase::ffmatch(covs$covariateId, table=ff::as.ff(lassocovs))
  covs <- covs[ffbase::ffwhich(t, !is.na(t)),]
  
  # create a transform function to replicate this for new data:
  transformData <- function(plpData2){
    writeLines('transform function running...')
    newcovs <- clone(plpData2$covariates)
    t <- ffbase::ffmatch(newcovs$covariateId, table=ff::as.ff(lassocovs))

    newcovs<- newcovs[ffbase::ffwhich(t, !is.na(t)),]

    metaData <- plpData2$metaData
    if(!is.null(metaData$featureInclude) )
      metaData$featureInclude <- c(metaData$featureInclude,lassocovs)
    if(is.null(metaData$featureInclude) )
      metaData$featureInclude <- lassocovs

    
    newData <- list(covariates = newcovs,
                    cohorts = ff::clone(plpData2$cohorts),
                    outcomes = ff::clone(plpData2$outcomes),
                    covariateRef = ff::clone(plpData2$covariateRef),
                    metaData = metaData
    )
    return(newData) 
  }
  
  # add the covariatesIncluded metadata list of featureSellect 
  metaDataUpdate <- plpData$metaData
  if(!is.null(metaDataUpdate$featureInclude))
    metaDataUpdate$featureInclude <- c(metaDataUpdate$featureInclude,lassocovs)
  if(is.null(metaDataUpdate$featureInclude))
    metaDataUpdate$featureInclude <- lassocovs
  
  featureSet <- list(method='lassoLR_select', 
                     variance = variance,
                     transform=transformData)
  if(!is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- list(metaDataUpdate$featureSettings, featureSet )
  if(is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- featureSet 
  
  result <- list(covariates = covs,
                 cohorts = ff::clone(plpData$cohorts),
                 outcomes = ff::clone(plpData$outcomes),
                 covariateRef = ff::clone(plpData$covariateRef),
                 metaData = metaDataUpdate
  )
  return(result)
}

glrm <- function(plpdata, clusterSize, ...){
  #convert to matrix
  
  #find topics
  
  #figure how to predict topics for new data and save this
  
  # save into covariate
  
  
}

wrapperGA <- function(plpData, model='gbm', varSize=300, iter=100, ...){
  #take small sample and clone data
  cohorts <-ff::clone(plpData$cohorts)
  outcomes <-ff::clone(plpData$outcomes)
  ppl.gen <- unique(as.ram(cohorts$rowId))
  ppl.out <- unique(as.ram(outcomes$rowId))
  if(length(ppl.gen>1500)){
    ppl.gen <- sample(ppl.gen, 1500)
    if(length(ppl.out)>1000)
      ppl.out <- sample(ppl.out, 1000)
  }
  
  covariates <- ff::clone(plpData$covariates)
  t <- ffbase::ffmatch(covariates$rowId, table=ff::as.ff(c(ppl.out, ppl.gen)))
  covariates.samp<- covariates[ffbase::ffwhich(t, !is.na(t)),]
  
  covMat <- reshape2::dcast(ff::as.ram(covariates.samp), rowId~covariateId, 
                            value.var='covariateValue', fill=0)
  covMat <- merge(covMat, as.ram(outcomes[,c('rowId','outcomeCount')]), by='rowId', all.x=T)
  writeLines(paste(covMat$outcomeCount[1:20], sep=',',collapse=','))
  labs <- rep(1, nrow(covMat))
  labs[is.na(covMat$outcomeCount)] <- 0
  covMat$outcomeCount <- labs
  writeLines(paste(covMat$outcomeCount[1:20], sep=',',collapse=','))
  #writeLines(paste0('dim: ', dim(covMat)))
  covMat <- h2o::as.h2o(covMat[,-1]) # remove rowId
  #writeLines(paste(as.data.frame(covMat$outcomeCount[1:20]), sep=',',collapse=','))
  
  covMat$outcomeCount <- h2o::as.factor(covMat$outcomeCount)
  # do genetic algorithm to search for subset of size varSize
  # eval
  evalFun <- function(x1){
    data <- covMat[,c(x1==1,T)]
    mod <- h2o::h2o.gbm(x=1:sum(x1==1), y=sum(x1==1)+1, training_frame=data, max_depth=4, ntrees=50, nfolds=5)
    return(1-mod@model$cross_validation_metrics@metrics$AUC)
  }
    
  require(genalg)
  ga.res <- genalg::rbga.bin(size=ncol(covMat)-1, 
                             popSize=20, iters=iter,
                             zeroToOneRatio=ncol(covMat)/varSize, mutationChance=0.01,
                             evalFunc = evalFun,
                             elitism=4, verbose = T)
  
  plot(1-ga.res$best)
  minEval = min(ga.res$evaluations);
  writeLines(paste0('best auc:',1-minEval))
  filter = ga.res$evaluations == minEval;
  ind <- ga.res$population[filter,]
  if(sum(filter)>1)
    ind <- ga.res$population[filter,][1,]
  selFeat <- colnames(covMat[,-ncol(covMat)])[ind==1]

  t <- ffbase::ffmatch(covariates$covariateId, table=ff::as.ff(as.double(selFeat)))
  covariates <- covariates[ffbase::ffwhich(t,!is.na(t)),]
  
  # create the transformData function
  transformData <- function(plpData2){
    writeLines('transform function running...')
    newcovs <- ff::clone(plpData2$covariates)
    t <- ffbase::ffmatch(newcovs$covariateId, table=ff::as.ff(as.double(selFeat)))
    
    newcovs<- newcovs[ffbase::ffwhich(t, !is.na(t)),]
    
    metaData <- plpData2$metaData
    if(!is.null(metaData$featureInclude) )
      metaData$featureInclude <- c(metaData$featureInclude,as.double(selFeat))
    if(is.null(metaData$featureInclude) )
      metaData$featureInclude <- as.double(selFeat)
    
    
    newData <- list(covariates = newcovs,
                    cohorts = ff::clone(plpData2$cohorts),
                    outcomes = ff::clone(plpData2$outcomes),
                    covariateRef = ff::clone(plpData2$covariateRef),
                    metaData = metaData
    )
    return(newData) 
  }
  
  
  # update metaData
  # add the covariatesIncluded metadata list of featureSellect 
  metaDataUpdate <- plpData$metaData
  if(!is.null(metaDataUpdate$featureInclude))
    metaDataUpdate$featureInclude <- c(metaDataUpdate$featureInclude,as.double(selFeat))
  if(is.null(metaDataUpdate$featureInclude))
    metaDataUpdate$featureInclude <- as.double(selFeat)
  
  featureSet <- list(method='wrapperGA', 
                     varSize = varSize,
                     iters=iter,
                     transform=transformData)
  if(!is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- list(metaDataUpdate$featureSettings, featureSet )
  if(is.null(metaDataUpdate$featureSettings))
    metaDataUpdate$featureSettings <- featureSet 
  
  result <- list(covariates = covariates,
                 cohorts = cohorts,
                 outcomes = outcomes,
                 covariateRef = ff::clone(plpData$covariateRef),
                 metaData = metaDataUpdate
  )
  return(result)
  
  
}

