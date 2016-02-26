# @file censorPlpData.R
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

#' Filters the plpData based on user specified criteria
#'
#' @description
#' Filters the data based on user censoring specifications for classification or survival analysis
#'
#' @details
#' Users can define a risk period of interest for the prediction of the outcome relative to index or use 
#' the cohprt dates.  The user can then specify whether they wish to exclude patients who are not observed
#' during the whole risk period, cohort period or experienced the outcome prior to the risk period.
#'
#' @param plpData                          An object of type \code{plpData} - the patient level prediction 
#'                                         data extracted from the CDM.
#' @param outcomeIds                       a vector of integers (corresponding to outcome ids) or NULL 
#' @param outcomeTime                      An integer - if you need to edit the time from index where you want to predict the
#'                                         outcome use this parameter.  For example, if you created the outcome table 
#'                                         by finding the occurrence of the outcome 30 days after cohort start but
#'                                         you wish to conduct a sensitivity test for this and reduce it to 20 days,
#'                                         then set outcomeTime=20. 
#' @param newOutcome                       A vector of existing outcomeIds - constructs a new outcomeId based on people having all the specified outcomeIds.
#'                                         For example, you may wish to predict the subset of the people who have the outcome
#'                                         who also get given a secific treatment.  If you create the outcome with outcomeId 1
#'                                         and the treatment with outcomeId 2, the set newOutcome =c(1,2) to find all the people
#'                                         who have outcomeIds 1 and 2, they are then assigned a new outcomeId -1.               
#' @param predictionPeriod                 A vector of length 2 with the first value corresponding to the
#'                                         number of days after index to define the start of the risk prediction period
#'                                         and the second value corresponding to the number of days after index
#'                                         defining the end of the risk prediction period.  If this is NULL the the cohort
#'                                         start and end date will define the risk prediction period. 
#' @param dateInterval                     (a vector of 2 dates or NULL) corresponding to the inclusion dates. For
#'                                         example, if the user inputs c(1990-01-01,2000-01-01) then any people with an index
#'                                         prior to 1990 or after (Jan 1st 2000 minus minPriorObservaton) will be excluded. If
#'                                         NULL all dates are included.
#' @param minPriorObservation              an integer - people with the time in days between their observation start and cohort start 
#'                                         less than this number will be removed from the data.
#' @param minCohortTime                    an integer - people with the time in days between the cohort start and cohort end
#'                                         less than this number will be removed from the data.                                               
#' @param excludeOutcomeOccurrence         A list containing named list members and vectors of two integer values, where the name corresponds
#'                                         to the outcomeId and the integer vector corresponds to an interval in days for filtering people 
#'                                         who had the outcomeId recorded during this interval.  For example the list:
#'                                         list('1'=c(180,40), '4'=c('inf',0)) would find all the people who had the outcomeId 1 recorded 
#'                                         180 days prior to index up to 40 days after index and filter these people from the data, it would
#'                                         also find the people who have outcomeId 4 anytime prior to index and filter these people from the data.
#' @param classificationCensor             A list detailing the exclusion criteria for classification.  The list contains:
#'                                         \itemize{
#'                                         \item{insufficientCohortObservation - }{a character vector of length 2 with each element either 'include' 
#'                                         or 'exclude' - indicating whether to include or exclude patients whose observation period ends 
#'                                         before their cohort end.  The first element is applied to people with the outcome (class 1) and the second element is
#'                                         applied to people without the outcome (class 0)}
#'                                         \item{insufficientPredictionPeriod - }{a character vector of length 2 with each element either 'include' 
#'                                         or 'exclude' - indicating whether to include or exclude patients whose predictionPeriod falls outside of 
#'                                         their observation period.  The first element is applied to people with the outcome (class 1) and the second element is
#'                                         applied to people without the outcome (class 0)}
#'                                         \item{minPostObservation - }{An integer - specifying the required minimum number of days after index (Used by insufficientPostObservation) }
#'                                         \item{insufficientPostObservation - }{a character vector of length 2 with each element either 'include' 
#'                                         or 'exclude' - indicating whether to include or exclude patients with an index date plus minPostObservation greater than
#'                                         their observationEndDate.  The first element is applied to people with the outcome (class 1) and the second element is
#'                                         applied to people without the outcome (class 0)}
#'                                         
#'                                         }
#' @param  survivalCensor                 A list containing the criteria for censoring the data...                       
#' @examples 
#' # Filter any patients with an index before 2008-01-01 or after 2011-01-01
#' # and who have less than 365 days observation prior to index
#' plpData.censor censorPlpData(plpData, minPriorObservation = 365, 
#' dateInterval = c('2008-01-01','2011-01-01'))
#' 
#' # Filter patients with less than 100 days observtion prior to index
#' # also filter all people who are not observed for
#' # at least 100 days post index 
#' plpData.censor censorPlpData(plpData, minPriorObservation= 100, 
#' minCohortTime=NULL, 
#' classificationCensor=list( minPostObservation=100,
#'                          insufficientPostObservation = c('exclude','exclude')
#'                          )
#'                          )
#'                          
#' # Filter patients with less than 100 days observtion prior to index
#' # also filter people who do not have the outcome who are not observed for
#' # at least 100 days post index 
#' plpData.censor censorPlpData(plpData, minPriorObservation= 100, 
#' minCohortTime=NULL, 
#' classificationCensor=list(minPostObservation=100,
#'                          insufficientPostObservation = c('include','exclude')
#'                          )
#'                          )
#' # Filter people with an outcomeId 2 that occurs within 180 days before index
#' # until 5 days after index, also filter people with less than 365 days 
#' # observation prior to index and without a minimum of 365 days after index                         
#' plpData.censor censorPlpData(plpData, minPriorObservation= 365, 
#' excludeOutcomeOccurrence=list('2'=c(180,5)),
#' classificationCensor=list(minPostObservation=365,
#'                          insufficientPostObservation = c('exclude','exclude')
#'                          )
#'                          )                         
#' 
#' @return
#' An object of type \code{plpData} containing information on the prediction problem that only contains the
#' data satisfying the user's specified censoring options. This object will
#' contain the following data:
#' \item{cohorts}{An ffdf object listing all persons and their prediction periods. This
#' object will have these fields: row_id (a unique ID per period), person_id, cohort_start_date,
#' cohort_id, time (number of days in the window).} 
#' \item{outcomes}{An ffdf object listing all
#' outcomes per period. This object will have these fields: row_id, outcome_id, outcome_count,
#' time_to_event.} 
#' \item{exclude}{Either NULL or an ffdf object listing per outcome ID which windows
#' had the outcome prior to the window. This object will have these fields: rowId, outcomeId.}
#' \item{covariates}{An ffdf object listing the baseline covariates per person in the cohorts. This is
#' done using a sparse representation: covariates with a value of 0 are omitted to save space. The
#' covariates object will have three columns: rowId, covariateId, and covariateValue. }
#' \item{covariateRef}{An ffdf object describing the covariates that have been extracted.}
#' \item{metaData}{A list of objects with information on how the plpData object was constructed
#' and censoring details.  The list member named 'excluded' contains a ffdf of the excluded 
#' people and reason for exclusion.} 
#'
#' @export


censorPlpData <- function(plpData, outcomeIds=NULL, outcomeTime=NULL, newOutcome=NULL,
                          predictionPeriod =NULL,  dateInterval=NULL,
                          minPriorObservation= 365 #washoutWindow
                          , minCohortTime=NULL,
                          excludeOutcomeOccurrence=list('1'=c('inf',0)),
                          classificationCensor=list(insufficientCohortObservation = c('include','include'),
                                                    insufficientPredictionPeriod = c('include','include'),
                                                    minPostObservation=NULL,
                                                    insufficientPostObservation = c('include','include')
                          ),
                          survivalCensor=list(useCohortObservation = F,
                                              usePredictionPeriod = T,
                                              maxPostObservation=NULL,
                                              useMaxPostObservation = F)){
  
  cen.start <- Sys.time()
  covariates <- ff::clone(plpData$covariates)
  outcomes <- ff::clone(plpData$outcomes)
  cohorts <- ff::clone(plpData$cohorts)
  exclude.main <- NULL
  
  # create a new outcome based on people having all the newOutcome
  if(!is.null(newOutcome)){
    ppl <- NULL
    for (i in 1:length(newOutcome)){
      t <- outcomes$outcomeId==newOutcome[i]
      if(sum(t)>0){
        if(is.null(ppl))
          ppl <- ff::as.ram(unique(outcomes$rowId[ffwhich(t, t==T)]))
        if(!is.null(ppl))
          ppl <- intersect(ppl, ff::as.ram(unique(outcomes$rowId[ffbase::ffwhich(t, t==T)])))
      }
    }
    
    if(!is.null(ppl)){
      newOut <- ff::ffdf(rowId=ff::as.ff(ppl), outcomeId=ff::as.ff(rep(-1, length(ppl))), 
                     outcomeCount = ff::as.ff(rep(1, length(ppl))),
                     timeToEvent = ff::as.ff(rep(1, length(ppl))))
      writeLines(paste0('Added: ', nrow(newOut) , ' outcomes for -1'))
      outcomes <- ffbase::ffdfappend(outcomes, newOut)
    }
    
  }
  
  
  # filter the people w,ho don't have minimum cohort time 
  if(!is.null(minCohortTime)){
    writeLines(paste0('Filtering patients without sufficient cohort time of ', minCohortTime, ' days'))
    t <- cohorts$time < minCohortTime
    excluded < list()
    if(sum(t)>0)
      excluded <- cohorts[ffbase::ffwhich(t, t==T),]
    excluded$reason <- ff::ff(as.factor(rep('Insufficient cohort time',length(excluded$rowId))))
    if(is.null(exclude.main))
      exclude.main <- excluded
    if(!is.null(exclude.main) && !is.null(excluded) )
      exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
    if(sum(t==F)==0)
      stop('Minimum cohort time criteria has excluded everyone')
    if(sum(t==F)>0)
      cohorts <- cohorts[ffbase::ffwhich(t,t==F),]
    writeLines(paste0('Excluded ',sum(t), ' cohort rows'))
  }
  
  # filter the people w,ho don't have minPriorObservation
  if(!is.null(minPriorObservation)){
    writeLines(paste0('Filtering patients without sufficient observation of ', minPriorObservation, ' days'))
    t <- cohorts$priorIndexObs < minPriorObservation
    excluded < list()
    if(sum(t)>0)
      excluded <- cohorts[ffbase::ffwhich(t, t==T),]
    excluded$reason <- ff::ff(as.factor(rep('Insufficient history',length(excluded$rowId))))
    if(is.null(exclude.main))
      exclude.main <- excluded
    if(!is.null(exclude.main) && !is.null(excluded))
      exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
    if(sum(t==F)==0)
      stop('Criteria has excluded everyone after removing minimum observation time people - please choose different criteria')
    if(sum(t==F)>0)
      cohorts <- cohorts[ffbase::ffwhich(t,t==F),]
    writeLines(paste0('Excluded ',sum(t), ' cohort rows'))
  }
  
  # filter the people who are outside the dateInterval
  if(!is.null(dateInterval)){
    writeLines(paste0('Excluding people outside specified date of:', dateInterval[1],'-', dateInterval[2] ))
    dateInterval <- as.Date(dateInterval, format = "%Y-%m-%d")
    t <- as.Date(ff::as.ram(cohorts$cohortStartDate), format = "%Y-%m-%d") <= dateInterval[2] & 
         as.Date(ff::as.ram(cohorts$cohortStartDate), format = "%Y-%m-%d")  >= dateInterval[1]     
    t <- ff::as.ff(t)
    if(sum(t==F)>0) {
      excluded <- cohorts[ffbase::ffwhich(t, t==F),]
      excluded$reason <- ff::ff(as.factor(rep('Outside date rate',length(excluded$rowId))))
      if(is.null(exclude.main))
        exclude.main <- excluded
      if(!is.null(exclude.main))
        exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
      if(sum(t==T)>0)
        cohorts <- cohorts[ffbase::ffwhich(t,t==T),]
      if(sum(t==T)==0)
        cohorts <- NULL
      writeLines(paste0('Excluded ', nrow(excluded), ' cohort rows'))
    }
    if(sum(t)==0){
      stop('No People left after specified filtering- error at date interval - please revise exclusion criteria')
      #return(0)
    }
  }
  
  exclude <- ff::clone(plpData$exclude)
  # add people from exclusion using new time
  if(!is.null(outcomeTime)){
    writeLines(paste0('Adding outcomes that occured after ',outcomeTime, ' days'))
    # find the people in the exclusion with the certian outcomes after a specified time and add to 
    # plpData$outcomes
    
    t <- -1*exclude$time >= outcomeTime
    if(sum(t)>0){
      newOut <- exclude[ffbase::ffwhich(t, t==T),]
      
      newOut  <- aggregate(x = ff::as.ram(newOut$time), by = list(ff::as.ram(newOut$rowId), ff::as.ram(newOut$outcomeId)), FUN = "max")
      colnames(newOut ) <- c('rowId','outcomeId','time')
      newOut  <- ff::as.ffdf(test)
      newOut <- ff::ffdf(rowId=newOut$rowId, outcomeId=newOut$outcomeId,
                     outcomeCount = ff::ff(rep(1, length(newOut$rowId))), 
                     timeToEvent = newOut$time*-1 )
      
      outcomes <- ffbase::ffdfappend(outcomes, newOut)
      if(sum(t==F)>0)
        exclude <- exclude[ffbase::ffwhich(t, t==F),]
    }
  }
  
  
  # filter the people who have the outcome during the exclusion:
  # use exclude in plpData with rowId, time 

  if(length(excludeOutcomeOccurrence)>0 & !is.null(exclude)){
    writeLines('Excluding people based on user input time parameters...')
    #find the people who have the outcome in the excludeOutcomeOccurrence
    # for each cohort_concept in exclude ...
    for (conceptName in names(excludeOutcomeOccurrence)){
      upper <- ifelse(excludeOutcomeOccurrence[[conceptName]][2]=='inf', 100*365, excludeOutcomeOccurrence[[conceptName]][2] )
      lower <- ifelse(excludeOutcomeOccurrence[[conceptName]][1]=='inf', -100*365, -1*excludeOutcomeOccurrence[[conceptName]][1] )
      
      t <- (-1*exclude$time <= as.double(upper) & 
        -1*exclude$time >= as.double(lower)) &
        exclude$outcomeId == conceptName
      if(sum(t==T)>0){
        rowIds.exclude <- unique(exclude$rowId[ffbase::ffwhich(t,t==T)])
        t <- ffbase::ffmatch(cohorts$rowId, table=rowIds.exclude)
        excluded < list()
        if(sum(!is.na(t))>0)
          excluded <- cohorts[ffbase::ffwhich(t, !is.na(t)),]
        excluded$reason <- ff::ff(as.factor(rep('Prior history',length(excluded$rowId))))
      if(is.null(exclude.main))
        exclude.main <- excluded
      if(!is.null(exclude.main) && !is.null(excluded))
        exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
      if(sum(is.na(t))==0)
        stop('excluded all people - please revise exclusion criteria')
      cohorts <- cohorts[ffbase::ffwhich(t, is.na(t)),]
      
      # t <- ffmatch(exclude$rowId, table=rowIds.exclude)
      # exclude[ffwhich(t, is.na(t)),]
      writeLines(paste0('Removed ', length(rowIds.exclude), ' cohort rows with a history of outcomeId ',
                        conceptName,' between ',lower,' to ',upper, ' days relative to index'))
      }
    }
  }
  
  # now filter the outcomes based on unfiltered cohorts
  t <- ffbase::ffmatch(outcomes$rowId, table=cohorts$rowId)
  if(sum(!is.na(t))>0) outcomes <- outcomes[ffbase::ffwhich(t, !is.na(t)),]
  if(sum(!is.na(t))==0){
    #outcomes <- NULL
    stop('All outcomes are filtered... Please revise exclusion criteria')
    #return(0)
  }
  
  # now filter the classification
  if(length(classificationCensor)>0){
    if (!is.null(outcomeIds)){
      t <- ffbase::ffmatch(outcomes$outcomeId, table=outcomeIds)
      if(sum(!is.na(t))==0)
        stop('No outcome found for the specified outcomeIds')
      outcomes <- outcomes[ffbase::ffwhich(t, !is.na(t)),]
    }
    
    outcomeIds <- unique(outcomes$rowId)
    t <- ffbase::ffmatch(cohorts$rowId, table=outcomeIds)
    nonOutcomeIds <- c()
    if(sum(is.na(t))>0)
      nonOutcomeIds <- unique(cohorts$rowId[ffbase::ffwhich(t, is.na(t))])
    Ids <- list(outcomeIds,nonOutcomeIds)
    
    
    # do censoring for people without/with outcome
    for (ind in 1:2){
      
      # exclude based on cohort period outside observation
      if(classificationCensor$insufficientCohortObservation[ind]=='exclude'){
        # find people not observed for the whole cohort period
        t1 <- ffbase::ffmatch(cohorts$rowId, table=Ids[[ind]])
        t2 <- cohorts$cohortObsIncomplete==1
        t <- !is.na(t1)&t2
        if(sum(t)>0){
          excluded <- cohorts[ffbase::ffwhich(t, t==T),]
          excluded$reason <- ff::ff(as.factor(rep('Insufficient cohort observation',length(excluded$rowId))))
          if(is.null(exclude.main))
            exclude.main <- excluded
          if(!is.null(exclude.main))
            exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
        }
        if(sum(!t)==0)
          stop('excluded all non-outcome people - please revise exclusion criteria')
        
        if(sum(!t)>0){
          cohorts <- cohorts[ffbase::ffwhich(t, t==F),]
          writeLines(paste0('Excluded ',sum(t),
                            ' cohort rows (',ifelse(ind==1,'with','without'),' outcome) because insufficient cohort observation'))
        }
      }
      
      # exclude based on prediction period outside observation
      if(classificationCensor$insufficientPredictionPeriod[ind]=='exclude'){
        if(is.null(predictionPeriod) | length(predictionPeriod)!=2){warning('Using insufficientPredictionPeriod without correctly specifying predictionPeriod')}
        
        if(!is.null(predictionPeriod) & length(predictionPeriod)==2){
          # find people not observed for the whole risk period
          t1 <- ffbase::ffmatch(cohorts$rowId, table=Ids[[ind]])
          t2 <- cohorts$postIndexObs < predictionPeriod[2]
          t <- !is.na(t1)&t2
          
          if(sum(t)>0){
            excluded <- cohorts[ffbase::ffwhich(t, t==T),]
            excluded$reason <- ff::ff(as.factor(rep('Insufficient risk prediction period',length(excluded$rowId))))
            if(is.null(exclude.main))
              exclude.main <- excluded
            if(!is.null(exclude.main))
              exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
  
          writeLines(paste0('Excluded ',sum(t),
                            ' cohort rows (',ifelse(ind==1,'with','without'),' outcome) because insufficient risk prediction period'))
          
          }
          if(sum(!t)==0)
            stop('excluded all outcome people - please revise exclusion criteria')
          cohorts <- cohorts[ffbase::ffwhich(t, t==F),]
        }
      }
      
      # exclude based on minimum post observation
      if(classificationCensor$insufficientPostObservation[ind]=='exclude'){
        if(is.null(classificationCensor$minPostObservation)){warning('Using insufficientPostObservation without correctly specifying minPostObservation')}
        
        if(!is.null(classificationCensor$minPostObservation) & length(classificationCensor$minPostObservation)==1){
          # find people not observed for the whole risk period
          t1 <- ffbase::ffmatch(cohorts$rowId, table=Ids[[ind]])
          t2 <- cohorts$postIndexObs < classificationCensor$minPostObservation
          t <- !is.na(t1)&t2
          
          if(sum(t)>0){
            excluded <- cohorts[ffbase::ffwhich(t, t==T),]
            excluded$reason <- ff::ff(as.factor(rep('Insufficient post observation',length(excluded$rowId))))
            if(is.null(exclude.main))
              exclude.main <- excluded
            if(!is.null(exclude.main))
              exclude.main <- ffbase::ffdfappend(exclude.main, excluded)
          
          
          writeLines(paste0('Excluded ',sum(t),
                            ' cohort rows (',ifelse(ind==1,'with','without'),' outcome) because insufficient post observation'))
          }
          if(sum(!t)==0)
            stop('excluded all outcome people - please revise exclusion criteria')
          cohorts <- cohorts[ffbase::ffwhich(t, t==F),]
        }
      }
      
    }
    
    
  } # end classification filtering
  
  #now do survival censoring    
  if(length(survivalCensor)>0){
    #TODO
    
  }
  
 
  
  # Now the final filtering of the covarites and outcomes:
  
  includedIds <- unique(cohorts$rowId)
  t <- ffbase::ffmatch(covariates$rowId, table=includedIds)
  if(sum(!is.na(t))==0)
    stop('excluded all people with records - please revise exclusion criteria')
  covariates <- covariates[ffbase::ffwhich(t, !is.na(t)),]
  t <- ffbase::ffmatch(outcomes$rowId, table=includedIds)
  if(sum(!is.na(t))==0)
    stop('excluded all outcomes - please revise exclusion criteria')
  outcomes <- outcomes[ffbase::ffwhich(t, !is.na(t)),]
  
  metaData <- plpData$metaData
  metaData$censoring <- as.list(match.call())
  metaData$censoring$classificationCensor <- classificationCensor
  metaData$excluded <- exclude.main
  result <- list(cohorts = cohorts,
                 outcomes = outcomes,
                 covariates = covariates,
                 covariateRef = ff::clone(plpData$covariateRef),
                 metaData = metaData)
  
  class(result) <- "plpData"
  cen.time <-  Sys.time()-cen.start
  writeLines(paste0('Censoring took: ', format(cen.time, digits=3)))
  return(result)
  
  
 }



