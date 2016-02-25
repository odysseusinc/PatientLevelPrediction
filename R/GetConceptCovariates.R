# @file GetConceptCovariates.R
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

#' Get HDPS covariate information from the database
#'
#' @description
#' Constructs the set of covariates for one or more cohorts using data in the CDM schema. This
#' implements the covariates based on user input concept sets.  This is a way to incorporate known risk factors.
#'
#' @param covariateSettings   An object of type \code{covariateSettings} as created using the
#'                            \code{\link{createConceptCovariateSettings}} function.
#'
#' @template GetCovarParams
#'
#' @export
getDbConceptCovariateData <- function(connection,
                                   oracleTempSchema = NULL,
                                   cdmDatabaseSchema,
                                   cdmVersion = "5",
                                   cohortTempTable = "cohort_person",
                                   rowIdField = "subject_id",
                                   covariateSettings) {
  writeLines("Constructing concept set covariates")
  if (substr(cohortTempTable, 1, 1) != "#") {
    cohortTempTable <- paste("#", cohortTempTable, sep = "")
  }
  cdmDatabase <- strsplit(cdmDatabaseSchema, "\\.")[[1]][1]
  
  if (cdmVersion == "4") {
    cohortDefinitionId <- "cohort_concept_id"
    conceptClassId <- "concept_class"
    measurement <- "observation"
  } else {
    cohortDefinitionId <- "cohort_definition_id"
    conceptClassId <- "concept_class_id"
    measurement <- "measurement"
  }
  
    # insert conceptSets table (covariateId,conceptId, daysPrior) 
    DatabaseConnector::insertTable(connection,
                                   tableName = "#risk_factors",
                                   data = covariateSettings$conceptSets,
                                   dropTableIfExists = TRUE,
                                   createTable = TRUE,
                                   tempTable = TRUE,
                                   oracleTempSchema = oracleTempSchema)
  
  renderedSql <- SqlRender::loadRenderTranslateSql("GetConceptCovariates.sql",
                                                   packageName = "PatientLevelPrediction",
                                                   dbms = attr(connection, "dbms"),
                                                   oracleTempSchema = oracleTempSchema,
                                                   cdm_database = cdmDatabase,
                                                   cdm_version = cdmVersion,
                                                   cohort_temp_table = cohortTempTable,
                                                   row_id_field = 'row_id',
                                                   cohort_definition_id = cohortDefinitionId,
                                                   concept_class_id = conceptClassId,
                                                   measurement = measurement,
                                                   use_measurement=cdmVersion=="5",
                                                   use_covariate_demographics=covariateSettings$useDemo
                                                   )
  
  DatabaseConnector::executeSql(connection, renderedSql)
  writeLines("Done")
  
  writeLines("Fetching data from server")
  start <- Sys.time()
  covariateSql <- "SELECT row_id, covariate_id, covariate_value FROM #cov ORDER BY covariate_id, row_id"
  covariateSql <- SqlRender::renderSql(covariateSql, cohort_definition_id = cohortDefinitionId)$sql
  covariateSql <- SqlRender::translateSql(covariateSql,
                                          "sql server",
                                          attr(connection, "dbms"),
                                          oracleTempSchema)$sql
  covariates <- DatabaseConnector::querySql.ffdf(connection, covariateSql)
  covariateRefSql <- "SELECT covariate_id, covariate_name, analysis_id, concept_id  FROM #cov_ref ORDER BY covariate_id"
  covariateRefSql <- SqlRender::translateSql(covariateRefSql,
                                             "sql server",
                                             attr(connection, "dbms"),
                                             oracleTempSchema)$sql
  covariateRef <- DatabaseConnector::querySql.ffdf(connection, covariateRefSql)
  
  sql <- "SELECT COUNT_BIG(*) FROM @cohort_temp_table"
  sql <- SqlRender::renderSql(sql, cohort_temp_table = cohortTempTable)$sql
  sql <- SqlRender::translateSql(sql,
                                 targetDialect = attr(connection, "dbms"),
                                 oracleTempSchema = oracleTempSchema)$sql
  populationSize <- DatabaseConnector::querySql(connection, sql)[1, 1]
  
  delta <- Sys.time() - start
  writeLines(paste("Loading took", signif(delta, 3), attr(delta, "units")))
  
  renderedSql <- SqlRender::loadRenderTranslateSql("RemoveCovariateTempTables.sql",
                                                   packageName = "PatientLevelPrediction",
                                                   dbms = attr(connection, "dbms"),
                                                   oracleTempSchema = oracleTempSchema)
  DatabaseConnector::executeSql(connection,
                                renderedSql,
                                progressBar = FALSE,
                                reportOverallTime = FALSE)
  
  colnames(covariates) <- SqlRender::snakeCaseToCamelCase(colnames(covariates))
  colnames(covariateRef) <- SqlRender::snakeCaseToCamelCase(colnames(covariateRef))
  
  metaData <- list(sql = renderedSql,
                   call = match.call(),
                   conceptSets=covariateSettings$conceptSets)
  result <- list(covariates = covariates, covariateRef = covariateRef, metaData = metaData)
  class(result) <- "covariateData"
  return(result)
}


#' Create Concept covariate settings
#'
#' @details
#' creates an object specifying how covariates should be contructed from data in the CDM model.
#'
#' @param conceptList   A list of lists - each inner list contains two objects: conceptSet a vector of conceptSetIds 
#'                      and prior an integer specifying the number of days prior to index to search for the concepts in the set
#'                        
#' @return
#' An object of type \code{conceptCovariateSettings}, to be used in other functions.
#'
#' @export
createConceptCovariateSettings <- function(conceptList,useDemo=TRUE) {

  #extract table of concept sets
  #CONCEPT_ID, CONCEPT_NAME, STANDARD_CONCEPT, INVALID_REASON,
  #CONCEPT_CODE, DOMAIN_ID, VOCABULARY_ID, CONCEPT_CLASS_ID,
  #STANDARD_CONCEPT_CAPTION, INVALID_REASON_CAPTION,
  #isExcluded, includeDescendants, includeMapped
  
  # extract conceptsets and create table
  tempConcept <- lapply(conceptList, function(x) loadConceptSetIncludedConcepts(x$concepts)  )
  conceptTable <- c()
  for(i in 1:length(tempConcept)){
    conceptSets <- data.frame(COVARIATE_ID=rep(-i, nrow(tempConcept[[i]])), 
                              CONCEPT_ID=tempConcept[[i]][,'CONCEPT_ID'],
                              CONCEPT_NAME=tempConcept[[i]][,'CONCEPT_NAME'], 
                              isExcluded=tempConcept[[i]][,'isExcluded'], 
                              includeDescendants=tempConcept[[i]][,'includeDescendants'], 
                              includeMapped=tempConcept[[i]][,'includeMapped'], 
                              priorTime=rep(conceptList[[i]]$time,nrow(tempConcept[[i]])),
                              covariate_name=rep(conceptList[[i]]$name,nrow(tempConcept[[i]]))
    )
    conceptTable <- rbind(conceptTable, conceptSets)
  }
  rownames(conceptTable) <- NULL
  
  covariateSettings <- list(conceptSets = conceptTable,
                            useDemo=useDemo)
  
  attr(covariateSettings, "fun") <- "getDbConceptCovariateData"
  class(covariateSettings) <- "covariateSettings"
  return(covariateSettings)
}


#exclusion concept IDs
loadConceptSetIncludedConcepts <- function(conceptSetIdentifiers) {
  allConcepts <- c()
  for(conceptSetIdentifier in conceptSetIdentifiers){
    conceptSetExpressionUrl = paste("http://hix.jnj.com:8080/WebAPI/conceptset/", conceptSetIdentifier, "/expression", sep = "")
    expression <- httr::content(httr::GET(conceptSetExpressionUrl),"text")
    expression <- RJSONIO::fromJSON(expression)
    concepts <- do.call(rbind, lapply(expression$items, function(x) c(x$concept$CONCEPT_ID, x$concept$CONCEPT_NAME,
      x$isExcluded,x$includeDescendants,x$includeMapped)))
    colnames(concepts) <- c('CONCEPT_ID','CONCEPT_NAME','isExcluded','includeDescendants','includeMapped')
    
    allConcepts <- rbind(allConcepts, concepts)
  }
  return(allConcepts)
}

