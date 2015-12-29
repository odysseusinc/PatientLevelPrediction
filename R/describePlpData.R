#' describePlpData
#'
#' @details
#' Describes the data
#'
#' @param plpData               An object of type \code{plpData}.
#' @param covariateVals         A vector of the covariateIds to focus on
#' @param cdmDatabase           The database the data comes from
#' @param outcomeId             The outcome of interest
#' @param agehist               Whether to plot a histogram of ages
#' @param plot                  Whether to plot results
#' @param plotFile              Location to save plot 
#' @param saveTable             Whether to save results
#' @param tableFile             Location to save results table 
#' @param outcomeName           Outcome name to be used on plots
#' @param cohortName            Cohort name to be used on plots
#' @param perYear               Whether to do descriptive stats per year
#' 
#' @export

# plots description of censored data and saves into csv
describePlpData <- function(plpData, covariateVals=NULL, cdmDatabase,
                            outcomeId =2,
                            agehist=TRUE,
                            plot=T, plotFile=NULL, saveTable=T, tableFile=NULL,
                            outcomeName='Pregnancy with GDM',
                            cohortName='Total pregnancy',
                            perYear=T){
  require(ggplot2)
  outcome <- ff::clone(plpData$outcomes)
  t <- outcome$outcomeId==outcomeId
  outcome <- outcome[ffbase::ffwhich(t,t==T),]
  # TODO add age mean/st_dev into yearSummary
  covSummary <- c()
  # summarise the data
  yearSummary <- data.frame(
  'year' = c('All years:'),
  'Total Cohort' = c(nrow(plpData$cohorts)),
  'Outcome Count'=length(unique(outcome$rowId)),
  'Prevalance'=length(unique(outcome$rowId))/nrow(plpData$cohorts))
  
  
  # calculate: year, people in cohort, outcome count, prevalance
  if(perYear==T){
    all <- merge(plpData$cohorts, outcome, by='rowId', all.x=T)
    
    years <- sapply(X=ff::as.ram(all$cohortStartDate), FUN=function(x) substring(as.character(x),1,4))
    years <- as.double(years)
    
    lab <- rep(1, length(ff::as.ram(all$outcomeCount)))
    lab[is.na(ff::as.ram(all$outcomeCount))] <- 0
    
    out.years <- data.frame(years, lab)
 
    GDM <- aggregate(out.years$lab, list(year= out.years$years), sum)
    colnames(GDM)[2] <- outcomeName
    Total <- aggregate(rep(1,length(out.years$lab)), list(year= out.years$years), sum)
    colnames(Total)[2] <- cohortName
    out.years <- merge(GDM, Total, by='year', all=T)
    out.years$prevalance <- out.years[,2]/out.years[,3]
  
    yearSummary2 <- data.frame(
      'year' = as.character(out.years$year),
      'Total Cohort' = out.years[,3],
      'Outcome Count'= out.years[,2],
      'Prevalance'=out.years[,4]
    )
    
    yearSummary <- rbind(yearSummary, yearSummary2)
      
  }
  
  # plot prevalence per year
  if(plot){
    if(!is.null(plotFile)) 
      pdf(file=file.path(plotFile, paste0(cdmDatabase, 'prev_byyear.pdf')), width = 14, height = 9)
    print(ggplot2::ggplot(data=out.years, aes(x=as.factor(year), y=prevalance)) +
            geom_bar(stat="identity") +
            xlab("Year") + ylab("Prevalence") +
            ggtitle(paste0("Outcome prevalence by year for ",cdmDatabase )) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) )
    if(!is.null(plotFile)) 
      dev.off()
    }
    
  
  
  # plot covariates frequencies in total, outcome and non-outcome groups
  if(!is.na(covariateVals)){
    covref <- ff::clone(plpData$covariateRef)
    if('riskFactors'%in%covariateVals) # add covariates
      covariateVals <- c(covariateVals[!covariateVals%in%'riskFactors'],-2:-16)
    if('age'%in%covariateVals) # add the ages
      covariateVals <- c(covariateVals[!covariateVals%in%'age'], covref$covariateId[substring(ff::as.ram(covref$covariateName),1,4)=='Age '])
      

    gdmCount <- nrow(outcome)
    allCount <- nrow(plpData$cohorts)
    
    outcomes <- merge(plpData$cohorts, outcome, by='rowId', all.x=T)
    t <- is.na(outcomes$outcomeCount)
    gdm <- outcomes[ffbase::ffwhich(t, t==F),]
    
    covariates <- ff::clone(plpData$covariates)
    t <- ffbase::ffmatch(covariates$covariateId, 
                 table= ff::ff(as.double(covariateVals))
                 )
    covariates <- covariates[ffbase::ffwhich(t, !is.na(t)),]
    covariates <- merge(covariates, ff::as.ffdf(covref), by='covariateId', all.x=T)
    
    all.gdm <- merge(covariates, gdm, by='rowId')
    all.ppl <- merge(covariates, outcomes, by='rowId')
  
    GDM <- aggregate(ff::as.ram(all.gdm$outcomeCount>0), list(covariateId= ff::as.ram(all.gdm$covariateName)), sum)
    colnames(GDM)[2] <- outcomeName
    GDM[,2] <- GDM[,2]/gdmCount
    allPlp <- aggregate(rep(1,length(all.ppl$outcomeCount)), list(covariateId= ff::as.ram(all.ppl$covariateName)), sum)
    colnames(allPlp)[2] <- 'All People'
    allPlp[,2] <- allPlp[,2]/allCount
    out.dat <- merge(GDM, allPlp, by='covariateId', all=T)
    out.dat[is.na(out.dat)] <- 0
    out.dat2 <- reshape2::melt(out.dat, 'covariateId')
    
    
    # plot the results:
    if(plot){
      if(!is.null(plotFile))
        pdf(file=file.path(plotFile, paste0(cdmDatabase, '_covariates.pdf')), width = 14, height = 9)
      print(ggplot2::ggplot(data=out.dat2, aes(x=as.factor(covariateId), y=value, fill=variable)) +
              geom_bar(stat="identity", position=position_dodge()) +
              xlab("") + ylab("Frequency") +
              facet_wrap(~covariateId, scales="free") +
              theme(axis.line=element_blank(),axis.text.x=element_blank()) +
              ggtitle(paste0("Covariate risk prevalence for ",cdmDatabase ))  )
      if(!is.null(plotFile))
        dev.off()
    }
    
    # now age distribution
    if(agehist==T){ # get continuous age
      covariateVals <- covref$covariateId[ff::as.ram(covref$covariateName)=='Age']
      covariates <- ff::clone(plpData$covariates)
      t <- ffbase::ffmatch(covariates$covariateId, 
                   table= ff::ff(as.double(covariateVals))
      )
      covariates <- covariates[ffbase::ffwhich(t, !is.na(t)),]
      covariates <- merge(covariates, ff::as.ffdf(covref), by='covariateId', all.x=T)
      
      all.gdm <- merge(covariates, gdm, by='rowId')
      all.ppl <- merge(covariates, outcomes, by='rowId')
      
      writeLines(paste0('gdm:',nrow(all.gdm), '- all:',nrow(all.ppl) ))
      
      ageHist <- rbind(
      data.frame(age=ff::as.ram(all.gdm$covariateValue), 
                 type=rep('Gestational Diabetes', length(all.gdm$covariateValue))),
      data.frame(age=ff::as.ram(all.ppl$covariateValue), 
                 type=rep('All people', length(all.ppl$covariateValue)))
      )
      writeLines(paste0('ageHist row:',nrow(ageHist), '- col:',ncol(ageHist) ))
      
      
      if(!is.null(plotFile))
        pdf(file=file.path(plotFile, paste0(cdmDatabase, '_age_hist.pdf')), width = 14, height = 9)
      print(ggplot2::ggplot(data=ageHist, aes(x=age, fill=type)) + 
        #geom_histogram() +
        #facet_wrap(~type, scales="free") +
        labs(x="Standardized Age", y="Density") +
        #geom_histogram(aes(y = ..density..), binwidth = 0.02, position="dodge") + 
        geom_density() + geom_density(alpha = 0.3)
        #stat_bin(aes(y=..count../sum(..count..)))
      )
      if(!is.null(plotFile))
        dev.off()
    
    }
    
  }
  
  # save table
  if(saveTable & !is.null(tableFile)){
    write.csv(yearSummary, file.path(tableFile,paste0(cdmDatabase,'yearSummary.csv')))
    write.csv(out.dat, file.path(tableFile,paste0(cdmDatabase,'covSummary.csv')))
  }
  # return table
  results <- list(yearSummary=yearSummary,
                  covSummary = out.dat
                  )
  return(results)
}