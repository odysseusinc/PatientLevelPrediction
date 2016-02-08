#' evaluatePlp
#'
#' @description
#' Evaluates the performance of the patient level prediction model
#' @details
#' The function calculates various metrics to measure the performance of the model
#' @param plpPredict                         The patient level prediction model's prediction
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param sparse                           (boolean) Whether the metrics should be calculated in sparse format
#'
#' @return
#' A list containing the performance values
#'

#' @export
evaluatePlp <- function(plpPredict, plpData, sparse=T ){
  predLab <- merge(as.data.frame(plpPredict)[, c('rowId', 'value')], ff::as.ram(plpData$outcomes), by='rowId', all.x=T )
  predLab[is.na(predLab)] <- 0
  
  if(sparse==F | nrow(plpPredict) <1000){
    roc <- pROC::roc(predLab$outcomeCount,predLab$value)
    
    lab.order <- predLab$outcomeCount[order(-predLab$val)]
    n <- nrow(predLab)
    P <- sum(predLab$outcomeCount>0)
    N <- n - P
    TP <- sapply(1:n, function(x) sum(lab.order[1:x]>0))
    FP <- sapply(1:n, function(x) sum(lab.order[1:x]==0))
    TN <- N-FP
    FN <- P-TP
    
    TPR <- TP/P
    FPR <- FP/N
    accuracy <- (TP+TN)/n
    PPV<- TP/(TP+FP)
    FOR <- FN/(FN+TN)
    Fmeasure <- 2*(PPV*TPR)/(PPV+TPR)
    
    metrics <- data.frame(TP,FP,TN,FN, TPR, FPR,PPV,FOR, accuracy, Fmeasure)
    
    #if(aveP==F)
    #  aveP <- NULL
    #if(aveP==T){
      P <- sum(predLab$outcomeCount>0)
      val <- rep(0, nrow(predLab))
      val[predLab$outcomeCount[order(-predLab$value)]>0] <- 1:P
      aveP <- sum(val/(1:nrow(predLab)))/P
    #}
    
    
    calPlot <- plotCalibration(plpPredict,
                               plpData,
                               removeDropoutsForLr = T,
                               numberOfStrata = 10,
                               truncateFraction = 0.01,
                               fileName = NULL)
    
    eval <- list(auc=roc$auc, aveP=aveP,
                 roc=roc,
                 raw = metrics[,c('TP','FP','TN','FN','FOR','accuracy')],
                 precision.recall = metrics[,c('TPR','PPV')],
                 F.measure = metrics[,c('Fmeasure')],
                 calPlot=calPlot)
    class(eval) <- 'metric.full'
  }
  if(sparse==T & nrow(plpPredict) >=1000 ){
    eval <- sparseMetric(plpPredict, plpData, predLab, aveP=T)
  }
  
 return(eval) 
}

sparseMetric <- function(prediction,plpData,predLab, aveP=T){
  writeLines(paste0('na: ', sum(is.na(prediction$value))))
  auc <- computeAuc(prediction,
                    plpData,
                    removeDropoutsForLr = T,
                    confidenceInterval = T)
  
  calPlot <- plotCalibration(prediction,
                             plpData,
                             removeDropoutsForLr = T,
                             numberOfStrata = 5,
                             truncateFraction = 0.01,
                             fileName = NULL)
  
  # now calculate 100 point tpr, fpr
  lab.order <- predLab$outcomeCount[order(-predLab$val)]
  n <- nrow(predLab)
  P <- sum(predLab$outcomeCount>0)
  N <- n - P
  change.ind <- which(lab.order==0)
  ind.oi <- change.ind[c(seq(1,length(change.ind)-1,
                             floor((length(change.ind)-1)/99)), length(change.ind))]
  TP <- sapply(ind.oi, function(x) sum(lab.order[1:x]>0))
  FP <- sapply(ind.oi, function(x) sum(lab.order[1:x]==0))
  TN <- N-FP
  FN <- P-TP
  
  TPR <- TP/P
  FPR <- FP/N
  accuracy <- (TP+TN)/n
  PPV<- TP/(TP+FP)
  FOR <- FN/(FN+TN)
  Fmeasure <- 2*(PPV*TPR)/(PPV+TPR)
  
  roc.sparse <- data.frame(TP,FP,TN,FN, TPR, FPR,PPV,FOR, accuracy, Fmeasure)
  
  aveP.val <- NULL
  if(aveP==T){
    val <- rep(0, n)
    val[lab.order>0] <- 1:P
    aveP.val <- sum(val/(1:n))/P
  }
  
  result <- list(auc=auc,aveP=aveP.val,
                 roc = roc.sparse[,c('FPR','TPR')],
                 raw = roc.sparse[,c('TP','FP','TN','FN','FOR','accuracy')],
                 precision.recall = roc.sparse[,c('TPR','PPV')],
                 F.measure = roc.sparse[,c('Fmeasure')],
                 calPlot =calPlot
  )
  class(result) <- 'metric.sparse'
  return(result)
  
}