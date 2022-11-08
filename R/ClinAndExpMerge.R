#' MainRun
#' find signature and survival analysis 
#'
#' @param expMat survival and expression data, col is status,time,g1,g2,...
#' @param clinData clinical survival data, must contain two columns: survival status and survival time
#' @param surv.col columns of status and time, default c("status", "OS")
#' @param ... other parameter 
#'
#'
#' @return data.frame, Contains survival and expression data
#' @export
#'


ClinAndExpMerge <- function(expMat, clinData, surv.col = c("status", "OS"),...){
  
  clinData <- clinData[,surv.col]
  sample <- rownames(clinData)
  clinData <- apply(clinData,2,as.numeric)
  rownames(clinData) <- sample
  colnames(clinData) <- c("status", "time")
  # Selecting tumor samples 
  tumorID <- getTumorID(sample = rownames(clinData))
  fsample <- intersect(colnames(expMat), tumorID)
  tumor.mat <- expMat[,fsample] %>% t()
  tumor.survival <- clinData[rownames(clinData) %in% fsample,]
  colnames(tumor.survival)[1:2] <- c("status", "time")
  
  # integrating expression and clinical information
  initiData <- data.frame(tumor.survival, tumor.mat[rownames(tumor.survival),], check.names = FALSE)
  initiData <- initiData %>% filter(!is.na(time)) %>% filter(!is.na(status))
  return(initiData)
}