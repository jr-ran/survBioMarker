#' MainFeatureChoose
#' Prognostic markers
#' 
#' @param data clinical and expression data
#' @param geneS candidate gene set
#' @param od output dirname
#' @param ... other parameters 
#' @return data.frame
#' 
#'
#' @rdname MainFeatureChoose
#' @export 

MainFeatureChoose <- function(data, geneS, od){
  inputData <- data[,colnames(data) %in% c("status", "time", geneS)]
  mgs <- survivalBioMarker::OptimizeFeature(inputData)
  write.table(file = file.path(od, "final_feature_set.txt"), mgs, quote = FALSE, sep = "\t", row.names = FALSE)
}

