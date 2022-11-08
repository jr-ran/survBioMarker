#' auc
#' model auc
#'
#' @param response response value
#' @param predictor predictor value
#' @param od output dirname
#' @param ... other parameter 
#'
#'
#' @return roc object
#' @export
#'

auc <- function(response, predictor, od,...){
  pdf(file.path(od, "roc.pdf"))
  rocobj <- roc(response = response, predictor = predictor, plot = TRUE, print.auc=TRUE)
  dev.off()
  return(rocobj)
}
