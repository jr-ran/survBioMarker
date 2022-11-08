#' MainRiskSurv
#' 
#' 
#' define risk group and survival for choose feature
#' 
#' @param data clinical and expression data
#' @param choose_feature signature gene set
#' @param od output dirname
#' @param ... other parameters 
#' @return data.frame
#' 
#'
#' @rdname MainRiskSurv
#' @export 

MainRiskSurv <- function(data, choose_feature, od){
  if(length(choose_feature) < 2){
    print("Error:: choose_feature must contains at least two genes!")
    return(NA)
  }
  sdata <- data[,colnames(data) %in% c("status", "time", choose_feature)]
  PI_data <- RiskScore(sdata, feature = choose_feature, method = "median")
  print(choose_feature)
  print(head(PI_data))
  finalData <- cbind(sdata[,1:2], PI_data)
  var <- "risk_group"
  KMSurvival(finalData, var, od, key)
  coxData <- cox_data_process(finalData, var.names = var)
  uniValid <- uniCox(coxData, var = var, method = "univariate")
  uniValid$coef <- NULL
  write.table(file = file.path(od, "univariate_HR.xls"), uniValid, quote = FALSE, sep = "\t", row.names = FALSE)
}

