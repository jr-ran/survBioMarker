#' RiskScore
#' 
#' calculate prognostic index(PI) or risk score
#' risk_score = Beta(1)Exp(1)+Beta(2)Exp(2)+...+Beta(p)Exp(p)
#' 
#' @param data TCGA expression data
#' @param feature choose feature,example g1,g2,g3
#' @param ... other parameters 
#' @return data.frame
#' 
#'
#' @rdname risk_score
#' @export 

RiskScore <- function(data, feature, method = "median"){
  feature <- intersect(feature, colnames(data))
  data1 <- data[,colnames(data) %in% c("status", "time", feature)]
  coxre <- uniCox(data1)
  coefL <- coxre$coef %>% as.numeric
  PI <- apply(data1[,feature], 1, function(x){
    val <- sum(x*coefL)
  })
  if(method == "median"){
    cutVal <- median(PI)
  }
  highIndex <- which(PI > cutVal)
  lowIndex <- which(PI <= cutVal)
  # risk group, high "1", low "0"
  risk_score <- PI
  risk_score[highIndex] <- "high"
  risk_score[lowIndex] <- "low"
  names(risk_score) <- rownames(data1)
  re <- data.frame(risk_score = PI, risk_group = risk_score)
  return(re)
}

