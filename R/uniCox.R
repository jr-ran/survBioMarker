#' Cox proportional hazards regression analysis
#'
#' @param data row is TCGA sample,col is variable
#' @param var variable name,example sex,age,...
#' @param method: univariate or multivariate
#' @param ... other parameter 
#'
#' @export
#'

uniCox <- function(data, var = NULL, method = "univariate"){
  if(method == "univariate"){
    ifelse(is.null(var),mvar <- colnames(data)[3:ncol(data)],mvar <- var)
    out <- sapply(mvar, function(x){
      mydata <- data[!is.na(data[,x]),]
      n <- paste0("`", x, "`")
      tmp <- coxph(as.formula(paste0('Surv(time = time, event = status)~', n)), data = mydata)
      tmp2 <- summary(tmp)
      t1 <- tmp2$coef[,c(1,5,2)]
      t2 <- tmp2$conf.int[,c("lower .95", "upper .95")]
      o1 <- c(rownames(tmp2$coef), t1,t2)
      o1
    })
    out <- t(out) %>% as.data.frame
    colnames(out) <- c("variable", "coef", "pvalue", "HR", "lower.95", "upper.95")
  }else{
    var <- paste0("`", var, "`")
    mformula <- as.formula(paste0("Surv(time = time, event = status)~", paste(var, collapse = "+")))
    tmp <- coxph(mformula,data)
    tmp2 <- summary(tmp)
    t1 <- tmp2$coef[,c(5,2)]
    t2 <- tmp2$conf.int[,c("lower .95", "upper .95")]
    out <- cbind(t1,t2)
    colnames(out) <- c("pvalue", "HR", "lower.95", "upper.95")
  }
  out
}

