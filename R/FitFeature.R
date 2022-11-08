#' LASSO fit feature
#' use LASSO regression optimize feature
#' 
#' @param x feature matrix
#' @param y observations, survival status
#' @param od output dirname
#' @param time run times
#' default 1
#' @param ... other parameters 
#'
#' @return choose_feature
#'
#' @export
LfitFeature <- function(x, y, od = NULL, time = 1, ...){
  model_lasso <- glmnet(x, y, family = "binomial", alpha = 1)
  cv_fit <- cv.glmnet(x=x, y=y, nlambda = 100, alpha = 1)
  if(!is.null(od)){
    pdf(file.path(od, paste(time, "lasso_model.pdf", sep = ".")))
    plot(model_lasso, xvar = "lambda", label = FALSE)
    dev.off()
    pdf(file.path(od, paste(time, "lambda_model.pdf", sep = ".")))
    plot(cv_fit)
    dev.off()
  }
  optimal_model_lasso <- glmnet(x=x, y=y, family = "binomial", alpha = 1, lambda = cv_fit$lambda.min)
  choose_feature <- rownames(optimal_model_lasso$beta)[as.numeric(optimal_model_lasso$beta)!=0]
  return(choose_feature)
}
