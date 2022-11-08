#' KMSurvival 
#' Run Kaplan-Meier survival analysis for biomarker
#' 
#' @param object survival data,first column is status,second column is time,...
#' @param var.name colnames for survfit
#' @param od output dirname
#' @param ... other parameter
#'
#'
#' @return data.frame
#' @export

KMSurvival <- function(object, var.name, od, xname = "OS",...){
  f <- as.formula(paste("Surv(time = time, event = status)~", var.name, sep = ""))
  surv_plot <- surv_fit(f, data = object)
  #surv_plot <- survfit(Surv(time = time, event = status)~risk_group,
  #                   data=data)
  PlotSurvivalCurve(surv_plot, dir = od, xname = "OS")
  return("K-M done")
}

