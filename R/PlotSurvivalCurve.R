#' PlotSurvivalCurve 
#' 
#' plot survival-curve picture
#' 
#' @param survival_plot survfit object
#' @param dir output dirname
#' @param xname xlab name
#' default OS
#' @param title title 
#' @return p,plot object
#' @export

PlotSurvivalCurve <- function(survival_plot, dir, xname = "OS", title = NULL){
  p1 <- ggsurvplot(survival_plot, pval = T,
                  legend.title = "",
                  font.main = c(16, "bold", "darkblue"),
                  #censor.shape = 26,
                  palette = c("red","blue"),
                  #censor=F,
                  xlab = xname,
                  title = title,
                  ggtheme = theme_bw(),
                  #risk.table.col = "riskGroup",
                  risk.table = T) 

  ifelse(is.null(title),filen <- "survival.pdf",filen <- paste(title, "survival.pdf", sep = "."))
  pdf(file.path(dir, filen), width=9)
  print(p1)
  dev.off()
  
  #png(paste(dir, paste(title, "survival.png", sep = "."), sep = "/"), width=680, height=480)
  #print(p1)
  #dev.off()
  return(p1)
}