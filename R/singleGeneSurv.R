#' singleGS 
#' 
#' plot single g survival analysis,high:T3 ,low/in:T1-2
#' 
#' @param data clinical ana expresion data.frame
#' @param gene g name
#' @param out output dirname
#' @param method "quartile" or "median", use defined high/low
#' @param high.quantile if methos = "quartile", the quantile of high group
#' 
#' @return 
#' @export


singleGS <- function(data, gene, out, method = c("quartile", "median"), high.quantile = 2/3){
  lapply(gene, function(g){
    fdata <- data[,c("status", "time", g)]
    if(method == "quartile"){
      cut <- quantile(fdata[,g], high.quantile)
    }
    if(method == "median"){
      cut <- median(fdata[,g])

    }
    fdata$group <- NA
    fdata$group[fdata[,g] > cut] <- "high"
    fdata$group[fdata[,g] <= cut] <- "low"
    coxre <- uniCox(fdata)
    rownames(coxre)[rownames(coxre) == "group"] <- paste(g, "-low", sep = "")
    write.table(file = file.path(out, paste(g, "HR.xls", sep = ".")), coxre, quote = F, row.names = F, sep = "\t")
    m.surv <- surv_fit(Surv(time = time, event = status)~group,
                      data = fdata)
    
    p <- ggsurvplot(m.surv, pval = T,
                    legend.title = "",
                    title = g,
                    size = 1,
                    font.main = c(16, "bold", "darkblue"),
                    xlab="survival time",
                    ggtheme = theme_bw(),
                    risk.table = T,
                    risk.table.col = "group") 
    pdf(file.path(out, paste(g, "surv.pdf", sep = ".")))
    print(p)
    dev.off()
    return(paste(g, "is done", sep = " "))
  })
}