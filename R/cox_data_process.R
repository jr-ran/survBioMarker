#' cox_data_process
#' data processing
#' cox input data
# 
#' @param data row is sample, col is phelotype variate,example age and AJCC stage
#' @param var.names phelotype variate
#' @return data.frame
#' @export
cox_data_process <- function(data, var.names){
  survival <- data[,c("status", "time")]
  out <- sapply(var.names[! var.names %in% c("status", "time")], function(x){
    if(!x %in% colnames(data)){
      return(NULL)
    }
    
    mydata <- data[,x]
    if((table(data[,x]) %>% length) == 2){
      return(mydata)
    }else if(is.numeric(data[,x])){
      return(mydata)
    }else{
      value <- vector()
      if(grepl("\\bage", x, ignore.case = T)){
        ind1 <- which(mydata >= 60)
        ind2 <- which(mydata < 60)
        value[ind1] <- "1"
        value[ind2] <- "0"
        return(value)
      }
      if(grepl("\\bajcc", x, ignore.case = T)){
        ind1 <- which(grepl(paste(prefix, "I{1,2}[A|B|C]?$", sep = ""), mydata, ignore.case = T, perl = T))
        ind2 <- which(grepl(paste(prefix, "(III[A|B|C]?)|(IV[A|B|C]?)$", sep = ""), mydata, ignore.case = T, perl = T))
        value[ind1] <- "Stage_I/II"
        value[ind2] <- "Stage_III/IV"
        return(value)
      }
      
    }
    
  })
  rownames(out) <- rownames(data)
  o <- cbind(survival, out)
  return(o)
}
