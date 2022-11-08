#' getTumorID: tumor sample ID list
#' TCGA sample name, "-" split
#'
#' @param sample is vector
#' @param ... other parameters  
#' @return a vector
#' @rdname getTumorID
#' @export

getTumorID <- function(sample, ...){
  tumor.sample <- sapply(sample, function(x){
    if(!grepl("^TCGA", x, perl = T)){
      return(x)
    }else{
      number <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x),split="-"))[4], split = "[A-Z]"))[1])
      if(number <= 9){
        id <- as.character(x)
        return(id)
      }
    }
    
  })
  tumor.sample <- unlist(tumor.sample)
  return(tumor.sample)
}
