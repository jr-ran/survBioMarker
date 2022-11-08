#' GetFilterExpression
#' filter sample by gene expression percentage
#'
#' @param data TCGA expression data, row is gene,column is sample
#' @param percent gene expression percentage
#' default 0.1
#' @param ... other parameters  
#'
#'
#' @return expression matrix
#' @rdname GetFilterExpression
#' @export

GetFilterExpression <- function(data, percent = 0.1){
  gratio <- .GetSampleRate(data)
  gnames <- .GetFilternames(gratio, x = percent)
  filter_exp <- data[rownames(data) %in% gnames,]
  return(filter_exp)
}

# statistic gene expression ratio
# data: expression matrix
.GetSampleRate <- function(data){
  # sample
  n <- ncol(data)
  filter.gene <- function(x){
    index <- which(x > 0) 
    rate <- length(index)/n 
    return(rate)
  }
  en <- apply(data, 1, filter.gene)
  return(en)
}

# data: function .GetSampleRate() result
# x: Minimum expression percentage 
# return genes, expression ratio > x%
.GetFilternames <- function(data, x = 0.1){
  data1 <- which(data > x)
  gnames <- names(data1)
  return(gnames)
}

