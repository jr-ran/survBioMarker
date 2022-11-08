#' ExpDataProcess
#' convert gene symbol and merge data
#'
#'
#' @param data1 data.frame,geneID and gene expression 
#' @param data2 data.frame, gene information
#' @param id.col set column, example:1,2 or "ID,symbol"
#' @param ... other parameter 
#'
#' @return expression data.frame
#' @export
#'

ExpDataProcess <- function(data1, data2 = NULL, id.col = NULL,...){
  .ConvertID <- function(data1, EnseID, id.col,...){
    if(!is.null(EnseID) && is.null(id.col)){
      print("Error:: please set id.col !!!")
      q()
    }
    EnsemblID <- EnseID[,unlist(strsplit(id.col, split = ",")) %>% as.numeric()]
    colnames(EnsemblID) <- c("ID", "symbol")
    colnames(data1)[1] <- "ID"
    mat <- inner_join(EnsemblID, data1, by = "ID")
    mat <- mat[,-1]
    return(mat)
  }
  rawmat <- data1
  if(!is.null(data2)){
    rawmat <- .ConvertID(data1, EnsemblID, id.col)
  }
  if(length(unique(rawmat[,1])) == nrow(rawmat)){
    uMat <- rawmat
    rownames(uMat) <- uMat[,1]
    uMat <- uMat[,-1]
  }else{
    uMat <- aggregate(rawmat[,-1], by = list(rawmat[,1]), FUN = mean, na.rm = T) 
    rownames(uMat) <- uMat$Group.1
    uMat <- uMat[,-1]
  }
  uMat <- uMat %>% as.matrix()
  return(uMat)
}