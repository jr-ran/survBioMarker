#' getGEOProbeID
#' get gene information form GEO
#'
#'
#' @param mfile GEO family.soft
#' @param ... other parameter 
#'
#' 
#' @export
#'

getGEOProbeID <- function(mfile){
  ID <- readLines(file(mfile, "r"), 1000)
  titleL <- grep("^ID\t", ID, perl = T) 
  countL <- grep("!Platform_data_row_count", ID, perl = T) 
  len <- strsplit(ID[countL], split = "=")[[1]][2] %>% as.numeric()
  IDt <- readLines(file(mfile, "r"), n = titleL+len)
  #endL <- grep("!platform_table_end", ID, perl = T)
  IDr <- IDt[titleL:length(IDt)]
  od <- dirname(mfile)
  print(IDr[length(IDr)])
  write.table(file = file.path(od, "gene_ID.txt"), IDr, quote = F, row.names = F, col.names = F)
  cat("final file: ", file.path(od, "gene_ID.txt"), "\n")
  return("Done")
}