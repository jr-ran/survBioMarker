#' getGEOclinical
#' get clinical information form GEO
#'
#'
#' @param mfile GEO series_matrix.txt
#' @param ... other parameter 
#'
#' 
#' @export
#'

getGEOclinical <- function(mfile){
  od <- dirname(mfile)
  keyLine <- readLines(file(mfile, "r"), 1000)
  index <- grep("!Sample_geo_accession|!Sample_characteristics", keyLine, perl = T)
  if(length(index) <= 1){
    cat("Warnning:: there is no clinical information in file:", mfile, "\n")
    return(NULL)
  }
  rawClin <- keyLine[index]
  write.table(file = file.path(od, "raw_clinical_info.txt") , rawClin, quote = F, sep = "\t", row.names = F, col.names = F)
  rawC <- read.delim(file.path(od, "raw_clinical_info.txt"))
  rawC <- rawC[,-1] %>% t()
  chs <- list()
  for (i in 1:ncol(rawC)) {
    x <- rawC[,i]
    x <- gsub("[(]|[)]", "", x)
    ch <- strsplit(x[1], split = ": ") %>% unlist %>% .[1]
    r <- gsub(paste(ch, ": ", sep = ""),"", x)
    if((grep("^[0-9]+$", r, perl = T, ignore.case = T) %>% length) > 1){
      r <- as.numeric(r)
    }
    chs[[ch]] <- r
  }
  clinical <- do.call(data.frame,chs)
  write.table(file = file.path(od, "clinical_info.txt"), clinical, quote = F, sep = "\t", row.names = T)
  system(paste("rm ", file.path(od, "raw_clinical_info.txt")))
  cat("clinical information path:", file.path(od, "clinical_info.txt"), "\n")
  return("Done")
}