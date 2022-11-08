#' gdcSampleID
#' Mapping between file names and sample names
#'
#'
#' @param mjson GDC metadata.cart json
#' @param ... other parameter 
#'
#' 
#' @export
#'

gdcSampleID <- function(mjson){
  metadata <- jsonlite::fromJSON(mjson)
  sids <- sapply(metadata$associated_entities, function(x){x[,c(1,3)]}) 
  filenames <- metadata$file_name 
  sids2 <- apply(t(sids), 2, function(x) {unlist(x)})
  id.data <- data.frame(sids2, file = filenames) 
  return(id.data)
}