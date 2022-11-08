#' gdcExp
#' get expression matrix form GDC
#'
#'
#' @param idir GDC case expreesion path
#' @param mjson GDC metadata.cart json
#' @param ... other parameter 
#'
#' 
#' @export
#'
gdcExp <- function(idir, mjson){
  
  sid.data <- gdcSampleID(mjson)
  # read expression
  allfiles <- dir(idir, pattern = "*.rna_seq.augmented_star_gene_counts.tsv$", recursive = T, full.names = TRUE)
  exps <- lapply(allfiles, function(f){
    n <- basename(f)
    exp <- read.delim(f, header  = T, comment.char = "#")
    exp <- exp[grep("protein_coding", exp$gene_type),]
    exp <- exp[match(unique(exp$gene_name), exp$gene_name),]
    r <- data.frame(exp[,"tpm_unstranded"])
    colnames(r) <- sid.data[which(sid.data$file == n),"entity_submitter_id"]
    rownames(r) <- exp[,"gene_name"]
    return(r)
  })                   
  
  Expdata <- do.call(cbind, exps) 
  saveRDS(file = file.path(idir, "../gdc_expression.Rds"), Expdata)
  cat("output file is:", file.path(idir, "../gdc_expression.Rds"), "\n")
  return(Expdata)
}
