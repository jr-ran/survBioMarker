#' gdcClinOS
#' get overall survival information form GDC
#'
#'
#' @param mfile GDC clinical.cases_selection file
#' @param mjson GDC metadata.cart json
#' @param col.n select mfile colnames ,defaut:c("entity_submitter_id", "days_to_death", "days_to_last_follow_up", "vital_status")
#' @param ... other parameter 
#'
#' 
#' @export
#'

gdcClinOS <- function(mfile, mjson, col.n = c("entity_submitter_id", "days_to_death", "days_to_last_follow_up", "vital_status")){
  
  sid.data <- gdcSampleID(mjson)
  # clinical
  od <- dirname(mfile)
  clinical <- read.delim(mfile)
  clin <- inner_join(sid.data[,1:2], clinical, by = "case_id")
  clinD <- clin[,col.n]
  clinD <- clinD[match(unique(clinD$entity_submitter_id), clinD$entity_submitter_id),]
  # trans status format
  alive.index <- which(clinD$vital_status == "Alive")
  dead.index <- which(clinD$vital_status == "Dead")
  status <- vector()
  status[alive.index] <- 0
  status[dead.index] <- 1
  # get OS time
  os.time <- vector()
  os.time[alive.index] <- clinD$days_to_last_follow_up[alive.index]
  os.time[dead.index] <- clinD$days_to_death[dead.index]
  os.data <- data.frame(sample = clinD$entity_submitter_id, status = status, OS = os.time)
  write.table(file = file.path(od, "clinical_os.txt"), os.data, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("output file is:", file.path(od, "clinical_os.txt"), "\n")
  return(os.data)
}
