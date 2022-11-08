#' MainRun
#' find signature and survival analysis 
#'
#'
#' @param initiData survival status/time and expression data,column is status,time,g1,g2,...,gn
#' @param geneS vector,genes
#' @param out output dirname
#' @param pvalcut unicox model pvalue, filter gene, default 0.05
#' @param key OS or DFS
#' default OS
#' @param phenoData phelotype dataset, rowname is sample
#' @param var.pheno vector phelotype name,example:age,sex,AJCC
#' @param ... other parameter 
#'
#'
#' @export
#'


MainRun <- function(initiData, geneS, out, pvalcut = 0.05, key = "OS", phenoData = NULL, var.pheno = NULL,...){
  # ------- module biomarker                       ------------------------------
  # feature choose
  print("optimize gene.")
  inputData <- initiData[,colnames(initiData) %in% c("status", "time", geneS)]
  if(file.exists(file.path(out, "optimize_gene_set.txt"))){
    optimize_gene_set <- read.delim(file.path(out, "optimize_gene_set.txt"))
  }else{
    optimize_gene_set <- OptimizeFeature(inputData, pvalcut = pvalcut, od = out)
    write.table(file = file.path(out, "optimize_gene_set.txt"), optimize_gene_set, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  print("Choose prognostic markers.")
  weightg <- optimize_gene_set %>% filter(weight >= 0.95) %>% .[, 'genes'] %>% as.character
  traing <- optimize_gene_set %>% filter(weight < 0.95) %>% .[, 'genes'] %>% as.character
  if(nrow(optimize_gene_set) < 1){
    cat("Warnning::optimize gene set is zero, please adjust the parameters of pvalue!!!")
    return(NULL)
  }
  seed <- NULL
  if(length(weightg) >= 1){
    seed <- weightg
  }
  if(length(traing) >= 1){
    cat("weightg is:", head(weightg), "\n")
    cat("traing is:", head(traing), "\n")
    choose_feature <- StableModel(inputData, covariates = traing, seed = seed, od = out)
  }else{
    choose_feature <- seed
  }
  coxFeature <- uniCox(inputData[,colnames(inputData) %in% c("status", "time", choose_feature)], method = "univariate")
  write.table(file = file.path(out, "final_genes_coff.xls"), coxFeature, quote = FALSE, sep = "\t", row.names = FALSE)

  # ------- survival analysis ------------------
  ####
  print("Calculate risk score.")
  sdata <- inputData[,colnames(inputData) %in% c("status", "time", choose_feature)]
  PI_data <- RiskScore(sdata, feature = choose_feature, method = "median")
  finalData <- cbind(sdata[,1:2], PI_data)
  var <- "risk_group"
  print("Survival analysis.")
  KMSurvival(finalData, var, out, key)
  
  # ------ model validate -------------------------
  # other phelotype 
  if(!is.null(var.pheno)){
    if(is.null(phenoData)){
      print("Error:: please set parameter phenoData!!!")
      return(NA)
    }
    if(length(intersect(var.pheno, colnames(phenoData))) == 0){
      print("Error:: please check colnames of parameter phenoData!!!")
      return(NA)
    }
    if(length(intersect(rownames(finalData), rownames(phenoData))) == 0 || length(intersect(rownames(finalData), rownames(phenoData))) < nrow(finalData)){
      print("Error:: please check the samples of phenoData!!!")
      return(NA)
    }
    var.name <- intersect(var.pheno, colnames(phenoData))
    finalData <- cbind(finalData, phenoData[rownames(finalData),var.name])
    var <- c(var, var.name)
  }
  # model access
  print("univariate cox model.")
  coxData <- cox_data_process(finalData, var.names = var)
  uniValid <- uniCox(coxData, var = var, method = "univariate")
  uniValid$coef <- NULL
  write.table(file = file.path(out, "univariate_HR.xls"), uniValid, quote = FALSE, sep = "\t", row.names = FALSE)

  if(length(var) > 1){
    print("multivariate cox model.")
    multiValid <- uniCox(coxData, var = var, method = "multivariate")
    write.table(file = file.path(out, "multivariate_HR.txt"), multiValid, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  # ROC
  print("Calculate model AUC.")
  rocObj <- auc(finalData$status, finalData$risk_score, out)
  return("Done")
}