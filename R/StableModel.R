#' StableModel
#' select gene that reached the highest C-index to the signature
#'
#'
#' @param Data survival status/time and expression data,column is status,time,g1,g2,...,gn
#' @param covariates vector,genes
#' @param od output dirname
#' @param ... other parameter 
#'
#'
#' @export
#'


StableModel <- function(Data, covariates, seed = NULL, od){
  
  .CoxIndex <- function(Data, covariates, feature = NULL){
    colnames(Data) <- sub("-", "_", colnames(Data))
    covariates <- sub("-", "_", covariates)
    if(is.null(feature)){
      univ_formulas <- sapply(covariates,
                              function(x) as.formula(paste0('Surv(time = time, event = status)~', x)))
    }else{
      univ_formulas <- sapply(covariates,
                              function(x) as.formula(paste0("Surv(time = time, event = status)~", paste(x, feature, sep = "+"))))
    }
    
    univ_models <- lapply( univ_formulas, function(x){coxph(x, data = Data)})
    
    univ_results <- lapply(univ_models,
                           function(x){ 
                             x <- summary(x)        
                             HR <- signif(x$coef[2], digits=2)                      
                             cindex <- signif(x$concordance[1], digits=2)  
                             #beta <- signif(x$coef[1], digits=2)
                             res <- c(cindex, HR)
                             names(res) <- c("C-index", "HR")
                             return(res)
                           })
    res <- t(as.data.frame(univ_results, check.names = FALSE))
    res <- as.data.frame(res)
    mres <- res[which(res$`C-index` == max(res$`C-index`)),]
    mres <- mres[which(mres$HR == max(mres$HR)),]
    return(mres)
  }
  .calculateCindex <- function(data, use_feature){
    PI_data <- RiskScore(data, feature = use_feature, method = "median")
    cond <- concordance.index(PI_data$risk_score, data$time, data$status)
    c.index <- cond$c.index
    #cat(paste(use_feature, collapse = ","), " cindex is:", c.index, "\n")
    return(c.index)
  }
  .AddFeatureSelect <- function(Addgene){
    recovariates <- setdiff(covariates, Addgene)
    mcindex <- NULL
    for(a in recovariates){
      tmp.add <- c(Addgene, a)
      tmp.cindex <- .calculateCindex(Data, tmp.add)
      names(tmp.cindex) <- a
      mcindex <- c(mcindex, tmp.cindex)
    }
    selec <- mcindex[which(mcindex == max(mcindex))]
    return(selec)
  }
  cindex <- 0
  model.gene <- NULL
  step <- 1
  if(is.null(seed)){
    if(step == 1){
      seedM <- .CoxIndex(Data, covariates)
      print(seedM)
      seed.gene <- rownames(seedM)
      cindex <- seedM[,1] %>% as.numeric %>% unique
      model.gene <- c(model.gene, seed.gene)
    }
    cat("seed gene:", seed.gene, ";;", "cindex:", cindex, ";;", "HR:", seedM[,2], "\n")
  }else{
    if(length(seed) == 1){
     seedM <- .CoxIndex(Data, seed)
      print(seedM)
      seed.gene <- rownames(seedM)
      cindex <- seedM[,1] %>% as.numeric %>% unique
      model.gene <- c(model.gene, seed.gene)
    }else{
    model.gene <- seed
    cindex <- .calculateCindex(Data, model.gene)
    }
    cat("seed gene:", model.gene, ";;", "cindex:", cindex, "\n")
  }
  for(s in 1:length(covariates)){
    tcindex <- .AddFeatureSelect(Addgene = model.gene)
    #print(tcindex)
    add.gene <- names(tcindex)
    increVal <- tcindex-cindex
    
    cat("module gene:", model.gene, "current cindex:",cindex, ";;", "new cindex:", tcindex, ";;", "add gene:", add.gene, ";;", "increVal:", increVal, "\n")
    if(increVal < 0.0001){
      break
    }else{
      model.gene <- c(model.gene, add.gene)
      cindex <- tcindex
    }
    model.gene <- c(model.gene, add.gene)
    cindex <- tcindex
    step <- step+1
  }
  return(model.gene)
}
