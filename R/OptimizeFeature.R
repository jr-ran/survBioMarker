#' final feature
#' use univariate COX and LASSO regression optimize feature
#'
#' Random Perturbation
#' description: random n(1000) times OptimizeFeature, return gene set and gene weight
#'
#'
#' @param data  clinical and expression data, fist column is OS status; second column is OS time; feature expression...
#' @param pvalcut univariate cox pvalue cutoff
#' default pvalcut is 0.05
#' @param od output dirname
#' @param n random times
#' default 1000
#' @param ... other parameter
#'
#' @rdname OptimizeFeature
#' @return data.frame genes and weight
#' @export

OptimizeFeature <- function(data, pvalcut = 0.05, od, n = 1000, ...){
  Status <- data[,"status"]
  names(Status) <- rownames(data)
  unicox.re <- uniCox(data, method = "univariate")
  cat(file.path(od, "s1_unixcox_nofilter.xls"), "\n")
  write.table(file = file.path(od, "s1_unixcox_nofilter.xls"), unicox.re, quote = FALSE, sep = "\t", row.names = FALSE)
  #unicox.re <- unicox.re[which(unicox.re$pvalue < pvalcut | unicox.re$HR > 1.02 | unicox.re$HR < 1.02),]
  unicox.re <- unicox.re[which(unicox.re$pvalue < pvalcut),]
  write.table(file = file.path(od, "s1_unixcox_filter.xls"), unicox.re, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("S1 candidate gene number is:", nrow(unicox.re), "\n")
  data1 <- data[,rownames(unicox.re)]
  tlen <- nrow(data1)
  feature.set <- sapply(1:n, function(i){
    nfolds <- 5
    foldid <- sample(rep(seq(nfolds), length = tlen))
    testSet <- which(foldid != 1)
    x <- data1[testSet,] %>% as.matrix()
    y <- Status[testSet]
    #feature <- LfitFeature(x, y, od, time = i)
    feature <- LfitFeature(x, y, od = NULL, time = i)
    feature
  })
  gset <- unlist(feature.set) %>% table
  gweight <- gset/n %>% as.numeric()
  gs <- sort(gweight, decreasing = T) %>% data.frame()
  colnames(gs) <- c("genes", "weight")
  write.table(file = file.path(od, "s2_optimize_nofilter.xls"), gs, quote = FALSE, sep = "\t", row.names = FALSE)
  fgs <- gs %>% filter(weight > 0.01)
  cat("S2 optimize gene number is:", nrow(fgs), "\n")
  return(fgs)
}

