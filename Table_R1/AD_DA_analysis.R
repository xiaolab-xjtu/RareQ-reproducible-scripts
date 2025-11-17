setwd('/home/rstudio/Projects/Rare_cell/data/GSE303823_AD_DLB_PDD/')

library(RareQ)
library(miloR)
library(ggplot2)

AD.obj <- readRDS('AD_obj.RDS')

AD.obj$condition = as.character(AD.obj$group)
AD.obj$condition <- factor(AD.obj$condition, levels=c('CTRL','ADD'), ordered = T)


run_nbglm <- function(
    obj,
    condition_col,
    sample_col,
    batch_col=NULL,
    norm.method="TMM"
){
  sce <- as.SingleCellExperiment(obj)
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)])
  design_df <- distinct(design_df)
  #design_df <- dplyr::rename(sample=sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }

  clust.ids <- colData(sce)[['cluster_sort']]

  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  clust.df <- data.frame("cell_id"=colnames(sce), "Louvain.Clust"=as.character(clust.ids))
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec

  louvain.count <- table(clust.df$Louvain.Clust, clust.df$Sample)
  attributes(louvain.count)$class <- "matrix"

  ## Test with same NB-GLM model as the Milo
  if(norm.method %in% c("TMM")){
    message("Using TMM normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
    dge <- calcNormFactors(dge, method="TMM")
  } else if(norm.method %in% c("logMS")){
    message("Using logMS normalisation")
    dge <- DGEList(counts=louvain.count,
                   lib.size=colSums(louvain.count))
  }

  model <- model.matrix(design, data=design_df)
  rownames(model) <- design_df$sample
  model <- model[colnames(louvain.count), ]

  dge <- estimateDisp(dge, model)
  fit <- glmQLFit(dge, model, robust=TRUE)
  n.coef <- ncol(model)
  louvain.res <- as.data.frame(topTags(glmQLFTest(fit, coef=n.coef), sort.by='none', n=Inf))

  clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, 'logFC']
  clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, 'FDR']
  return(list(sc=clust.df, mc=louvain.res))
}




# DA.res <- run_nbglm(obj = AD.obj, sample_col = 'sample', condition_col = 'condition')
# saveRDS(DA.res, file = 'AD_DA_res.RDS')

DA.res <- readRDS('AD_DA_res.RDS')

plot.data = DA.res$mc
plot.data$cluster = rownames(plot.data)
plot.data$cluster = factor(plot.data$cluster, levels=as.character(1:dim(plot.data)[1]), ordered = T)

p.DA.plot <- ggplot(data=plot.data) + geom_point(aes(y=logFC, x=cluster, size=-log10(PValue)))


plot.data.sort = plot.data[order(as.integer(plot.data$cluster)),]
write.csv(plot.data.sort, file = 'AD_DA_res.csv')



