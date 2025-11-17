## Sensitivity on Jurkat simulated data with Q Ablation
library(anndata)
library(Seurat)
library(RareQ)
library(magrittr)
setwd('/home/rstudio/Projects/Rare_cell/data/Jurkat/Sensitivity/')



source('FindRare_Q_Ablation.R', encoding='utf-8')



gene.sample.mat <- readRDS('sample_matrix.RDS')
gene.group <- readRDS('Gene_group.RDS')
nonDE.idx <- which(gene.group == 'NonDE')

GEX <- read_mtx('Gene_Cell.mtx')
count.all <- GEX$X

for(k in 1:dim(gene.sample.mat)[1]){

  ind.k <- which(gene.sample.mat[k,]==1)
  count.k <- count.all[c(ind.k, nonDE.idx),]

  sc_object <- CreateSeuratObject(count=count.k, project = "sc_object", min.cells = 3)
  sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
  sc_object <- subset(sc_object, percent.mt<20)
  sc_object <- NormalizeData(sc_object) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData()
  sc_object <- RunPCA(sc_object, features = VariableFeatures(object = sc_object))
  sc_object <- FindNeighbors(object = sc_object,
                             k.param = 20,
                             compute.SNN = F,
                             prune.SNN = 0,
                             reduction = "pca",
                             dims = 1:50,
                             force.recalc = F, return.neighbor = T)
  # sc_object <- RunUMAP(sc_object, dims=1:50)
  cluster = FindRare_Ablation_Q(sc_object = sc_object)
  saveRDS(cluster, file=paste0('OUR/DE_', k, '_Q_Ablation.RDS'))

  rm(list=c('sc_object', 'GEX', 'cluster', 'ind.k', 'count.k'))
  gc()

}

