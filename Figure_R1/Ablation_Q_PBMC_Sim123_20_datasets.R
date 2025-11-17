library(magrittr)
library(anndata)
library(Seurat)
library(RareQ)
library(random)

setwd('/home/rstudio/Projects/Rare_cell/data/PBMCs/')



source('FindRare_Q_Ablation.R', encoding='utf-8')



scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10")
repeats <- c(0:49)


for(sc in scenarios){
  for(k in repeats){

    path <- paste0('site1/',sc,'/',k,'/R_input/')
    GEX <- read_mtx(paste0(path, 'sub_gene_cell.mtx'))
    sc_object <- CreateSeuratObject(count=GEX$X, project = "sc_object", min.cells = 3)
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
    #sc_object <- RunUMAP(sc_object, dims=1:50)

    cluster = FindRare_Ablation_Q(sc_object = sc_object)
    saveRDS(cluster, file=paste0(path, 'OUR_result_Ablation_Q_1.RDS'))

  }
}





## Analysis real data
library(anndata)

setwd('/home/rstudio/Projects/Rare_cell/data/')
scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                "UUOkidney", "Zelsel", "Choroid_plexus")


for(sc in scenarios){
  path <- paste0(sc,'/')
  print(sc)
  GEX <- read_mtx(paste0(path, 'Gene_Cell.mtx'))
  sc_object <- CreateSeuratObject(count=GEX$X, project = "sc_object", min.cells = 3)
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
  saveRDS(cluster, file=paste0(path, 'OUR_result_Ablation_Q_1.RDS'))

  rm(list=c('sc_object', 'GEX', 'cluster'))
  gc()

}



