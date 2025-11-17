library(anndata)
library(Seurat)
library(RareQ)
library(magrittr)
library(random)

setwd('/home/rstudio/Projects/Rare_cell/data/PBMCs/')


## generate seeds
#secure_seeds <- randomNumbers(n = 30, min = 1, max = 1e6, col = 1)
#table(duplicated(secure_seeds[,1]))
#FALSE
#30

# Save it
# save(secure_seeds, file = "../secure_seeds.RData")  # 务必保存

load("../secure_seeds.RData")
#set.seed(secure_seeds[1])


scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10")
repeats <- c(0:49)

for(sc in scenarios){
  for(rep in repeats){
    path <- paste0('site1/',sc,'/',rep,'/R_input/')
    GEX <- read_mtx(paste0(path, 'sub_gene_cell.mtx'))
    GEX.mat <- GEX$X
    for(k in 1:length(secure_seeds)){
      set.seed(secure_seeds[k])
      shuffled_cols <- sample(ncol(GEX.mat))
      saveRDS(shuffled_cols, file = paste0(path, 'shuffled_cols_seed_',k,'.RDS'))

      sc_object <- CreateSeuratObject(count=GEX.mat[,shuffled_cols], project = "sc_object", min.cells = 3)
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
      cluster = FindRare(sc_object = sc_object)
      saveRDS(cluster, file=paste0(path, 'OUR_result_shuffled_cols_seed_',k,'.RDS'))

    }
  }
}








## Extract top features of the rare cell types
scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10")
repeats <- c(0:49)

for(sc in scenarios){
  for(rep in repeats){

    path <- paste0('site1/',sc,'/',rep,'/R_input/')
    if(file.exists(paste0(path, 'Top_marker_ref.RDS'))){next}

    GEX <- read_mtx(paste0(path, 'sub_gene_cell.mtx'))
    GEX.mat <- GEX$X

    lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
    labs <- lab.df$V1
    type.cnt <- table(labs)

    if(sc %in% c('Rare1_Ordinary4', 'Rare1_Ordinary9')){
      rare_types <- names(type.cnt)[which.min(type.cnt)]
    }else{
      rare_types <- names(type.cnt)[order(type.cnt)[1:5]]
    }

    sc_object <- CreateSeuratObject(count=GEX.mat, project = "sc_object", min.cells = 3)
    sc_object$percent.mt <- PercentageFeatureSet(sc_object, pattern = "^MT-")
    sc_object <- subset(sc_object, percent.mt<20)
    sc_object <- NormalizeData(sc_object) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData()
    sc_object$cell_type = labs
    Idents(sc_object) = 'cell_type'

    type.vec = c()
    marker.vec = c()
    for(rare_type in rare_types){
      mk = FindMarkers(sc_object, ident.1=rare_type, only.pos=T, logfc.threshold=0.01)
      type.id <- which(rare_types == rare_type)
      type.vec = c(type.vec, rep(paste0('R',type.id), 50))
      marker.vec = c(marker.vec, rownames(mk)[1:50])
    }
    marker.df = data.frame(type = type.vec, marker=marker.vec)
    saveRDS(marker.df, file=paste0(path, 'Top_marker_ref.RDS'))

    for(k in 1:length(secure_seeds)){

      shuffled_cols <- readRDS(paste0(path, 'shuffled_cols_seed_',k,'.RDS'))
      shuffled_res <- readRDS(paste0(path, 'OUR_result_shuffled_cols_seed_',k,'.RDS'))
      ordered_res = shuffled_res[order(shuffled_cols)]

      type.vec.pred = c()
      marker.vec.pred = c()
      for(rare_type in rare_types){
        rare.cluster <- tapply(labs, ordered_res, function(x){
          type.cnt <- table(x)
          if(names(type.cnt)[which.max(type.cnt)]==rare_type){
            return(T)
          }else{
            return(F)
          }
        })

        type.id <- which(rare_types == rare_type)
        pred.lab <- ifelse(ordered_res %in% names(rare.cluster)[rare.cluster], 1, 0)
        if(any(pred.lab==1)){
          sc_object$prediction = pred.lab
          Idents(sc_object) = 'prediction'

          mk.pred = FindMarkers(sc_object, ident.1=1, only.pos=T, logfc.threshold=0.01)
          type.vec.pred = c(type.vec.pred, rep(paste0('R',type.id), 50))
          marker.vec.pred = c(marker.vec.pred, rownames(mk.pred)[1:50])
        }else{
          type.vec.pred = c(type.vec.pred, rep(paste0('R',type.id), 50))
          marker.vec.pred = c(marker.vec.pred, rep(NA, 50))
        }

      }
      marker.df.pred = data.frame(type = type.vec.pred, marker=marker.vec.pred)
      saveRDS(marker.df.pred, file=paste0(path, 'Top_marker_shuffled_cols_seed_',k,'.RDS'))

    }
  }
}







