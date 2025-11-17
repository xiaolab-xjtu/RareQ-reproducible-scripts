library(anndata)
library(Seurat)
library(RareQ)
library(magrittr)


## Analysis real data


setwd('/home/rstudio/Projects/Rare_cell/data/')
scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                "UUOkidney", "Zelsel", "Choroid_plexus")

load("secure_seeds.RData")



for(sc in scenarios){
  path <- paste0(sc,'/')
  print(sc)
  GEX <- read_mtx(paste0(path, 'Gene_Cell.mtx'))
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
    # sc_object <- RunUMAP(sc_object, dims=1:50)
    cluster = FindRare(sc_object = sc_object)
    saveRDS(cluster, file=paste0(path, 'OUR_result_shuffled_cols_seed_',k,'.RDS'))

    rm(list=c('sc_object', 'cluster'))
    gc()
  }
}








## Extract top features of the rare cell types
scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                "UUOkidney", "Zelsel", "Choroid_plexus")


for(sc in scenarios){

  path <- paste0(sc,'/')
  if(file.exists(paste0(path, 'Top_marker_ref.RDS'))){next}

  GEX <- read_mtx(paste0(path, 'Gene_Cell.mtx'))
  GEX.mat <- GEX$X

  labs <- readRDS(paste0(path, 'label.RDS'))
  type.cnt <- sort(table(labs))

  rare_types <- names(type.cnt)[type.cnt <= sum(type.cnt) * 0.01]
  if(length(rare_types) < 1){
    rare_types <- names(type.cnt)[which.min(type.cnt)]
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




