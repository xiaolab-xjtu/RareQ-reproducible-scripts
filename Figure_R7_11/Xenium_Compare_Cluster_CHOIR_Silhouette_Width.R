setwd('/home/rstudio/Projects/Rare_cell/data/Xenium/Test_cluster/')
# ref: https://www.nature.com/articles/s41587-022-01455-3

library(anndata)
library(Seurat)
library(RareQ)
library(CHOIR)
library(ggplot2)
library(ggpubr)

adata = read_h5ad('../Brain/ms_brain_multisection1.h5ad')
mat <- t(adata$X)
all(dimnames(mat)[[1]][1:248]==adata$var$`Ensembl ID`[1:248])
#dimnames(mat)[[1]] <- adata$var$`Ensembl ID`


######################################################################################################################
### For RareQ analysis and CHOIR-based testing  ######################################################################
######################################################################################################################
# Then, we preprocessed the dataset using the Scanpy pipeline70: normalized the total count of each cell to 1,000,
# log1p transformed the counts and scaled the transformed counts to Z scores (PMID: 38092912)

seu_obj <- CreateSeuratObject(count=mat, project = "sc_object", min.cells = 3)
seu_obj <- NormalizeData(seu_obj, scale.factor = 80) %>% ScaleData()  #
seu_obj <- RunPCA(seu_obj, features = dimnames(seu_obj)[[1]], npcs = 50)
seu_obj <- FindNeighbors(object = seu_obj,
                           k.param = 20,
                           compute.SNN = F,
                           prune.SNN = 0,
                           reduction = "pca",
                           dims = 1:20,
                           force.recalc = F, return.neighbor = T)
seu_obj <- RunUMAP(seu_obj, dims=1:50)
cluster = FindRare(sc_object = seu_obj)
seu_obj$cluster = cluster
seu_obj$X = adata$obs$x_centroid
seu_obj$Y = adata$obs$y_centroid
#saveRDS(seu_obj, file = '../Seu_obj.RDS')





seu_obj <- readRDS('../Seu_obj.RDS')
seu_obj_standard <- seu_obj
seu_obj_standard <- FindNeighbors(object = seu_obj_standard,
                                  k.param = 20,
                                  reduction = "pca",
                                  dims = 1:20)


norm.data = seu_obj@assays$RNA@data

cluster.count <- sort(table(seu_obj$cluster))

## Rename cluster according to the size of cluster
cluster.ordered <- (1:length(cluster.count))[match(seu_obj$cluster, as.integer(names(cluster.count)))]

ave.list <- tapply(1:dim(norm.data)[2], cluster.ordered, function(x){
  rowMeans(norm.data[,x])
})
ave.mat <- (matrix(unlist(ave.list), nrow=length(ave.list), ncol=dim(norm.data)[1], byrow = T))
dimnames(ave.mat) <- list(names(ave.list), rownames(norm.data))


seu_obj$cluster_ordered = cluster.ordered
Idents(seu_obj) <- 'cluster_ordered'

dist.mat <- as.matrix(dist(ave.mat))

test.res <- data.frame(result=NA, mean_accuracy=NA, var_accuracy=NA, mean_permuted_accuracy=NA, var_permuted_accuracy=NA)
for(kk in 1:dim(dist.mat)[1]){
  ident1 = rownames(dist.mat)[kk]
  ident2 = as.character(order(dist.mat[kk,])[2])
  test <- compareClusters(seu_obj, ident1 = ident1, ident2 = ident2, group_by = 'cluster_ordered', feature_set='all', nn_matrix = seu_obj_standard@graphs$RNA_nn)
  print(paste0(kk, ':', test$comparison_result))
  test.res[kk,'result'] <- test$comparison_result
  test.res[kk,2:5] <- c(test$comparison_records$mean_accuracy,
                        test$comparison_records$var_accuracy,
                        test$comparison_records$mean_permuted_accuracy,
                        test$comparison_records$var_permuted_accuracy)
}

#saveRDS(test.res, file = 'Test_result.RDS')

#test.res <- readRDS('Test_result.RDS')

df <- rbind(data.frame(cluster=1:dim(test.res)[1], mean_accuracy=test.res$mean_accuracy, var_accuracy=test.res$var_accuracy, Group=rep('Prediction accuracy', dim(test.res)[1])),
            data.frame(cluster=1:dim(test.res)[1], mean_accuracy=test.res$mean_permuted_accuracy, var_accuracy=test.res$var_permuted_accuracy, Group=rep('Permuted accuracy', dim(test.res)[1])))
df$std_accuracy = sqrt(df$var_accuracy)
df$Group <- factor(df$Group, levels = c('Prediction accuracy', 'Permuted accuracy'), ordered = T)


p.test <- ggplot(df, aes(x = cluster, y = mean_accuracy, group = Group, color=factor(Group))) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +

  geom_errorbar(
    aes(ymin = mean_accuracy - std_accuracy, ymax = mean_accuracy + std_accuracy),
    width = 0.1,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c('#059E75', '#F0426B')) + theme_bw() +
  theme(panel.border = element_rect(fill = NA), axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.title = element_blank(), legend.position = 'top',
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
  ) + scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x='Cluster (with cluster size small to large)', y='Accuracy (mean with standard deviation)')

# ggsave(p.test, filename = 'CHOIR_Cluster_Test_plot.pdf', width = 180, height = 90, units = 'mm')


df.ind <- data.frame(cluster=1:dim(test.res)[1],
                     ind=test.res$result)

df.ind$res <- ifelse(df.ind$ind=='split', 'Valid', 'Invalid')
df.ind$res <- factor(df.ind$res, levels=c('Valid', 'Invalid'), ordered = T)

p.test.ind <- ggplot(df.ind, aes(x = cluster, y = res, group = factor(res), color=factor(res))) + geom_point(size=0.5) +
  scale_color_manual(values = c('black', 'grey')) + theme_bw() +
  theme(panel.border = element_rect(fill = NA), axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.title = element_blank(), legend.position = 'none',
        panel.grid = element_blank()
  ) + scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x='Cluster (with cluster size small to large)', y=NULL)



p.comb = ggarrange(p.test, p.test.ind, ncol=1, align = 'hv', heights = c(1, 0.3))

# ggsave(p.comb, filename = 'CHOIR_Cluster_Test_plot_comb.pdf', width = 180, height = 100, units = 'mm')



df$Test = df.ind$ind
df$Result = df.ind$res
# write.csv(df, file='Xenium_CHOIR_test.csv')




## Heatmaps of DEGs between small clusters with their nearest clusters
norm.count = seu_obj@assays$RNA@scale.data
p.comb = do.call(ggarrange,lapply(1:30, function(kk){
  ident1 = rownames(dist.mat)[kk]
  ident2 = as.character(order(dist.mat[kk,])[2])
  mk.12 = FindMarkers(seu_obj, ident.1 = ident1, ident.2 = ident2, min.pct=0, logfc.threshold = -Inf)
  features = rownames(mk.12)[abs(mk.12$avg_log2FC) > 0.15]
  if(length(features) < 1){
    features = rownames(mk.12)[abs(mk.12$avg_log2FC) > 0.05]
  }
  # Subsampling cells in ident2 when the cells are too many relative the cells in ident1

  idx.ident1 = which(seu_obj$cluster_ordered==ident1)
  idx.ident2 = which(seu_obj$cluster_ordered==ident2)
  if(length(idx.ident2)/length(idx.ident1) >= 5){
    set.seed(2025)
    idx.ident2.sample = sample(idx.ident2, size=5*length(idx.ident1))
  }else{
    idx.ident2.sample = idx.ident2
  }

  lab1 = paste0('C',ident1)
  lab2 = paste0('C',ident2)
  cols = c('#FF0099', '#C3EF00')
  names(cols) = c(lab1, lab2)
  ha1 = HeatmapAnnotation('Group' = factor(c(rep(lab1, length(idx.ident1)),
                                             rep(lab2, length(idx.ident2.sample))),
                                           levels=c(lab1,lab2), ordered = T),
                          col=list('Group' = cols))

  ht = Heatmap(cbind(norm.count[features,idx.ident1], norm.count[features, idx.ident2.sample]),
               top_annotation = ha1, show_heatmap_legend = F,
               cluster_columns = F, show_row_names = F, show_column_names = F, show_row_dend = F)
  plt = ggplotify::as.ggplot(ht)
  return(plt)
  # plt.list[kk] = plt


}))

ggsave(p.comb, filename = 'Xenium_cluster_neighbor_comparison_heatmap.pdf', width = 10, height = 12)







######################################################################################################################
### For Seurat analysis and CHOIR-based testing  ######################################################################
######################################################################################################################
## Seurat benchmarking
resolutions = c(11, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12)
for(reso in resolutions){
  seu_obj_standard <- FindClusters(seu_obj_standard, resolution = reso)
}


#saveRDS(seu_obj_standard, file = 'Seu_obj_standard_for_testing.RDS')
seu_obj_standard <- readRDS('Seu_obj_standard_for_testing.RDS')


cluster.num = c()
for(reso in  resolutions){
  cluster.num = c(cluster.num, length(unique(seu_obj_standard@meta.data[,paste0('RNA_snn_res.', as.character(reso))])))
}

which(cluster.num==length(unique(seu_obj$cluster_ordered)))
# 3 6
reso.select = paste0('RNA_snn_res.', as.character(resolutions[3]))

norm.data = seu_obj_standard@assays$RNA@data

cluster.count <- sort(table(seu_obj_standard@meta.data[,reso.select]))

## Rename cluster according to the size of cluster
cluster.ordered <- (1:length(cluster.count))[match(seu_obj_standard@meta.data[,reso.select], as.integer(names(cluster.count)))]

ave.list <- tapply(1:dim(norm.data)[2], cluster.ordered, function(x){
  rowMeans(norm.data[,x])
})
ave.mat <- (matrix(unlist(ave.list), nrow=length(ave.list), ncol=dim(norm.data)[1], byrow = T))
dimnames(ave.mat) <- list(names(ave.list), rownames(norm.data))


seu_obj_standard$cluster_ordered = cluster.ordered
Idents(seu_obj_standard) <- 'cluster_ordered'


dist.mat <- as.matrix(dist(ave.mat))

test.res <- data.frame(result=NA, mean_accuracy=NA, var_accuracy=NA, mean_permuted_accuracy=NA, var_permuted_accuracy=NA)
for(kk in 1:dim(dist.mat)[1]){
  ident1 = rownames(dist.mat)[kk]
  ident2 = as.character(order(dist.mat[kk,])[2])
  test <- compareClusters(seu_obj_standard, ident1 = ident1, ident2 = ident2, group_by = 'cluster_ordered', feature_set='all', nn_matrix = seu_obj_standard@graphs$RNA_nn)
  print(paste0(kk, ':', test$comparison_result))
  test.res[kk,'result'] <- test$comparison_result
  test.res[kk,2:5] <- c(test$comparison_records$mean_accuracy,
                        test$comparison_records$var_accuracy,
                        test$comparison_records$mean_permuted_accuracy,
                        test$comparison_records$var_permuted_accuracy)
}

#saveRDS(test.res, file = 'Test_result_Seurat.RDS')

#test.res <- readRDS('Test_result_Seurat.RDS')

df <- rbind(data.frame(cluster=1:dim(test.res)[1], mean_accuracy=test.res$mean_accuracy, var_accuracy=test.res$var_accuracy, Group=rep('Prediction accuracy', dim(test.res)[1])),
            data.frame(cluster=1:dim(test.res)[1], mean_accuracy=test.res$mean_permuted_accuracy, var_accuracy=test.res$var_permuted_accuracy, Group=rep('Permuted accuracy', dim(test.res)[1])))
df$std_accuracy = sqrt(df$var_accuracy)
df$Group <- factor(df$Group, levels = c('Prediction accuracy', 'Permuted accuracy'), ordered = T)


p.test <- ggplot(df, aes(x = cluster, y = mean_accuracy, group = Group, color=factor(Group))) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +

  geom_errorbar(
    aes(ymin = mean_accuracy - std_accuracy, ymax = mean_accuracy + std_accuracy),
    width = 0.1,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c('#059E75', '#F0426B')) + theme_bw() +
  theme(panel.border = element_rect(fill = NA), axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.title = element_blank(), legend.position = 'top',
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
  ) + scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x='Cluster (with cluster size small to large)', y='Accuracy (mean with standard deviation)')

# ggsave(p.test, filename = 'CHOIR_Cluster_Test_plot.pdf', width = 180, height = 90, units = 'mm')


df.ind <- data.frame(cluster=1:dim(test.res)[1],
                     ind=test.res$result)

df.ind$res <- ifelse(df.ind$ind=='split', 'Valid', 'Invalid')
df.ind$res <- factor(df.ind$res, levels=c('Valid', 'Invalid'), ordered = T)

p.test.ind <- ggplot(df.ind, aes(x = cluster, y = res, group = factor(res), color=factor(res))) + geom_point(size=0.5) +
  scale_color_manual(values = c('black', 'grey')) + theme_bw() +
  theme(panel.border = element_rect(fill = NA), axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), legend.title = element_blank(), legend.position = 'none',
        panel.grid = element_blank()
  ) + scale_x_continuous(expand = c(0.02,0.02)) +
  labs(x='Cluster (with cluster size small to large)', y=NULL)



p.comb = ggarrange(p.test, p.test.ind, ncol=1, align = 'hv', heights = c(1, 0.3))

ggsave(p.comb, filename = 'CHOIR_Cluster_Test_plot_comb_Seurat.pdf', width = 180, height = 100, units = 'mm')






## Test cluster using silhouette coefficient

library(Rcpp)

# Compile Rcpp file(make sure silhouette_rcpp.cpp in current working directory)
sourceCpp("silhouette_rcpp.cpp")

silhouette_fast <- function(clusters, data, block_size = 1000) {
  # 1. as matrix
  if (!is.matrix(data)) {
    data <- as.matrix(data[, seq_len(ncol(data)), drop = FALSE])
    message("Data has been transformed into matrix")
  }

  # 2. Check input
  n <- nrow(data)
  if (length(clusters) != n) {
    stop("clusters length (", length(clusters), ") and sample size (", n, ") do not match")
  }
  if (block_size < 100) {
    warning("block_size should be ≥100 to boost efficiency, current value is", block_size)
  }

  # 3. Transform cluster labels as integers
  clusters_int <- as.integer(clusters)

  # 4. Use Rcpp function
  sil_width <- silhouette_rcpp_final(
    data = data,
    clusters = clusters_int,
    block_size = block_size
  )

  return(sil_width)
}





library(dplyr)

embed = seu_obj@reductions$pca@cell.embeddings
dist_matrix <- (dist(embed, method = "euclidean"))


sil_width_RareQ <- silhouette_fast(
  clusters = as.integer(seu_obj$cluster_ordered),
  data = embed
)

sil_width_Seurat <- silhouette_fast(
  clusters = as.integer(seu_obj_standard$cluster_ordered),
  data = embed
)


# saveRDS(sil_width_RareQ, file = 'sil_width_RareQ.RDS')
# saveRDS(sil_width_Seurat, file = 'sil_width_Seurat.RDS')

sil_width_RareQ <- readRDS('sil_width_RareQ.RDS')
sil_width_Seurat <- readRDS('sil_width_Seurat.RDS')


sil.df <- data.frame(sil_width_RareQ = sil_width_RareQ,
                     cluster_RareQ = seu_obj$cluster_ordered,
                     sil_width_Seurat = sil_width_Seurat,
                     cluster_Seurat = seu_obj_standard$cluster_ordered)

p.sil.RareQ <- ggplot(data=sil.df, aes(x=factor(cluster_RareQ), y=sil_width_RareQ)) + geom_boxplot(outlier.shape = NA, color='grey50', width=0.2) +
  stat_summary(fun.y = median, color='red', size=0.1, shape=15) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(yintercept = 0) +
  labs(x='RareQ cluster (with cluster size small to large)', y='Silhouette width') + scale_x_discrete(breaks = c(1, 20, 40,60, 80,100,120,140,160,180,200))
p.sil.RareQ

p.sil.Seurat <- ggplot(data=sil.df, aes(x=factor(cluster_Seurat), y=sil_width_Seurat)) + geom_boxplot(outlier.shape = NA, color='grey50', width=0.2) +
  stat_summary(fun.y = median, color='red', size=0.1, shape=15) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(yintercept = 0) +
  labs(x='Seurat cluster (with cluster size small to large)', y='Silhouette width') + scale_x_discrete(breaks = c(1, 20, 40,60, 80,100,120,140,160,180,200))
p.sil.Seurat

p.sil.comb = ggarrange(p.sil.RareQ,
          p.sil.Seurat, ncol=1)
ggsave(p.sil.comb, filename = 'Xenium_Sil_boxplot_comb.pdf', width = 8, height = 5)









