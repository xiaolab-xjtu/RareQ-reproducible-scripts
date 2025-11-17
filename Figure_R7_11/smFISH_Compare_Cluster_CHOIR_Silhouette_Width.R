setwd('/home/rstudio/Projects/Rare_cell/data/smFish/')


library(anndata)
library(Seurat)
library(RareQ)
library(CHOIR)
library(ggplot2)
library(ggpubr)

seu <- readRDS('Seu_obj.RDS')
seu_obj <- readRDS('Seu_obj_resolutions_Standard_Pipeline.RDS')
norm.data = seu_obj@assays$RNA@data
seu_obj$cluster = seu$Cluster

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
  test <- compareClusters(seu_obj, ident1 = ident1, ident2 = ident2, group_by = 'cluster_ordered', feature_set='all', nn_matrix = seu_obj@graphs$RNA_nn)
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

ggsave(p.test, filename = 'CHOIR_Cluster_Test_plot.pdf', width = 180, height = 90, units = 'mm')


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



p.test.publication = ggarrange(p.test, p.test.ind, ncol=1, align = 'hv', heights = c(1, 0.3))




######################################################################################################################
### For Seurat analysis and CHOIR-based testing  ######################################################################
######################################################################################################################
## Seurat benchmarking
resolutions = c(17.8, 17.9, 18.1, 18.3, 18.5, 18.7, 18.9)
for(reso in resolutions){
  seu_obj <- FindClusters(seu_obj, resolution = reso)
}


# saveRDS(seu_obj, file = 'Seu_obj_for_testing_Seurat.RDS')
seu_obj <- readRDS('Seu_obj_for_testing_Seurat.RDS')

cluster.num = c()
for(reso in  resolutions){
  cluster.num = c(cluster.num, length(unique(seu_obj@meta.data[,paste0('RNA_snn_res.', as.character(reso))])))
}

which(cluster.num==length(unique(seu$Cluster)))
# 3
reso.select = paste0('RNA_snn_res.', as.character(resolutions[3]))

norm.data = seu_obj@assays$RNA@data

cluster.count <- sort(table(seu_obj@meta.data[,reso.select]))

## Rename cluster according to the size of cluster
cluster.ordered <- (1:length(cluster.count))[match(seu_obj@meta.data[,reso.select], as.integer(names(cluster.count)))]

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
  test <- compareClusters(seu_obj, ident1 = ident1, ident2 = ident2, group_by = 'cluster_ordered', feature_set='all', nn_matrix = seu_obj@graphs$RNA_nn)
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








## side by side reconciliation with the original labels
## Rename cluster according to the size of cluster
cluster.count <- sort(table(seu$cluster))
cluster.ordered <- (1:length(cluster.count))[match(seu$cluster, as.integer(names(cluster.count)))]
seu$cluster_ordered = cluster.ordered

## Rename cluster according to the size of cluster
Cluster.count <- sort(table(seu$Cluster))
Cluster.ordered <- (1:length(Cluster.count))[match(seu$Cluster, as.integer(names(Cluster.count)))]
seu$Cluster_ordered = Cluster.ordered


uniq.cluster.RareQ = sort(unique(seu$cluster_ordered))
uniq.cluster.orig = sort(unique(seu$Cluster_ordered))


jaccard.mat <- matrix(0, nrow=length(uniq.cluster.RareQ), ncol=length(uniq.cluster.orig))
for(m in 1:dim(jaccard.mat)[1]){
  for(n in 1:dim(jaccard.mat)[2]){
    id.m = which(seu$cluster_ordered==uniq.cluster.RareQ[m])
    id.n = which(seu$Cluster_ordered==uniq.cluster.orig[n])

    jaccard.sim <- length(intersect(id.m, id.n))/length(union(id.m, id.n))
    jaccard.mat[m,n] <- jaccard.sim
  }
}
dimnames(jaccard.mat) <- list(uniq.cluster.RareQ, uniq.cluster.orig)
write.csv(jaccard.mat, file='Cluster_jaccard_matrix.csv')
col_fun = circlize::colorRamp2(c(0,1), c("#FFF5F0",'#A50F15'))
#col_fun = circlize::colorRamp2(c(0,1), c("grey95", "#ff349c"))
p.jaccard <- as.ggplot(ComplexHeatmap::Heatmap(jaccard.mat, cluster_rows = F, cluster_columns = F, name='Jaccard index',
                                               row_names_gp = gpar(fontsize = 9),
                                               column_names_gp = gpar(fontsize = 9),
                                               column_title_gp = gpar(fontsize = 11),
                                               col = col_fun,rect_gp = gpar(col= "white", lwd=0.3),row_names_side='left',column_names_rot = 0))
p.jaccard
ggsave(p.jaccard, filename = 'smFish_Cluster_type_Jaccard_index_heatmap.pdf', width=12, height =4)




## side by side reconciliation with the Seurat labels
## Rename cluster according to the size of cluster


uniq.cluster.RareQ = sort(unique(seu$cluster_ordered))
uniq.cluster.Seurat = sort(unique(seu_obj$cluster_ordered))


jaccard.mat <- matrix(0, nrow=length(uniq.cluster.RareQ), ncol=length(uniq.cluster.Seurat))
for(m in 1:dim(jaccard.mat)[1]){
  for(n in 1:dim(jaccard.mat)[2]){
    id.m = which(seu$cluster_ordered==uniq.cluster.RareQ[m])
    id.n = which(seu_obj$cluster_ordered==uniq.cluster.Seurat[n])

    jaccard.sim <- length(intersect(id.m, id.n))/length(union(id.m, id.n))
    jaccard.mat[m,n] <- jaccard.sim
  }
}
dimnames(jaccard.mat) <- list(uniq.cluster.RareQ, uniq.cluster.Seurat)
write.csv(jaccard.mat, file='Cluster_jaccard_matrix_RareQ_vs_Seurat_cluster.csv')
col_fun = circlize::colorRamp2(c(0,1), c("#FFF5F0",'#A50F15'))
#col_fun = circlize::colorRamp2(c(0,1), c("grey95", "#ff349c"))
p.jaccard <- as.ggplot(ComplexHeatmap::Heatmap(jaccard.mat, cluster_rows = F, cluster_columns = F, name='Jaccard index',
                                               row_names_gp = gpar(fontsize = 9),
                                               column_names_gp = gpar(fontsize = 9),
                                               column_title_gp = gpar(fontsize = 11),
                                               col = col_fun,rect_gp = gpar(col= "white", lwd=0.3),row_names_side='left',column_names_rot = 0))
p.jaccard
ggsave(p.jaccard, filename = 'smFish_Cluster_type_Jaccard_index_RareQ_vs_Seurat_heatmap.pdf', width=12, height =4)








## Test cluster using silhouette coefficient

library(Rcpp)

# Compile C++ code
sourceCpp("silhouette_rcpp.cpp")

silhouette_fast <- function(clusters, data, block_size = 1000) {
  # 1. Transform data to matrix
  if (!is.matrix(data)) {

    data <- as.matrix(data[, seq_len(ncol(data)), drop = FALSE])
    message("Data has been transformed into matrix")
  }

  # 2. Check input
  n <- nrow(data)
  if (length(clusters) != n) {
    stop("clusters size(", length(clusters), ") does not match with sample size(", n, ")")
  }
  if (block_size < 100) {
    warning("block_size should be ≥100 to boost efficiency, current value is", block_size)
  }

  # 3. Transform cluster labels into integers
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

Cluster.count <- sort(table(seu$Cluster))

## Rename cluster according to the size of cluster
Cluster.ordered <- (1:length(Cluster.count))[match(seu$Cluster, as.integer(names(Cluster.count)))]
seu$Cluster_ordered = Cluster.ordered

embed = seu_obj@reductions$pca@cell.embeddings



sil_width_Orig <- silhouette_fast(
  clusters = as.integer(seu$Cluster_ordered),
  data = embed
)

sil_width_Seurat <- silhouette_fast(
  clusters = as.integer(seu_obj$cluster_ordered),
  data = embed
)

# saveRDS(sil_width_Orig, file = 'sil_width_Orig.RDS')
# saveRDS(sil_width_Seurat, file = 'sil_width_Seurat.RDS')


sil_width_Orig <- readRDS('sil_width_Orig.RDS')
sil_width_Seurat <- readRDS('sil_width_Seurat.RDS')


sil.df <- data.frame(sil_width_Orig = sil_width_Orig,
                     cluster_Orig = seu$Cluster_ordered,
                     sil_width_Seurat = sil_width_Seurat,
                     cluster_Seurat = seu_obj$cluster_ordered)

p.sil.Orig <- ggplot(data=sil.df, aes(x=factor(cluster_Orig), y=sil_width_Orig)) + geom_boxplot(outlier.shape = NA, color='grey50', width=0.2) +
  stat_summary(fun.y = median, color='red', size=0.1, shape=15) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(yintercept = 0) +
  labs(x='Original cluster (with cluster size small to large)', y='Silhouette width') + scale_x_discrete(breaks = c(1, 20, 40,60, 80,100,120,140,160,180,200))
p.sil.Orig

p.sil.Seurat <- ggplot(data=sil.df, aes(x=factor(cluster_Seurat), y=sil_width_Seurat)) + geom_boxplot(outlier.shape = NA, color='grey50', width=0.2) +
  stat_summary(fun.y = median, color='red', size=0.1, shape=15) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_hline(yintercept = 0) +
  labs(x='Seurat cluster (with cluster size small to large)', y='Silhouette width') + scale_x_discrete(breaks = c(1, 20, 40,60, 80,100,120,140,160,180,200))
p.sil.Seurat

p.sil.comb = ggarrange(p.sil.Orig,
                       p.sil.Seurat, ncol=1)
ggsave(p.sil.comb, filename = 'smFish_Sil_boxplot_comb.pdf', width = 8, height = 5)









