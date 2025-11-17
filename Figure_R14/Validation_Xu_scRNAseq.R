setwd('/home/rstudio/Projects/Rare_cell/data/Choroid_plexus/')

library(Seurat)
library(RareQ)
library(ggplot2)
library(ggpubr)



seu_obj = readRDS('seu_obj.RDS')
gene.all = read.table('Gene_names.tsv')
rownames(seu_obj@assays$RNA@data) = gene.all$V1[as.integer(rownames(seu_obj))+1]

label = readRDS('label.RDS')
RareQ.pred = readRDS('OUR_result.RDS')

seu_obj$cell_type = label
seu_obj$cluster = RareQ.pred

seu_obj$Epi_subtype = label
seu_obj$Epi_subtype[!grepl('Epi', seu_obj$Epi_subtype)] = 'Other'
seu_obj$Epi_subtype = factor(seu_obj$Epi_subtype, levels=c('Epi-1','Epi-2','Chil1+ Epi','Chil1+ Icam1+ Epi',
                                                           'Pcp4+ Epi','Pcp4+ Cntn1+ Epi','Other'),ordered = T)

cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)


getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])

p.Epi = DimPlot(seu_obj, group.by = 'Epi_subtype', label = T) + NoLegend() +
  scale_color_manual(values = getPalette(length(unique(seu_obj$Epi_subtype)))) +
  labs(title = 'Choroid plexus (Xu et al., 2024)')
p.Epi
ggsave(p.Epi, filename = '../smFish/smFish_Choroid_plexus_validation_Xu_data.pdf', width = 4, height = 4)


cluster.type <- tapply(seu_obj$Epi_subtype, seu_obj$cluster, function(x){
  cnt <- table(x)
  return(names(cnt)[which.max(cnt)])
})

seu_obj$Epi_cluster = seu_obj$cluster
seu_obj$Epi_cluster[seu_obj$Epi_cluster %in% names(cluster.type)[cluster.type=='Other']] = 0
cluster.cnt <- sort(table(seu_obj$Epi_cluster))
seu_obj$cluster_sort = as.character(factor(as.character(seu_obj$Epi_cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T))



p.Epi.cluster = DimPlot(seu_obj, group.by = 'cluster_sort', label = T) + NoLegend() +
  scale_color_manual(values = getPalette(length(unique(seu_obj$cluster_sort)))) +
  labs(title = 'Choroid plexus (Xu et al., 2024)')
p.Epi.cluster
ggsave(p.Epi.cluster, filename = '../smFish/smFish_Choroid_plexus_validation_cluster_Xu_data.pdf', width = 4, height = 4)


p.Epi.cluster.dotplot = DotPlot(seu_obj, features = c('Foxj1', 'Otx2','Lratd2','Mbp','Aldoc','Gfap','Ccdc151', 'Cldn5', 'Ttr'), group.by = 'cluster_sort', scale = T) + coord_flip()
ggsave(p.Epi.cluster.dotplot, filename = '../smFish/smFish_Choroid_plexus_validation_cluster_dotplot_Xu_data.pdf', width = 7, height = 2.5)


