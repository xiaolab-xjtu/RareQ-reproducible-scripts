setwd("/home/rstudio/Projects/Rare_cell/data/Xenium_Human_Kidney_RNA_Protein/")


library(anndata)
library(Seurat)
library(RareQ)
library(dplyr)
Assays = Seurat::Assays
h5dat = Read10X_h5('cell_feature_matrix.h5')

cell.info = read.csv('cells.csv.gz')

# Then, we preprocessed the dataset using the Scanpy pipeline70: normalized the total count of each cell to 1,000,
# log1p transformed the counts and scaled the transformed counts to Z scores (PMID: 38092912)

# Filter cells with low feature counts
good.id <- which(colSums(h5dat$`Gene Expression`) >= 10 & colSums(h5dat$`Protein Expression`) >= 10)


sc_object <- CreateSeuratObject(count=h5dat$`Gene Expression`[,good.id], project = "sc_object", min.cells = 3)
sc_object[["ADT"]] <- CreateAssayObject(counts = h5dat$`Protein Expression`[,good.id])

DefaultAssay(sc_object) <- 'RNA'
sc_object <- NormalizeData(sc_object, scale.factor = 80) %>% ScaleData()  #
sc_object <- RunPCA(sc_object, features = dimnames(sc_object)[[1]], npcs = 50)

DefaultAssay(sc_object) <- 'ADT'
sc_object <- NormalizeData(sc_object, normalization.method = "CLR", margin = 2)
sc_object <- ScaleData(sc_object, features = rownames(sc_object[["ADT"]]))
sc_object <- RunPCA(sc_object,
  features = rownames(sc_object[["ADT"]]),
  reduction.name = "pca_adt",
  reduction.key = "ADT_",
  npcs = 10
)



sc_object <- FindMultiModalNeighbors(sc_object, reduction.list = list("pca", "pca_adt"), dims.list = list(1:50, 1:10))
sc_object <- RunUMAP(sc_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims = 1:20)
DefaultAssay(sc_object) = 'RNA'
sc_object <- RunSPCA(sc_object, assay = 'RNA', graph = 'wsnn', features = dimnames(sc_object)[[1]])

sc_object <- FindNeighbors(object = sc_object,
                           k.param = 20,
                           compute.SNN = F,
                           prune.SNN = 0,
                           reduction = "spca",
                           dims = 1:30,
                           force.recalc = F, return.neighbor = T)

cluster = FindRare(sc_object = sc_object)
sc_object$cluster = cluster
sc_object$X = cell.info$x_centroid[good.id]
sc_object$Y = cell.info$y_centroid[good.id]

cluster.cnt <- sort(table(sc_object$cluster))
sc_object$cluster_sort = factor(as.character(sc_object$cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T)

# saveRDS(sc_object, file = 'Seu_obj.RDS')
sc_object <- readRDS('Seu_obj.RDS')

count.RNA <- sc_object@assays$RNA@counts
count.ADT <- sc_object@assays$ADT@counts
dimnames(count.ADT)[[1]] <- paste0('ADT_', dimnames(count.ADT)[[1]])
count.comb <- rbind(count.RNA, count.ADT)

ad = AnnData(t(count.comb))
# write_h5ad(ad, filename='matrix.h5ad')



cluster.df <- sc_object@meta.data

cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,30,31)])

#cluster.df$cluster_mod <- cluster.df$cluster # mod for color distinguish
#cluster.df$cluster_mod[cluster.df$cluster==82936] = max(cluster.df$cluster) + 10 # mod for color distinguish
p.cluster <- ggplot(data = cluster.df) + geom_point(aes(x = X, y = Y, color=factor(cluster_sort)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(cluster.df$cluster_sort))))
p.cluster


X.range <- diff(range(cluster.df$X))
Y.range <- diff(range(cluster.df$Y))

ggsave(p.cluster, filename = 'Xenium_Human_Kidney_cluster_all.pdf', width = X.range/25, height = Y.range/(25), units = 'mm')






cls = as.integer(names(sort(table(cluster))))
for(cl in cls){
  plot(sc_object$X,sc_object$Y,cex=0.1,col=ifelse(cluster==cl,'red','grey'))
}

Idents(sc_object) <- 'cluster_sort'

mk = FindAllMarkers(sc_object)
saveRDS(mk, file = 'Cluster_Marker_RNA.RDS')
mk = readRDS('Cluster_Marker_RNA.RDS')

top_markers <- mk %>%
  dplyr::group_by(cluster) %>%  # 按集群分组
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%  # 按log2FC降序
  dplyr::slice_head(n = 5)  # 取每个集群的前10个基因

p1 <- DotPlot(sc_object, features = unique(top_markers$gene)) + coord_flip()
ggsave(p1, filename = 'cluster_marker_dotplot.pdf', height = 20, width = 16)


meta.info <- sc_object@meta.data
meta.info$UMAP1 <- sc_object@reductions$wnn.umap@cell.embeddings[,1]
meta.info$UMAP2 <- sc_object@reductions$wnn.umap@cell.embeddings[,2]

p.DimPlot1 <- ggplot(data = meta.info) + geom_point(aes(x = UMAP1, y = UMAP2, color=factor(cluster_sort)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(meta.info$cluster_sort))))
p.DimPlot1

p.DimPlot <- DimPlot(sc_object, reduction = 'wnn.umap', group.by = 'cluster_sort',label = T, label.size = 9, raster = T,pt.size = 1, raster.dpi = c(1024, 1024)) + NoLegend() +
  scale_color_manual(values = getPalette(length(unique(sc_object$cluster_sort))))
p.DimPlot
ggsave(p.DimPlot, file='Xenium_Human_Kidney_UMAP_plot.pdf', width = 15, height = 15)



## Feature Plot
gs <- c('MET','KRT18',
        'COL5A2','ACTA2',
        'CD68', 'CD63',
        'CD3E','FOXP3','CD8A',
        'CD79A','MS4A1',
        'PECAM1','VWF')
gs <- c('MET',
        'CD68',
        'CD3E',
        'CD79A',
        'ACTA2',
        'VWF')
p.feature <- FeaturePlot(sc_object, features = gs, raster = T)
p.feature
ggsave(p.feature, filename = 'Xenium_Human_Kidney_Feature_plot.pdf', width = 8, height = 10)



plot.type <- function(type){
  plot(meta.info$X, meta.info$Y, cex=0.1, col=ifelse(meta.info$cluster_sort==type, 'red', 'grey'))

}



DefaultAssay(sc_object) = 'ADT'
p.ADT.DotPlot = DotPlot(sc_object, features = rownames(sc_object), group.by = 'cluster_sort') + coord_flip()
p.ADT.DotPlot
ggsave(p.ADT.DotPlot, filename = 'Xenium_Human_Kidney_cluster_ADT_marker_dotplot.pdf', height = 8, width = 16)

# Cycling: 1: Tumor; 5: Fibro; 7: CD4 T; 12: Plasma; 14: Tumor
# 3: cDC c(ITGAX)
# 2/4: Mast c(CPA3, KIT, MS4A2)
# 6: Monocyte
# 8: Lymphatic Endothelial Cells  c(PROX1, LYVE1, MMRN1)
# 9: Endothelial  c(PECAM1, VWF)
# 10/21: pDC   c('IL3RA','LILRA4')
# 18: cDC (LAMP3, CCR7, CD83) PD-L1 (ADT)

Idents(sc_object) <- 'cluster_sort'
GS1 <- c('TOP2A','MKI67','MET','KRT18','CPA3','KIT','ITGAX','COL5A2','ACTA2','FCN1','S100A12',
         'FOXP3','CTLA4','PROX1','LYVE1','PECAM1','VWF','IL3RA','LILRA4')
p.DotPlot.GS1 <- DotPlot(sc_object, assay = 'RNA', features = GS1, idents = c('1','2','3','4','5','6','7','8','9','10')) + coord_flip()
p.DotPlot.GS1

plot(meta.info$X, meta.info$Y, cex=0.1, col=ifelse(meta.info$cluster_sort %in% c('18'), 'red',ifelse(meta.info$cluster_sort %in% c('56'),'green', 'grey')))


p.Ki67 <- FeaturePlot(sc_object, features = c('Ki-67'), raster=T)
ggsave(p.Ki67, filename = 'Xenium_Human_Kidney_Ki67_featureplot.pdf', width = 15, height = 15)



p.RNA.Cycling <- DotPlot(sc_object, features = c('TOP2A','MKI67','MET','KRT18','COL5A2','ACTA2','FOXP3','CTLA4','MZB1','SDC1'),
                         idents = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14'), scale = F) + coord_flip()
p.ADT.Cycling <- DotPlot(sc_object, features = c('Ki-67', 'PCNA.1','CD4.1',''), idents = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14'), assay = 'ADT', scale = F) + coord_flip()
p.Cycling <- ggpubr::ggarrange(p.RNA.Cycling, p.ADT.Cycling, ncol=1, align='hv', heights = c(3,1.3))
ggsave(p.Cycling, filename = 'Xenium_Human_Kidney_Cycling_features_RNA_ADT_DotPlot.pdf', width = 8, height = 6)



p.RNA.DC <- DotPlot(sc_object, features = c('CD68','MS4A4A','MS4A6A',
                                            'CPA3','KIT','IL1RL1','MS4A2',
                                            'ITGAX','PLA2G7','IL2RA','KCNMA1',
                                            'S100A12','VCAN','FCN1','AQP9',
                                            'GZMB','TCL1A','MPEG1','LILRA4',
                                            'LAMP3','CCR7','CD83','IL7R'),
                         idents = c('2','3','4','6','10','13','16','18','20','21','25','40','55','59'), scale = F) + coord_flip()
p.ADT.DC <- DotPlot(sc_object, features = c('CD68.1','CD45RA.1','CD11c',
                                            'HLA−DR','CD163.1','GranzymeB','CD16',
                                            'PD-L1'),
                    assay = 'ADT',
                    idents = c('2','3','4','6','10','13','16','18','20','21','25','40','55','59'), scale = F) + coord_flip()
p.DC <- ggpubr::ggarrange(p.RNA.DC, p.ADT.DC, ncol=1, align='hv', heights = c(3,1.2))
ggsave(p.DC, filename = 'Xenium_Human_Kidney_DC_features_RNA_ADT_DotPlot.pdf', width = 8, height = 7)



## TLS
p.cluster.TLS <- ggplot(data = cluster.df) + geom_point(aes(x = X, y = Y, color=factor(cluster_sort)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(cluster.df$cluster_sort)))) +
  scale_x_continuous(limits = c(6800, max(cluster.df$X))) + scale_y_continuous(limits = c(min(cluster.df$Y), 2500))
p.cluster.TLS


ADT.genes <- c('CD8A.1', 'CD20', 'CD11c', 'CD4.1', 'PD-L1','CD31')
ADT.mk = sc_object@assays$ADT@data[ADT.genes,]


col.subset = cols[c(3,5,13,19,20,24,30)[1:6]]

cluster.df$gene_group <- apply(ADT.mk, 2, which.max)
cluster.df$gene_group <- factor(cluster.df$gene_group, levels=1:length(ADT.genes), labels = ADT.genes)
cluster.df.TLS <- cluster.df[cluster.df$X >= 6800 & cluster.df$X <= max(cluster.df$X) &
                               cluster.df$Y >= min(cluster.df$Y) & cluster.df$Y <= 2500,]

p.cluster.TLS.gene.nolegend <- ggplot(data = cluster.df.TLS) + geom_point(aes(x = X, y = Y, color=factor(gene_group)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = col.subset) +
  scale_x_continuous(limits = c(6800, max(cluster.df.TLS$X))) + scale_y_continuous(limits = c(min(cluster.df$Y), 2500))
p.cluster.TLS.gene.nolegend

p.cluster.TLS.gene.legend <- ggplot(data = cluster.df.TLS) + geom_point(aes(x = X, y = Y, color=factor(gene_group)), size=0.001) +
  theme_bw() + theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'right',legend.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  scale_color_manual(values = col.subset) +
  scale_x_continuous(limits = c(6800, max(cluster.df.TLS$X))) + scale_y_continuous(limits = c(min(cluster.df$Y), 2500))

p.TLS.comb <- ggpubr::ggarrange(p.cluster.TLS, p.cluster.TLS.gene.nolegend, ncol=1, align='hv')

X.range.TLS <- diff(c(6800, max(cluster.df.TLS$X)))
Y.range.TLS <- diff(c(min(cluster.df$Y), 2500))
ggsave(p.TLS.comb, filename = 'Xenium_Human_Kidney_TLS_combine.pdf', width = X.range.TLS/25, height = Y.range.TLS/25 * 2, units = 'mm')
ggsave(p.cluster.TLS.gene.legend, filename = 'Xenium_Human_Kidney_TLS_gene_legend.pdf', width = X.range.TLS/25, height = Y.range.TLS/25, units = 'mm')

