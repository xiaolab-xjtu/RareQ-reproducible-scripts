setwd('/home/rstudio/Projects/Rare_cell/data/GSE303823_AD_DLB_PDD/')

library(Seurat)
library(RareQ)
library(dplyr)
library(anndata)
library(Matrix)
library(ggplot2)
library(magrittr)

seu_obj <- readRDS('GSE303823_rna_annotated.rds')
seu.obj.subset <- subset(seu_obj, subset=group %in% c('ADD','CTRL'))
seu.obj.subset <- subset(seu.obj.subset, subset=cell.type %in% c('Astrocyte', 'Endothelial cell', 'Ependymal cell', 'Microglia', 'OPC', 'Oligodendrocyte'))



AD.obj = seu.obj.subset
AD.obj <- FindNeighbors(object = AD.obj,
                      k.param = 20,
                      compute.SNN = F,
                      prune.SNN = 0,
                      reduction = "pca",
                      dims = 1:20,
                      force.recalc = F, return.neighbor = T)
AD.obj <- RunUMAP(AD.obj, dims = 1:20)
Assays = Seurat::Assays
cluster = FindRare(sc_object = AD.obj)
AD.obj$RareQ_cluster = cluster

cluster.cnt <- sort(table(AD.obj$RareQ_cluster))
AD.obj$cluster_sort = factor(as.character(AD.obj$RareQ_cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T)
Idents(AD.obj) <- 'cluster_sort'

# saveRDS(AD.obj, file = 'AD_obj.RDS')
# RNA_ad = AnnData(t(AD.obj@assays$RNA@counts))
# write_h5ad(RNA_ad, filename = 'RNA.h5ad')



AD.obj <- readRDS('AD_obj.RDS')

meta.info = AD.obj@meta.data[,c('sample','cell.type','group','sex','RareQ_cluster','cluster_sort')]
write.csv(meta.info, file = 'Meta_info.csv')


cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)

getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,30,31)])


p1 <- DimPlot(AD.obj, group.by='cell.type', label = T, label.size = 5) +
  scale_color_manual(values = getPalette(length(unique(AD.obj$cell.type)))) + labs(title='Main cell types') + NoLegend()
p2 <- DimPlot(AD.obj, group.by = 'group') +
  scale_color_manual(values = getPalette(length(unique(AD.obj$group)))) + labs(title='Group')
p3 <- DimPlot(AD.obj, label = T, label.size = 5) +
  scale_color_manual(values = getPalette(length(unique(AD.obj$RareQ_cluster)))) + labs(title='RareQ clusters') + NoLegend()
p1 + p2 + p3

p.umap <- p1 + p3
ggsave(p.umap, filename = 'UMAP_plot.pdf',  width = 20, height = 10)



cluster.type = tapply(AD.obj$cell.type, AD.obj$RareQ_cluster, function(x){
  cnt <- table(x)
  return(names(cnt)[which.max(cnt)])
})





# Microglia
MG.obj <- subset(AD.obj, subset=RareQ_cluster %in% as.integer(names(cluster.type)[cluster.type=='Microglia']))
MG.cluster.cnt <- sort(table(MG.obj$RareQ_cluster))
#MG.obj$cluster_sort = factor(as.character(MG.obj$RareQ_cluster), levels=names(MG.cluster.cnt), labels = 1:length(MG.cluster.cnt), ordered = T)
Idents(MG.obj) <- 'cluster_sort'
#MG.obj <- RunUMAP(MG.obj, reduction = 'pca', dims = 1:20)
DimPlot(MG.obj)

DefaultAssay(MG.obj) <- 'RNA'
MG.obj <- NormalizeData(MG.obj)
#mk.MG = FindAllMarkers(MG.obj)
#saveRDS(mk.MG, file = 'Microglia_Marker.RDS')

# saveRDS(MG.obj, file = 'Microglia_obj.RDS') #
MG.obj <- readRDS('Microglia_obj.RDS')
mk.MG <- readRDS('Microglia_Marker.RDS')

cluster.patient.cnt.MG <- table(as.character(MG.obj$cluster_sort), MG.obj$sample)
# write.csv(cluster.patient.cnt.MG, file = 'Microglia_cluster_patient_stat.csv')


MG_top_markers <- mk.MG %>%
  dplyr::group_by(cluster) %>%  
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%  
  dplyr::slice_head(n = 50)  
write.csv(MG_top_markers, file = 'MG_top50_markers.csv')
saveRDS(MG_top_markers, file = 'MG_top50_markers.RDS')

common.marker.MG <- c('CSF1R', 'MERTK','CYFIP1','C3','APOE','SPP1')
p.MG.marker = DotPlot(MG.obj, features = (c(common.marker.MG, unique(MG_top_markers$gene))), assay = 'RNA', cols = c('#FFF5F0', '#A50F15'), scale = F, group.by = 'cluster_sort') + coord_flip()
ggsave(p.MG.marker, filename = 'Microglia_top_marker.pdf', width = 4, height = 5)

MG.obj$group = as.character(MG.obj$group)
MG.cell.prop <- prop.table(table(as.character(MG.obj$cluster_sort), MG.obj$group), margin = 1)


MG.prop.df <- data.frame(Group=c(rep('AD', dim(MG.cell.prop)[1]), rep('CTRL', dim(MG.cell.prop)[1])),
                         Ratio=c(MG.cell.prop[,1], MG.cell.prop[,2]),
                         Cluster=as.character(rep(rownames(MG.cell.prop), 2)))
write.csv(MG.prop.df, file = 'AD_MG_proporation.csv')

p.prop.MG <- ggplot(data=MG.prop.df) + geom_bar(aes(x=(Cluster), y=Ratio, fill=factor(Group)), stat='identity',position = position_stack(), width = 0.5) +
  scale_fill_manual(values=cols[c(22,14)]) + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(),
        panel.grid=element_blank()) +
  scale_y_continuous(expand = c(0,0)) + labs(y='Proportion')
ggsave(p.prop.MG, filename = 'MG_cell_proportion.pdf', width = 3, height = 3)


col_fun = circlize::colorRamp2(seq(0,1, length.out=15), (cm.colors(30)[16:30]))
MG.cell.prop.pat <- prop.table(table(MG.obj$sample, as.character(MG.obj$cluster_sort)), margin = 2)
p.prop.ht.MG <- ggplotify::as.ggplot(Heatmap(MG.cell.prop.pat, cluster_rows = F, cluster_columns = F, col = col_fun, rect_gp = gpar(col = "white", lwd = 3), border_gp = gpar(col = "black"),
                                          name = 'Proportion'))
p.prop.ht.MG
write.csv(MG.cell.prop.pat, file = 'AD_MG_proportion_patnt.csv')

MG.cell.prop.pat.df <- data.frame(Prop=c(MG.cell.prop.pat[,1],MG.cell.prop.pat[,2], MG.cell.prop.pat[,3]),
                                  Patnt=rep(rownames(MG.cell.prop.pat), 3),
                                  Cluster=rep(colnames(MG.cell.prop.pat), each=dim(MG.cell.prop.pat)[1]))
p.prop.patnt.bar <- ggplot(data = MG.cell.prop.pat.df) + geom_bar(aes(x=Cluster, y=Prop, fill=factor(Patnt)), stat = 'identity', position = position_stack(), width = 0.5) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(name='RdBu', n=11)[c(5,4,3,2,7,8,9,10)]) + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(),
        panel.grid=element_blank()) +
  scale_y_continuous(expand = c(0,0)) + labs(y='Cell Proportion')
ggsave(p.prop.patnt.bar, filename = 'MG_proportion_patnt_barplot.pdf', width = 3, height = 3)

p.prop.comb.MG <- ggarrange(p.prop.MG, p.prop.ht.MG, ncol=1, align='v', heights = c(1,0.8))
ggsave(p.prop.comb.MG, filename = 'MG_proportion_comb.pdf', width = 3, height = 3)


##### GSEA
gsea.ana <- function(DE.res, category){

  m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = category) %>%
    dplyr::select(gs_name, entrez_gene)

  genename <- rownames(DE.res)
  gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("SYMBOL","ENTREZID"))
  colnames(gene_map)[1]<-"Gene"

  genelist_input <- data.frame(Gene=rownames(DE.res), logFC=DE.res$avg_log2FC)

  aaa<-inner_join(gene_map,genelist_input,by = "Gene")
  aaa<-aaa[,-1]
  aaa<-na.omit(aaa)
  aaa$logFC<-sort(aaa$logFC,decreasing = T)

  geneList = aaa[,1]
  names(geneList) = as.character(aaa[,1])
  geneList

  #Go_gseresult <- GSEA(geneList, TERM2GENE=m_t2g)

  Go_gseresult <- clusterProfiler::enricher(geneList, TERM2GENE=m_t2g)
  #Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
  #gseaplot2(Go_gseresult,1:5,pvalue_table = TRUE) 
  go.res <- Go_gseresult@result
  go.res <- go.res[go.res$p.adjust < 0.01,]
  #return(go.res)
  return(summary(Go_gseresult))
}


Micro.15.GO <- gsea.ana(mk.MG[mk.MG$cluster=='15' & mk.MG$avg_log2FC >= 1 & mk.MG$p_val_adj <= 0.01,], category = 'C5')
Micro.16.GO <- gsea.ana(mk.MG[mk.MG$cluster=='16' & mk.MG$avg_log2FC >= 1 & mk.MG$p_val_adj <= 0.01,], category = 'C5')
#Micro.21.GO <- gsea.ana(mk.MG[mk.MG$cluster=='21' & mk.MG$avg_log2FC >= 1 & mk.MG$p_val_adj <= 0.01,], category = 'C5')

Micro.GO <- rbind(Micro.15.GO, Micro.16.GO)
Micro.GO$Group = c(rep('15', dim(Micro.15.GO)[1]), rep('16', dim(Micro.16.GO)[1]))

plot.GO <- function(path){

  GO.dat <- Micro.GO[Micro.GO$Description %in% path,]
  GO.dat$Pathway = substr(GO.dat$Description, 6, nchar(GO.dat$Description))
  GO.dat$Pathway = stringr::str_to_sentence(gsub('_', ' ', GO.dat$Pathway, fixed = T))
  GO.dat1 <- GO.dat[GO.dat$Group=='15',]
  GO.dat2 <- GO.dat[GO.dat$Group=='16',]
  path.ord <- c(GO.dat2$Pathway[order(GO.dat2$qvalue, decreasing = T)],
                GO.dat1$Pathway[order(GO.dat1$qvalue, decreasing = T)])
  GO.dat$Pathway <- factor(GO.dat$Pathway, levels = path.ord, ordered = T)
  GO.dat$q_value = -log10(GO.dat$qvalue)
  write.csv(GO.dat, file = 'AD_MG_GO_dat.csv')

  p <- ggplot(data= GO.dat) + geom_bar(aes(x=q_value, y=Pathway, fill=factor(Group)), stat = 'identity', alpha=0.5, width = 0.75) + theme_bw() +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_manual(values = c('#8870ad', '#65A83E'), name='Microglia Cluster') +
    labs(x='-log10(q value)', y='GO terms') +
    geom_text(
      aes(label = Pathway, x=0, y = Pathway),  
      hjust = 0, 
      vjust = 0.5, 
      position = position_dodge(0.7),  
      color = "black", 
      size = 4  
    )
  p
}


path = c('GOBP_AXON_ENSHEATHMENT_IN_CENTRAL_NERVOUS_SYSTEM',
         'GOBP_ENSHEATHMENT_OF_NEURONS',
         'GOBP_GLIAL_CELL_DEVELOPMENT',
         'GOBP_GLIOGENESIS',
         'GOBP_MYELIN_MAINTENANCE',
         'GOBP_REGULATION_OF_PHAGOCYTOSIS',
         'GOBP_IRON_ION_TRANSMEMBRANE_TRANSPORT',
         'GOBP_CELLULAR_IRON_ION_HOMEOSTASIS',
         'GOBP_MYELOID_LEUKOCYTE_ACTIVATION',
         'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY')

p.MG.GO <- plot.GO(path = path)

ggsave(p.MG.GO, filename = 'MG_GO_barplot.pdf', width = 7, height = 4)



## MG (15) signature
DefaultAssay(AD.obj) = 'RNA'
Idents(AD.obj) = 'cluster_sort'
mk.MG.15 = FindMarkers(AD.obj, ident.1 = 15)
saveRDS(mk.MG.15, file = 'MG_15_marker.RDS')




# Astrocyte
AC.obj <- subset(AD.obj, subset=RareQ_cluster %in% as.integer(names(cluster.type)[cluster.type=='Astrocyte']))
Idents(AC.obj) <- 'cluster_sort'

DefaultAssay(AC.obj) <- 'RNA'
AC.obj <- NormalizeData(AC.obj)
# mk.AC = FindAllMarkers(AC.obj)
# saveRDS(mk.AC, file = 'AC_top_marker.RDS')
# saveRDS(AC.obj, file='Astrocyte_obj.RDS')

AC.obj = readRDS('Astrocyte_obj.RDS')
mk.AC <- readRDS('AC_top_marker.RDS')

cluster.patient.cnt.AC <- table(as.character(AC.obj$cluster_sort), AC.obj$sample)
# write.csv(cluster.patient.cnt.AC, file = 'Astrocyte_cluster_patient_stat.csv')


AC_top_markers <- mk.AC %>%
  dplyr::group_by(cluster) %>%  #
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::slice_head(n = 5)
write.csv(AC_top_markers, file = 'AC_top5_markers.csv')

p.AC.marker = DotPlot(AC.obj, features = c((unique(AC_top_markers$gene))), assay = 'RNA', scale = F, cols = c('#FFF5F0', '#A50F15'), group.by = 'cluster_sort') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
ggsave(p.AC.marker, filename = 'Astrocyte_top_marker.pdf', width = 10, height = 3)

AC.obj$group = as.character(AC.obj$group)
AC.cell.prop <- prop.table(table(as.character(AC.obj$cluster_sort), AC.obj$group), margin = 1)

AC.prop.df <- data.frame(Group=c(rep('AD', dim(AC.cell.prop)[1]), rep('CTRL', dim(AC.cell.prop)[1])),
                         Ratio=c(AC.cell.prop[,1], AC.cell.prop[,2]),
                         Cluster=as.character(rep(rownames(AC.cell.prop), 2)))
AC.prop.df$Cluster <- factor(AC.prop.df$Cluster, levels=c('5','7','12','18','19','20','24'), ordered = T)
write.csv(AC.prop.df, file = 'AD_AC_proportion.csv')


p.prop.AC <- ggplot(data=AC.prop.df) + geom_bar(aes(x=Cluster, y=Ratio, fill=factor(Group)), stat='identity',position = position_stack(), width = 0.6) +
  scale_fill_manual(values=cols[c(22,14)]) + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(),
        panel.grid=element_blank()) +
  scale_y_continuous(expand = c(0,0)) + labs(y='Proportion')
ggsave(p.prop.AC, filename = 'AC_cell_proportion.pdf', width = 3, height = 3)


col_fun = circlize::colorRamp2(seq(0,1, length.out=15), (cm.colors(30)[16:30]))
AC.cell.prop.pat <- prop.table(table(AC.obj$sample, as.character(AC.obj$cluster_sort)), margin = 2)
write.csv(AC.cell.prop.pat, file = 'AD_AC_proportion_patnt.csv')
p.prop.ht.AC <- ggplotify::as.ggplot(Heatmap(AC.cell.prop.pat[,c('5','7','12','18','19','20','24')], cluster_rows = F, cluster_columns = F, col = col_fun, rect_gp = gpar(col = "white", lwd = 3), border_gp = gpar(col = "black"),
                                             name = 'Proportion'))
p.prop.ht.AC


AC.cell.prop.pat.df <- data.frame(Prop=c(AC.cell.prop.pat[,1],AC.cell.prop.pat[,2], AC.cell.prop.pat[,3], AC.cell.prop.pat[,4], AC.cell.prop.pat[,5], AC.cell.prop.pat[,6],AC.cell.prop.pat[,7]),
                                  Patnt=rep(rownames(AC.cell.prop.pat), 7),
                                  Cluster=rep(colnames(AC.cell.prop.pat), each=dim(AC.cell.prop.pat)[1]))
AC.cell.prop.pat.df$Cluster <- factor(AC.cell.prop.pat.df$Cluster, levels=c('5','7','12','18','19','20','24'), ordered = T)
p.prop.patnt.bar.AC <- ggplot(data = AC.cell.prop.pat.df) + geom_bar(aes(x=Cluster, y=Prop, fill=factor(Patnt)), stat = 'identity', position = position_stack(), width = 0.6) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(name='RdBu', n=11)[c(5,4,3,2,7,8,9,10)]) + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(),
        panel.grid=element_blank()) +
  scale_y_continuous(expand = c(0,0)) + labs(y='Cell Proportion')
ggsave(p.prop.patnt.bar.AC, filename = 'AC_proportion_patnt_barplot.pdf', width = 3, height = 3)



p.prop.comb.AC <- ggarrange(p.prop.AC, p.prop.ht.AC, ncol=1, align='v', heights = c(1,0.8))
ggsave(p.prop.comb.AC, filename = 'AC_proportion_comb.pdf', width = 3, height = 3)




AC.5.GO <- gsea.ana(mk.AC[mk.AC$cluster=='5' & mk.AC$avg_log2FC >= 1 & mk.AC$p_val_adj <= 0.01,], category = 'C5')
AC.7.GO <- gsea.ana(mk.AC[mk.AC$cluster=='7' & mk.AC$avg_log2FC >= 1 & mk.AC$p_val_adj <= 0.01,], category = 'C5')
AC.12.GO <- gsea.ana(mk.AC[mk.AC$cluster=='12' & mk.AC$avg_log2FC >= 1 & mk.AC$p_val_adj <= 0.01,], category = 'C5')
AC.20.GO <- gsea.ana(mk.AC[mk.AC$cluster=='20' & mk.AC$avg_log2FC >= 1 & mk.AC$p_val_adj <= 0.01,], category = 'C5')
AC.24.GO <- gsea.ana(mk.AC[mk.AC$cluster=='24' & mk.AC$avg_log2FC >= 0.25 & mk.AC$p_val_adj <= 0.01,], category = 'C5')
AC.19.GO <- gsea.ana(mk.AC[mk.AC$cluster=='19' & mk.AC$avg_log2FC >= 0.25 & mk.AC$p_val_adj <= 0.01,], category = 'C5')

## GSEA analysis
data(c2BroadSets)

all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='C5')
#all_gene_sets = all_gene_sets[all_gene_sets$gs_subcat=='CP:KEGG',]

gcSample = split(all_gene_sets$gene_symbol,
                 all_gene_sets$gs_name)
names(gcSample)

gs = lapply(gcSample, unique)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
gsc


GSEA_res_AC <- gsva(as.matrix(AC.obj@assays$RNA@counts), gsc, min.sz=10, max.sz=500, method='ssgsea')
#saveRDS(GSEA_res_AC, file = 'AC_GSEA_data.RDS')
GSEA_res_AC <- readRDS('AC_GSEA_data.RDS')

sc_path <- CreateSeuratObject(count=GSEA_res_AC, project = "sc_path", min.cells = 0)
sc_path$cluster <- AC.obj$cluster_sort
Idents(sc_path) <- 'cluster'

path.astrocyte <- FindAllMarkers(sc_path, logfc.threshold = 0, min.pct = 0, return.thresh = 1.1, min.cells.group = 1)
top5.path.astrocyte <- path.astrocyte %>% group_by(cluster) %>% dplyr::top_n(n=20,wt=avg_log2FC)
#saveRDS(path.astrocyte, file = 'Astrocyte_pathway_top.RDS')
path.astrocyte =  readRDS('Astrocyte_pathway_top.RDS')

p.path.AC <- DotPlot(sc_path, features = unique(top5.path.astrocyte$gene)) + coord_flip()
ggsave(p.path.AC, filename = 'Astrocyte_pathway_dotplot.pdf', width = 15, height = 20)




plot.signature <- function(signature, title){

  plot.df <- data.frame(Score = unlist(GSEA_res_AC[which(dimnames(GSEA_res_AC)[[1]]==signature),]),
                        Type = sc_path$cluster)

  p <- ggplot(data=plot.df) + geom_boxplot(aes(x=Type, y=Score, fill=factor(Type))) + labs(title=title, x=NULL, y='Pathway enrichment score') +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.y = element_text(colour = 'black'), legend.position = 'none',
          plot.title = element_text(hjust=0.5), legend.title = element_blank()
    ) +
    scale_fill_manual(values = c('#B51D8D','#8870ad','#C72228','#ed8f84','#c19f70','#354E23','#FACB12'))
  p
}

p.path.AC.boxplot <- ggarrange(plot.signature(signature = 'GOBP_GLIAL_CELL_DERIVED_NEUROTROPHIC_FACTOR_RECEPTOR_SIGNALING_PATHWAY',
                         title='Glial cell-derived neurotrophic factor receptor signaling pathway'),
          plot.signature(signature = 'HP_ALZHEIMER_DISEASE',
                         title='Extracellular matrix structural constituent'),
          plot.signature(signature = 'GOBP_CHAPERONE_MEDIATED_AUTOPHAGY',
                         title='Chaperone-mediated autophagy'),
          #plot.signature(signature = 'HP_ALZHEIMER_DISEASE',
          #               title='Alzheimer disease'),
          plot.signature(signature = 'GOBP_REGULATION_OF_AMYLOID_FIBRIL_FORMATION',
                         title='Regulation of amyloid fibril formation'),
          nrow=1
)

ggsave(p.path.AC.boxplot, filename = 'AD_AC_path_boxplot.pdf',width = 15, height = 3)

gs <- c('CD44','OSMR','EMP1','SERPING1','HSPB1','CP','CD109','GPC5','CDH10') # Ref PMID: 28099414
p.AC.marker.Vln <- VlnPlot(AC.obj, features = gs, stack = T, flip = T,
                           cols = cols[1:length(gs)]) + NoLegend()
ggsave(p.AC.marker.Vln, filename = 'AC_Marker_Vln.pdf', width = 4, height = 3.5)


plot.df <- data.frame(Score1 = unlist(GSEA_res_AC[which(dimnames(GSEA_res_AC)[[1]]=='GOBP_CHAPERONE_MEDIATED_AUTOPHAGY'),]),
                      Score2 = unlist(GSEA_res_AC[which(dimnames(GSEA_res_AC)[[1]]=='GOBP_REGULATION_OF_AMYLOID_FIBRIL_FORMATION'),]),
                      Type = sc_path$cluster)

anova1 <- aov(Score1 ~ Type, data = plot.df)
anova2 <- aov(Score2 ~ Type, data = plot.df)
summary(anova1)
summary(anova2)





# Endothelial
EC.obj <- subset(AD.obj, subset=RareQ_cluster %in% as.integer(names(cluster.type)[cluster.type=='Endothelial cell']))
EC.cluster.cnt <- sort(table(EC.obj$RareQ_cluster))
Idents(EC.obj) <- 'cluster_sort'

DefaultAssay(EC.obj) <- 'RNA'
EC.obj <- NormalizeData(EC.obj)
#mk.EC = FindAllMarkers(EC.obj)
#saveRDS(mk.EC, file = 'Endothelial_Marker.RDS')

# saveRDS(EC.obj, file = 'Endothelial_obj.RDS') #
EC.obj <- readRDS('Endothelial_obj.RDS')
mk.EC <- readRDS('Endothelial_Marker.RDS')


EC_top_markers <- mk.EC %>%
  dplyr::group_by(cluster) %>%  # 按集群分组
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%  # 按log2FC降序
  dplyr::slice_head(n = 5)  # 取每个集群的前10个基因

write.csv(EC_top_markers, file = 'EC_top5_markers.csv')

p.EC.marker = DotPlot(EC.obj, features = (c(unique(EC_top_markers$gene))), assay = 'RNA', cols = c('#FFF5F0', '#A50F15'), scale = F, group.by = 'cluster_sort') + coord_flip()
ggsave(p.EC.marker, filename = 'Endothelial_top_marker.pdf', width = 5, height = 5)





# Oligodendrocyte (OLG)
OLG.obj <- subset(AD.obj, subset=RareQ_cluster %in% as.integer(names(cluster.type)[cluster.type=='Oligodendrocyte']))
OLG.cluster.cnt <- sort(table(OLG.obj$RareQ_cluster))
Idents(OLG.obj) <- 'cluster_sort'

DefaultAssay(OLG.obj) <- 'RNA'
OLG.obj <- NormalizeData(OLG.obj)
#mk.OLG = FindAllMarkers(OLG.obj)
#saveRDS(mk.OLG, file = 'OLG_Marker.RDS')

# saveRDS(OLG.obj, file = 'OLG_obj.RDS') #
OLG.obj <- readRDS('OLG_obj.RDS')
mk.OLG <- readRDS('OLG_Marker.RDS')


OLG_top_markers <- mk.OLG %>%
  dplyr::group_by(cluster) %>%  # 按集群分组
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%  # 按log2FC降序
  dplyr::slice_head(n = 5)  # 取每个集群的前10个基因

write.csv(OLG_top_markers, file = 'OLG_top5_markers.csv')

p.OLG.marker = DotPlot(OLG.obj, features = (c(unique(OLG_top_markers$gene))), assay = 'RNA', cols = c('#FFF5F0', '#A50F15'), scale = F, group.by = 'cluster_sort') + coord_flip()
ggsave(p.OLG.marker, filename = 'OLG_top_marker.pdf', width = 5, height = 5)






# OPC
OPC.obj <- subset(AD.obj, subset=RareQ_cluster %in% as.integer(names(cluster.type)[cluster.type=='OPC']))
OPC.cluster.cnt <- sort(table(OPC.obj$RareQ_cluster))
Idents(OPC.obj) <- 'cluster_sort'

DefaultAssay(OPC.obj) <- 'RNA'
OPC.obj <- NormalizeData(OPC.obj)
#mk.OPC = FindAllMarkers(OPC.obj)
#saveRDS(mk.OPC, file = 'OPC_Marker.RDS')

# saveRDS(OPC.obj, file = 'OPC_obj.RDS') #
OPC.obj <- readRDS('OPC_obj.RDS')
mk.OPC <- readRDS('OPC_Marker.RDS')


OPC_top_markers <- mk.OPC %>%
  dplyr::group_by(cluster) %>%  # 按集群分组
  dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%  # 按log2FC降序
  dplyr::slice_head(n = 5)  # 取每个集群的前10个基因

write.csv(OPC_top_markers, file = 'OPC_top5_markers.csv')

p.OPC.marker = DotPlot(OPC.obj, features = (c(unique(OPC_top_markers$gene))), assay = 'RNA', cols = c('#FFF5F0', '#A50F15'), scale = F, group.by = 'cluster_sort') + coord_flip()
ggsave(p.OPC.marker, filename = 'OPC_top_marker.pdf', width = 5, height = 5)



