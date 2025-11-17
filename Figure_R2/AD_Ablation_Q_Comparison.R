setwd('/home/rstudio/Projects/Rare_cell/data/')

library(Seurat)
library(ggplot2)
library(RareQ)


cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)


getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])



## AD data
sc.obj2 <- readRDS('GSE303823_AD_DLB_PDD/AD_obj.RDS')
source('FindRare_Q_Ablation.R', encoding = 'utf-8')

res.standard = readRDS('GSE303823_AD_DLB_PDD/OUR_result.RDS')
res.Q.ablation = FindRare_Ablation_Q(sc.obj2)


sc.obj2$Standard = sc.obj2$cluster_sort
sc.obj2$Q_Ablation = res.Q.ablation


standard.cluster.count <- sort(table(sc.obj2$Standard))
Standard.ordered <- (1:length(standard.cluster.count))[match(sc.obj2$Standard, as.integer(names(standard.cluster.count)))]
sc.obj2$Standard_ordered = Standard.ordered

Q.Ablation.cluster.count <- sort(table(sc.obj2$Q_Ablation))
Q.Ablation.ordered <- (1:length(Q.Ablation.cluster.count))[match(sc.obj2$Q_Ablation, as.integer(names(Q.Ablation.cluster.count)))]
sc.obj2$Q_Ablation_ordered = Q.Ablation.ordered


p.DimPlot.standard.2 = DimPlot(sc.obj2, group.by = 'Standard_ordered', label = T) + labs(title = 'Standard') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(sc.obj2$Standard_ordered))), guide=guide_legend(ncol=1, override.aes = list(size=2))) + NoLegend()
p.DimPlot.Q.ablation.2 = DimPlot(sc.obj2, group.by = 'Q_Ablation_ordered', label = T) + labs(title = 'Ablation (Q)') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_color_manual(values = getPalette(length(unique(sc.obj2$Q_Ablation_ordered))), guide=guide_legend(ncol=1, override.aes = list(size=2))) + NoLegend()

p.DimPlot.comb.2 = ggpubr::ggarrange(p.DimPlot.standard.2, p.DimPlot.Q.ablation.2, nrow=1)

ggsave(p.DimPlot.comb.2, filename = 'GSE303823_AD_DLB_PDD/AD_Standard_Q_Ablation_compare.pdf', width = 8, height = 4)







uniq.cluster.Standard = sort(unique(sc.obj2$Standard_ordered))
uniq.cluster.Q.Ablation = sort(unique(sc.obj2$Q_Ablation_ordered))


jaccard.mat <- matrix(0, nrow=length(uniq.cluster.Standard), ncol=length(uniq.cluster.Q.Ablation))
for(m in 1:dim(jaccard.mat)[1]){
  for(n in 1:dim(jaccard.mat)[2]){
    id.m = which(sc.obj2$Standard_ordered==uniq.cluster.Standard[m])
    id.n = which(sc.obj2$Q_Ablation_ordered==uniq.cluster.Q.Ablation[n])

    jaccard.sim <- length(intersect(id.m, id.n))/length(union(id.m, id.n))
    jaccard.mat[m,n] <- jaccard.sim
  }
}
dimnames(jaccard.mat) <- list(uniq.cluster.Standard, uniq.cluster.Q.Ablation)
write.csv(jaccard.mat, file='GSE303823_AD_DLB_PDD/Cluster_jaccard_matrix_Q_Ablation.csv')
col_fun = circlize::colorRamp2(c(0,1), c("#FFF5F0",'#A50F15'))
p.jaccard <- as.ggplot(ComplexHeatmap::Heatmap(jaccard.mat, cluster_rows = F, cluster_columns = F, name='Jaccard index',
                                               row_names_gp = gpar(fontsize = 9),
                                               column_names_gp = gpar(fontsize = 9),
                                               column_title_gp = gpar(fontsize = 11),
                                               col = col_fun,rect_gp = gpar(col= "white", lwd=0.3),row_names_side='left',column_names_rot = 0))
p.jaccard
ggsave(p.jaccard, filename = 'GSE303823_AD_DLB_PDD/AD_Cluster_Q_Ablation_Jaccard_index_heatmap.pdf', width=5, height =4)




