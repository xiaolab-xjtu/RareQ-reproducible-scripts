library(Seurat)
library(RareQ)

setwd('/home/rstudio/Projects/Rare_cell/data/Choroid_plexus/')


load('cor_matrix.Rdata')

c(4286, 12574)

seu_obj = readRDS('seu_obj.RDS')
label = readRDS('label.RDS')
RareQ.pred = readRDS('OUR_result.RDS')

seu_obj$cell_type = label
seu_obj$cluster = RareQ.pred

seu_obj$Epi_subtype = label
seu_obj$Epi_subtype[!grepl('Epi', seu_obj$Epi_subtype)] = 'Other'
seu_obj$Epi_subtype = factor(seu_obj$Epi_subtype, levels=c('Epi-1','Epi-2','Chil1+ Epi','Chil1+ Icam1+ Epi',
                                                           'Pcp4+ Epi','Pcp4+ Cntn1+ Epi','Other'),ordered = T)


cluster.type <- tapply(seu_obj$Epi_subtype, seu_obj$cluster, function(x){
  cnt <- table(x)
  return(names(cnt)[which.max(cnt)])
})

seu_obj$Epi_cluster = seu_obj$cluster
seu_obj$Epi_cluster[seu_obj$Epi_cluster %in% names(cluster.type)[cluster.type=='Other']] = 0
cluster.cnt <- sort(table(seu_obj$Epi_cluster))
seu_obj$cluster_sort = as.character(factor(as.character(seu_obj$Epi_cluster), levels=names(cluster.cnt), labels = 1:length(cluster.cnt), ordered = T))



smFish.obj = readRDS('../smFish/Seu_obj.RDS')

cluster.count <- sort(table(smFish.obj$cluster))

## Rename cluster according to the size of cluster
cluster.ordered <- (1:length(cluster.count))[match(smFish.obj$cluster, as.integer(names(cluster.count)))]
smFish.obj$cluster_ordered = cluster.ordered

df = smFish.obj@meta.data


plot.cluster = function(cl){

  cl.raw = unique(seu_obj$cluster[seu_obj$cluster_sort == cl])

  df$cor_eff = cor_matrix[,paste0('g',cl.raw)]

  p = ggplot(data=df) + geom_point(aes(x=X, y=Y, color=cor_eff), size=0.001) + scale_color_gradient(low=('#FFF5F0'), high=('#A50F15')) +
    labs(title = NULL, x=NULL, y=NULL) + theme_bw() +
    theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    annotate(
      "text", label=paste0('# ', cl),
      x = min(df$X)-5, y = max(df$Y-20), size = 5, colour = "black", hjust=0
    )


  return(p)
}



p.cor.comb = ggarrange(plot.cluster(1),
                   plot.cluster(2),
                   plot.cluster(3),
                   plot.cluster(4),
                   plot.cluster(5),
                   plot.cluster(6),
                   plot.cluster(7),
                   plot.cluster(8),
                   plot.cluster(9),
                   plot.cluster(10),
                   plot.cluster(11),
                   plot.cluster(12),
                   plot.cluster(13),
                   plot.cluster(14),
                   plot.cluster(15),
                   plot.cluster(16),
                   ncol=3, nrow=6)
p.cor.comb
X.range = diff(range(df$X))
Y.range = diff(range(df$Y))
ggsave(p.cor.comb, filename = 'smFish_Choroid_plexus_correlation_spatial_map.png',  width = X.range/(140 * 4) * 3, height = Y.range/(140 * 4) * 6, units = 'mm')



plot.cluster.legend = function(cl){

  cl.raw = unique(seu_obj$cluster[seu_obj$cluster_sort == cl])

  df$cor_eff = cor_matrix[,paste0('g',cl.raw)]

  p = ggplot(data=df) + geom_point(aes(x=X, y=Y, color=cor_eff), size=0.001) + scale_color_gradient(low=('#FFF5F0'), high=('#A50F15')) +
    labs(title = NULL, x=NULL, y=NULL) + theme_bw() +
    theme(legend.position = 'right', axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    annotate(
      "text", label=paste0('# ', cl),
      x = min(df$X)-5, y = max(df$Y-20), size = 5, colour = "black", hjust=0
    )


  return(p)
}

p.cor.show = plot.cluster.legend(3)
ggsave(p.cor.show, filename = 'smFish_Choroid_plexus_correlation_spatial_map_cluster3.pdf',  width = X.range/(140 * 4) * 1.2 , height = Y.range/(140 * 4), units = 'mm')


plot.dist = function(cl){

  cl.raw = unique(seu_obj$cluster[seu_obj$cluster_sort == cl])

  df$cor_eff = cor_matrix[,paste0('g',cl.raw)]
  df.subset = df[df$cluster_ordered %in% c(1, 11),]

  p = ggplot(data=df.subset, aes(x=factor(cluster_ordered), y=cor_eff, group=factor(cluster_ordered), fill=factor(cluster_ordered))) +
               geom_boxplot(alpha=0.5, outlier.shape = NA) +
    geom_jitter(size=0.1, alpha=1, width = 0.2, shape=21, stroke=0.3) +
    labs(title=paste0('# ', cl), x=NULL, y='Correlation') +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5), axis.text = element_text(colour = 'black'))

}

p.cor.dist.comb = ggarrange(plot.dist(1),
                            plot.dist(2),
                            plot.dist(3),
                            plot.dist(4),
                            plot.dist(5),
                            plot.dist(6),
                            plot.dist(7),
                            plot.dist(8),
                            plot.dist(9),
                            plot.dist(10),
                            plot.dist(11),
                            plot.dist(12),
                            plot.dist(13),
                            plot.dist(14),
                            plot.dist(15),
                            plot.dist(16), ncol=4, nrow=4, common.legend = T, legend='right')

p.cor.dist.comb




