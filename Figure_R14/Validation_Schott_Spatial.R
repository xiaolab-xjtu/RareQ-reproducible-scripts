setwd('/home/rstudio/Projects/Rare_cell/data/smFish/Validation/')

dat = read_h5ad('GSM7990097_e13_mouse_head.h5ad')


object = CreateSeuratObject(counts = t(dat$X))
object = NormalizeData(object)
object$X = -dat$obsm$spatial[,2]
object$Y = dat$obsm$spatial[,1]
object$cell_type = dat$obs$annotation

meta.info = object@meta.data

cols <- c("#532C8A","#c19f70","#f9decf","#c9a997","#B51D8D","#9e6762","#3F84AA","#F397C0",
          "#C594BF","#DFCDE4","#eda450","#635547","#C72228","#EF4E22","#f77b59","#989898",
          "#7F6874","#8870ad","#65A83E","#EF5A9D","#647a4f","#FBBE92","#354E23","#139992",
          "#C3C388","#8EC792","#0F4A9C","#8DB5CE","#1A1A1A","#FACB12","#C9EBFB","#DABE99",
          "#ed8f84","#005579","#CDE088","#BBDCA8","#F6BFCB"
)


getPalette = colorRampPalette(cols[c(2,3,1,5,6,7,8,9,11,12,13,14,16,18,19,20,22,23,24,27,28,29,30,31,34)])

p.cluster = ggplot(data=meta.info) + geom_point(aes(x=X, y=Y, color=factor(cell_type)), cex=0.1) +
  scale_color_manual(values = getPalette(length(unique(meta.info$cell_type))), guide=guide_legend(ncol=1, override.aes = list(size=2))) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
        legend.title = element_blank()) +
  labs(title='', x=NULL, y=NULL)

p.cluster

ggsave(p.cluster, filename = 'GSM7990097_cell_type_Dimplot.pdf', width = 9, height = 6)


plot.gene = function(gene, obj){

  df.base = data.frame(x=obj$X, y=obj$Y)
  X.range = diff(range(df.base$x))
  Y.range = diff(range(df.base$y))

  df.gene = df.base
  df.gene$gene = obj@assays$RNA@data[gene,]

  df.gene.sort = df.gene[order(df.gene$gene,decreasing = F),]

  p=ggplot(data=df.gene.sort) + (geom_point(aes(x=x,y=y,color=gene), size=0.1)) +
    scale_color_gradient(low=('#FFF5F0'), high=('#A50F15')) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = 'none') +
    annotate(
      "text", label =gene,
      x = min(df.gene.sort$x) + 2, y = max(df.gene.sort$y-2), size = 8, colour = "black", hjust=0
    )
  return(p)
}


p.gene.comb = ggarrange(plot.gene(gene='Foxj1', obj = object),
                        plot.gene(gene='Otx2', obj = object),
                        plot.gene(gene='Lratd2', obj = object),
                        plot.gene(gene='Ttr', obj = object),
                        plot.gene(gene='Cldn11', obj = object),
                        plot.gene(gene='Aqp1', obj = object),
                        nrow=2, ncol=3)

X.range = diff(range(meta.info$X))
Y.range = diff(range(meta.info$Y))

ggsave(p.gene.comb, filename = 'GSM7990097_Gene_plot.pdf', width = X.range/80*3, height = Y.range/80*2, units = 'mm')


