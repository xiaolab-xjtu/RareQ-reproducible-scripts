setwd('/home/rstudio/Projects/Rare_cell/data/Xenium/Validation/Visium_10X/Data2/')

library(arrow)
library(Seurat)
library(RareQ)
library(ggplot2)
library(ggrastr)

## Read in
#tissue_pos <- read_parquet("binned_outputs/square_008um/spatial/tissue_positions.parquet")

# write as .csv 文件
# write.table(
#   tissue_pos,
#   "binned_outputs/square_008um/spatial/tissue_positions.csv",
#   col.names = TRUE,
#   row.names = FALSE,
#   quote = FALSE,
#   sep = ","
# )


object <- Load10X_Spatial(data.dir = 'binned_outputs/square_008um/')

count.plot <- SpatialFeaturePlot(object, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size.factor=0.3) + theme(legend.position = "right")
count.plot

object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object)
# object <- FindNeighbors(object, dims = 1:50)
# object <- FindClusters(object, resolution = 0.5) # Very slow and over clustering
# object <- RunUMAP(object, dims = 1:50)

SpatialFeaturePlot(object, features = c("Rspo2",'Nwd2','Strip2','Stard5'), pt.size.factor = 0.01)
SpatialFeaturePlot(object, features = c("Rspo2"), pt.size.factor = 0.01)
SpatialFeaturePlot(object, features = c("Plcxd2"), pt.size.factor = 0.01)
SpatialFeaturePlot(object, features = c("Nwd2"), pt.size.factor = 0.01)
SpatialFeaturePlot(object, features = c("Strip2"), pt.size.factor = 0.01)
SpatialFeaturePlot(object, features = c("Stard5"), pt.size.factor = 0.01)



based.df <- data.frame(x=object@images$slice1@coordinates$row,
                       y=object@images$slice1@coordinates$col)

X.range = diff(range(based.df$x))
Y.range = diff(range(based.df$y))

plot.gene = function(gene){
  df = based.df
  df$gene = object@assays$Spatial@data[gene,]
  df.sort = df[order(df$gene,decreasing = F),]

  p=ggplot(data=df.sort) + (geom_point(aes(x=x,y=y,color=gene), size=0.5)) +
    scale_color_gradient(low=('#FFF5F0'), high=('#A50F15')) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = 'none') +
    annotate(
      "text", label =gene,
      x = min(df.sort$x) + 2, y = max(df.sort$y-2), size = 8, colour = "black", hjust=0
    )
  return(p)
}

p.comb = ggarrange(plot.gene('Cpne4'), plot.gene('Stard5'), plot.gene('Strip2'),
                   plot.gene('Rspo2'), plot.gene('Plcxd2'), plot.gene('Nwd2'),
                   plot.gene('Dner'), plot.gene('Sema5b'),plot.gene('Meis2'),
                   ncol=3, nrow=3)

ggsave(p.comb, filename = 'p_comb_gene.pdf', width = X.range/5 * 3, height = Y.range/(5) * 3, limitsize = F, units = 'mm')
ggsave(p.comb, filename = 'p_comb_gene.tiff', width = X.range/5 * 3, height = Y.range/(5) * 3, limitsize = F, units = 'mm')






