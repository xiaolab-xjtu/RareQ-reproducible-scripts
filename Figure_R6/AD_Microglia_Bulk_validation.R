setwd('/home/rstudio/Projects/Rare_cell/data/GSE303823_AD_DLB_PDD/Validation_data/')


normalize.data <- function(counts, meta, group){
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,       # Raw count
    colData = meta,           # Meta data
    design = ~ group          # Group column
  )
  # Calculate size factor
  dds <- estimateSizeFactors(dds)

  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  return(normalized_counts)
}



## GSE109887
GSE109887.raw <- read.table('GSE109887_MTG_non-normalized_GA_illumina_expression.txt/GSE109887_MTG_non-normalized_GA_illumina_expression.txt',sep='\t', header = T)

GSE109887.mat = as.matrix(GSE109887.raw[,grepl('Signal', colnames(GSE109887.raw))])
rownames(GSE109887.mat) = GSE109887.raw$ID_REF

GSE109887.group = ifelse(1:dim(GSE109887.mat)[2] %in% c(1:7,9,11,13:17,19:21,26,28,29,31,34,37,39:49,51,53,54,57,58,67,71:73,75,76,78), 'AD', 'Control')
GSE109887.meta <- data.frame(sample=colnames(GSE109887.mat), group=GSE109887.group)


GSE109887.norm.count = normalize.data(counts = round(GSE109887.mat), meta = GSE109887.meta, group = 'group')

saveRDS(GSE109887.mat, file = 'GSE109887_count.RDS')
saveRDS(GSE109887.norm.count, file = 'GSE109887_norm_count.RDS')
saveRDS(GSE109887.meta, file = 'GSE109887_meta.RDS')



library(GSVA)

mk.MG.15 <- readRDS('../MG_15_marker.RDS')
GSE109887.meta <- readRDS('GSE109887_meta.RDS')
GSE109887.norm.count <- readRDS('GSE109887_norm_count.RDS')
gsc = list('Path1'=rownames(mk.MG.15)[order(mk.MG.15$avg_log2FC, decreasing = T)[1:50]])


GSEA_MG_AD <- gsva(GSE109887.norm.count, gsc, min.sz=10, max.sz=500, method='ssgsea')
MG.sign.ES.df = data.frame(ES=GSEA_MG_AD[1,], group=GSE109887.meta$group)
wilcox.test(MG.sign.ES.df$ES~MG.sign.ES.df$group)

p.MG.ES.boxplot <- ggplot(data=MG.sign.ES.df, aes(x=factor(group), y=ES, fill=factor(group))) + geom_boxplot(outlier.shape = NA, alpha=0.6) + geom_jitter(size=1.5, alpha=1, width = 0.2, shape=21, stroke=0.3) +
  theme_bw() + theme(panel.border = element_rect(fill=NA), panel.grid = element_blank(), legend.position = 'none', axis.text = element_text(colour = 'black'),
                     plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(name='RdBu', n=11)[c(3,9)]) +
  labs(x=NULL, y='MG (15) signature enrichment score', title='p = 5.93 * 10-5')
p.MG.ES.boxplot
ggsave(p.MG.ES.boxplot, filename = 'MG_15_ES_boxplot_AD.pdf', width = 1.2, height = 3)



