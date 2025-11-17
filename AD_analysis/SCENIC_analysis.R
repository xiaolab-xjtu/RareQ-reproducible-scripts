setwd('/home/rstudio/Projects/Rare_cell/data/GSE303823_AD_DLB_PDD/')

library(Seurat)


## SCENIC
library(SCENIC)

MG.obj <- readRDS('Microglia_obj.RDS')

exprMat <- as.matrix(MG.obj@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
cellInfo <-  MG.obj@meta.data[,c('group','cluster_sort','sample')]

head(cellInfo)


### Initialize settings
# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="hgnc",
                                  dbDir="/home/rstudio/cell_type/SCENIC_db",
                                  nCores=4)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget")) # Toy run settings
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings`

## Save output
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

rm(list = ls())
library(Seurat)
library(SCENIC)
library(doParallel)

scenicOptions=readRDS(file="int/scenicOptions.Rds")

### Exploring output
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))

## Take IRF7 as example
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="IRF7"]
viewMotifs(tableSubset)

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="IRF7" & highConfAnnot==TRUE]
viewMotifs(tableSubset)


rm(list = ls())
gc()
library(Seurat)
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")




tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
AUCell::AUCell_plotTSNE(tSNE_scenic$Y,
                        exprMat,
                        aucell_regulonAUC, plots="Expression")
# pdf(file="AUCell_plotTSNE.pdf")
# par(mfrow=c(4,5))
# AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
# dev.off()


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cluster_sort),
                                     function(cells) rowMeans(AUCell:::getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


#





TF.subset <- c('STAT1 (30g)', 'FOXO3_extended (17g)','RUNX2 (37g)','FOXP1_extended (21g)',
               'ATF6_extended (21g)', 'RUNX2 (37g)', 'ETV6_extended (72g)','BRCA1 (16g)',
               'EZH2_extended (27g)', 'HDAC2_extended (17g)', 'NR2C2_extended (21g)',
               'CEBPG_extended (10g)','KDM5A_extended (13g)','JDP2_extended (12g)',
               'FOXO1_extended (11g)','NFIL3_extended (10g)', 'ELF1_extended (12g)',
               'ETS1_extended (17g)','KLF13_extended (13g)', 'BACH2_extended (20g)',
               'YY1 (10g)','NRF1 (28g)')

pdf(file="Microglia_regulonActivity_by_cluster_Scaled.pdf")
col = circlize::colorRamp2(breaks = seq(-1.5, 1.5, length.out=11), colors = rev(RColorBrewer::brewer.pal(name='RdBu', n=11)))
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[unique(TF.subset),c('15','16','21')], name="Regulon activity",col=col, )
dev.off()

dat = regulonActivity_byCellType_Scaled[unique(TF.subset),c('15','16','21')]
write.csv(dat, file = 'AD_SCENIC_HT_data.csv')

