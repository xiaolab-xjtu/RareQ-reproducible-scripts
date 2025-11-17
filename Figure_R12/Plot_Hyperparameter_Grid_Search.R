library(anndata)
library(reshape2)
library(ggplot2)

setwd('/home/rstudio/Projects/Rare_cell/data/')



eval.rare.MarsGT <- function(scenario, labsm, wd, lr){
  path <- paste0(scenario,'/')
  postfix = paste0('_labsm_', labsm, '_lr_', lr, '_wd_', wd)

  res <- read.csv(paste0(path, 'pred_GPU', postfix, '.csv'), header = F)[,1]
  id_marsgt <- read.csv(paste0(path, 'Node_Ids_GPU', postfix, '.csv'), header=F)[,1]
  labs <- readRDS(paste0(path, 'label.RDS'))
  type.cnt <- table(labs)
  res <- res[order(id_marsgt)]

  rare_types <- names(type.cnt)[type.cnt <= sum(type.cnt) * 0.01]
  F1.vec <- c()
  Precision.vec <- c()
  Recall.vec <- c()

  for(rare_type in rare_types){
    rare.cluster <- tapply(labs, res, function(x){
      type.cnt <- table(x)
      if(names(type.cnt)[which.max(type.cnt)]==rare_type){
        return(T)
      }else{
        return(F)
      }
    })

    true.lab <- ifelse(labs %in% rare_type, 1, 0)
    pred.lab <- ifelse(res %in% names(rare.cluster)[rare.cluster], 1, 0)

    TP <- sum(true.lab==1 & pred.lab==1)
    FP <- sum(true.lab==0 & pred.lab==1)
    FN <- sum(true.lab==1 & pred.lab==0)

    precision <- TP/(TP + FP)
    recall <- TP/(TP + FN)
    F1 <- TP/(TP + 0.5*(FN + FP))
    F1.vec <- c(F1.vec, F1)
    Precision.vec <- c(Precision.vec, precision)
    Recall.vec <- c(Recall.vec, recall)
  }

  res.df <- data.frame(Type=rare_types,
                       F1=F1.vec,
                       Precision=Precision.vec,
                       Recall=Recall.vec,
                       Data=scenario,
                       Method='MarsGT',
                       labsm=labsm,
                       lr=lr,
                       wd=wd)
  return(res.df)
}






methods <- c('RareQ_RNA', 'RareQ_ATAC', 'RareQ_WNN', 'MarsGT')

labsm_all = c('0', '0.1', '0.3')
wd_all = c('0', '0.1', '0.3')
lr_all = c('0.001', '0.0005')


## For one rare cluster
scenarios <-  c(#'ISSAAC_mCortex_RNA_ATAC',
                #'Human_retina_RNA_ATAC',
                'Human_retina_rpe_choroid_RNA_ATAC',
                #'Chen_2019_RNA_ATAC',
                #'Mouse_Kidney_RNA_ATAC',
                #'10x_Multiome_Pbmc10k_RNA_ATAC',
                'Mouse_gdT_RNA_ATAC',
                'Human_Gray_matter_RNA_ATAC',
                #'10x_Multiome_PBMC_Chromium_RNA_ATAC',
                'Mouse_colon_RNA_ATAC'
                )

cols <- c('#FF0099', '#C3EF00', '#007ED3', '#FF9D1E', '#7FD2FF')


res.df <- data.frame(Type=NA,
                     F1=NA,
                     Precision=NA,
                     Recall=NA,
                     Data=NA,
                     Method='MarsGT',
                     labsm=NA,
                     lr=NA,
                     wd=NA)
for(scenario in scenarios){
  for(labsm in labsm_all){
    for(wd in wd_all){
      for(lr in lr_all){

        acc = eval.rare.MarsGT(scenario = scenario, labsm = labsm, wd = wd, lr = lr)
        res.df <- rbind(res.df, acc)
      }
    }
  }
}


res.df = res.df[!is.na(res.df$Type),]
res.df$Precision[is.na(res.df$Precision)] = 0




p.acc.grid.search = ggplot(res.df, aes(x = lr, y = F1, fill = lr)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, color = "black") +
  facet_grid(
    labsm ~ wd,
    labeller = labeller(
      labsm = ~paste("labsm=", .x),
      wd = ~paste("wd=", .x)
    )
  ) +
  scale_fill_manual(values = c("#FF9999", "#66B2FF")) +
  labs(
    title = NULL,
    x = "wd", y = "F1 score",
    fill = "lr"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text.y = element_text(color = "darkred"),
    strip.text.x = element_text(color = "darkblue")
  )

p.acc.grid.search
ggsave(p.acc.grid.search, filename = 'Figure/MarsGT_Acc_Grid_Search_10_Multiome_data.pdf', width = 6, height = 4)



