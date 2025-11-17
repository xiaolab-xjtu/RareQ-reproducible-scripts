library(anndata)
library(reshape2)
library(ggplot2)

setwd('/home/rstudio/Projects/Rare_cell/data/')



eval.rare.OUR.RNA <- function(scenario){
  path <- paste0(scenario,'/')
  res <- readRDS(paste0(path, 'OUR_result.RDS'))
  labs <- readRDS(paste0(path, 'label.RDS'))
  type.cnt <- sort(table(labs))

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
                       Data=scenario)
  res.df$Method = 'RareQ_RNA'
  return(res.df)
}


eval.rare.OUR.ATAC <- function(scenario){
  path <- paste0(scenario,'/')
  res <- readRDS(paste0(path, 'OUR_result_ATAC.RDS'))
  labs <- readRDS(paste0(path, 'label.RDS'))
  type.cnt <- sort(table(labs))

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
                       Data=scenario)
  res.df$Method = 'RareQ_ATAC'
  return(res.df)
}




eval.rare.OUR.WNN <- function(scenario){
  path <- paste0(scenario,'/')
  res <- readRDS(paste0(path, 'OUR_result_WNN.RDS'))
  labs <- readRDS(paste0(path, 'label.RDS'))
  type.cnt <- sort(table(labs))

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
                       Data=scenario)
  res.df$Method = 'RareQ_WNN'
  return(res.df)
}





eval.rare.MarsGT <- function(scenario){
  path <- paste0(scenario,'/')
  res <- read.csv(paste0(path, 'pred_CPU_Optimal_Parameter.csv'), header = F)[,1]
  id_marsgt <- read.csv(paste0(path, 'Node_Ids_CPU_Optimal_Parameter.csv'), header=F)[,1]
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
                       Data=scenario)
  res.df$Method = 'MarsGT'
  return(res.df)
}




eval.rare.MarsGT.WNN <- function(scenario){
  path <- paste0(scenario,'/')
  res <- read.csv(paste0(path, 'pred_CPU_Parameter_Comparable.csv'), header = F)[,1]
  id_marsgt <- read.csv(paste0(path, 'Node_Ids_CPU_Parameter_Comparable.csv'), header=F)[,1]
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
                       Data=scenario)
  res.df$Method = 'MarsGT_WNN'
  return(res.df)
}




methods <- c('RareQ_RNA', 'RareQ_ATAC', 'RareQ_WNN', 'MarsGT', 'MarsGT_WNN')



## For one rare cluster
scenarios <-  c('ISSAAC_mCortex_RNA_ATAC',
                'Human_retina_RNA_ATAC',
                'Human_retina_rpe_choroid_RNA_ATAC',
                'Chen_2019_RNA_ATAC',
                'Mouse_Kidney_RNA_ATAC',
                '10x_Multiome_Pbmc10k_RNA_ATAC',
                'Mouse_gdT_RNA_ATAC',
                'Human_Gray_matter_RNA_ATAC',
                '10x_Multiome_PBMC_Chromium_RNA_ATAC',
                'Mouse_colon_RNA_ATAC'
                )

cols <- c('#FF0099', '#C3EF00', '#007ED3', '#FF9D1E', '#7FD2FF')

F1.summary <- function(k){
  labs <- readRDS(paste0(scenarios[k], '/label.RDS'))
  type.cnt <- sort(table(labs))
  rare_types <- names(type.cnt)[type.cnt <= sum(type.cnt) * 0.01]
  if(length(rare_types) < 1){
    rare_types <- names(type.cnt)[which.min(type.cnt)]
  }

  for(j in k){
    F1.OUR.RNA <- eval.rare.OUR.RNA(scenario = scenarios[j])
    F1.OUR.ATAC <- eval.rare.OUR.ATAC(scenario = scenarios[j])
    F1.OUR.WNN <- eval.rare.OUR.WNN(scenario = scenarios[j])
    F1.MarsGT <- eval.rare.MarsGT(scenario = scenarios[j])
    F1.MarsGT.WNN <- eval.rare.MarsGT.WNN(scenario = scenarios[j])
  }
  perf.mat <- rbind(F1.OUR.RNA, F1.OUR.ATAC, F1.OUR.WNN, F1.MarsGT, F1.MarsGT.WNN)
  #rownames(perf.mat) <- rare_types

  # plot.df <- melt(perf.mat)
  # names(plot.df) <- c('Type','Method','F1')

  return(perf.mat)
}

plot.data <- do.call(rbind, lapply(1:length(scenarios), F1.summary))
plot.data$Precision[is.na(plot.data$Precision)] <- 0
plot.data$Recall[is.na(plot.data$Recall)] <- 0



method.order <- c('RareQ_RNA','RareQ_ATAC','RareQ_WNN','MarsGT','MarsGT_WNN')
cols <- c('#FF0099', '#C3EF00', '#007ED3', '#FF9D1E', '#7FD2FF')

plot.data$Method <- factor(plot.data$Method, levels=method.order, ordered = T)

write.csv(plot.data, file = 'Figure/F1_precision_recall_10_multiome_Optimal_Parameter.csv')


p.F1 <- ggplot(data=plot.data, aes(x=Method, y=F1, fill=Method))+ geom_boxplot(width=0.5, outlier.shape = NA) + geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
  scale_fill_manual(values = cols) + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = 'black'),
        axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
        legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
  scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) + labs(title='F1 score',y=NULL) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))

p.Precision <- ggplot(data=plot.data, aes(x=Method, y=Precision, fill=Method)) + geom_boxplot(width=0.5, outlier.shape = NA) + geom_jitter(size=0.1, alpha=0.3, width = 0.1)+ #,outlier.shape=NA
  scale_fill_manual(values = cols) + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = 'black'),
        axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
        legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
  scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) +  labs(title='Precision',y=NULL) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))

p.Recall <- ggplot(data=plot.data, aes(x=Method, y=Recall, fill=Method)) +  geom_boxplot(width=0.5, outlier.shape = NA) + geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
  scale_fill_manual(values = cols) + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = 'black'),
        axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
        legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
  scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) +  labs(title='Recall',y=NULL) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))




p.10.multiome.datasets.list <- list(p.F1, p.Precision, p.Recall)


p.multiome.acc = ggpubr::ggarrange(p.F1, p.Precision, p.Recall, ncol=3, align='hv', common.legend = T, legend = 'right')
ggsave(p.multiome.acc, filename = 'Figure/10_Multiome_F1_Precision_Recall_plot_comb_Parameter_Comparable.pdf', width = 7, height = 2.5)







## clustering accuracy for abundant and rare cell types
methods <- c('RareQ_RNA','RareQ_ATAC','RareQ_WNN','MarsGT')


##
scenarios <-  c('ISSAAC_mCortex_RNA_ATAC',
                'Human_retina_RNA_ATAC',
                'Human_retina_rpe_choroid_RNA_ATAC',
                'Chen_2019_RNA_ATAC',
                'Mouse_Kidney_RNA_ATAC',
                '10x_Multiome_Pbmc10k_RNA_ATAC',
                'Mouse_gdT_RNA_ATAC',
                'Human_Gray_matter_RNA_ATAC',
                '10x_Multiome_PBMC_Chromium_RNA_ATAC',
                'Mouse_colon_RNA_ATAC'
)


NMI.mat <- matrix(NA, length(scenarios), 5)
for(j in 1:length(scenarios)){
  path <- paste0(scenarios[j],'/')
  lab.df <- readRDS(paste0(path, 'label.RDS'))
  OUR.res <- readRDS(paste0(path, 'OUR_result.RDS'))
  OUR.ATAC.res <- readRDS(paste0(path, 'OUR_result_ATAC.RDS'))
  OUR.WNN.res <- readRDS(paste0(path, 'OUR_result_WNN.RDS'))

  marsgt.res <- read.csv(paste0(path, 'pred_CPU_Optimal_Parameter.csv'), header = F)[,1]
  id_marsgt <- read.csv(paste0(path, 'Node_Ids_CPU_Optimal_Parameter.csv'), header=F)[,1]

  MarsGT.res <- marsgt.res[order(id_marsgt)]

  marsgt.res.wnn <- read.csv(paste0(path, 'pred_CPU_Parameter_Comparable.csv'), header = F)[,1]
  id_marsgt.wnn <- read.csv(paste0(path, 'Node_Ids_CPU_Parameter_Comparable.csv'), header=F)[,1]

  MarsGT.res.wnn <- marsgt.res.wnn[order(id_marsgt.wnn)]

  NMI.mat[j,] <- c(aricode::NMI(lab.df, OUR.res),
                   aricode::NMI(lab.df, OUR.ATAC.res),
                   aricode::NMI(lab.df, OUR.WNN.res),
                   aricode::NMI(lab.df, MarsGT.res),
                   aricode::NMI(lab.df, MarsGT.res.wnn))
}

colnames(NMI.mat) <- c('RareQ_RNA', 'RareQ_ATAC', 'RareQ_WNN', 'MarsGT', 'MarsGT_WNN')
rownames(NMI.mat) <- scenarios

write.csv(NMI.mat, file = 'Figure/NMI_10_Multiome_Optimal_Parameter.csv')

NMI.df <- data.frame(NMI=c(NMI.mat[,1],NMI.mat[,2],NMI.mat[,3],NMI.mat[,4],NMI.mat[,5]),
                     Method=c(rep(colnames(NMI.mat)[1], length(scenarios)),
                              rep(colnames(NMI.mat)[2], length(scenarios)),
                              rep(colnames(NMI.mat)[3], length(scenarios)),
                              rep(colnames(NMI.mat)[4], length(scenarios)),
                              rep(colnames(NMI.mat)[5], length(scenarios))),
                     Dataset=rep(colnames(NMI.mat), length(scenarios)))
NMI.df$Method <- factor(NMI.df$Method, levels=c('RareQ_RNA', 'RareQ_ATAC', 'RareQ_WNN', 'MarsGT', 'MarsGT_WNN'), ordered = T)



cols <- c('#FF0099', '#C3EF00', '#007ED3', '#FF9D1E', '#7FD2FF')
p.NMI <- ggplot(data=NMI.df, aes(x=Method, y=NMI, fill=Method)) + geom_boxplot(width=0.5, outlier.shape = NA) + geom_jitter(size=0.5, alpha=0.3, width = 0.1) + #,outlier.shape=NA
  scale_fill_manual(values = cols) + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, color = 'black'),
        axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
        legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
  scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0,0)) +  labs(y='NMI') +
  guides(color=guide_legend(nrow=1, byrow=TRUE))
p.NMI


ggsave(p.NMI, filename = 'Figure/10_Multiome_NMI_PBMC_ATAC_MarsGT_boxplot_Parameter_Comparable.pdf', width=3, height = 3.5)






