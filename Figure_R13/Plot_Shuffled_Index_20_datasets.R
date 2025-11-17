library(anndata)
library(ggplot2)
library(ggpubr)

setwd('/home/rstudio/Projects/Rare_cell/data/')


Evaluate.parameter <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                  "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                  "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                  "UUOkidney", "Zelsel", "Choroid_plexus")

  eval.rare.OUR <- function(scenario, parameter.symbol, parameter.value){

    path <- paste0(scenario,'/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }

    labs <- readRDS(paste0(path, 'label.RDS'))
    type.cnt <- sort(table(labs))

    rare_types <- names(type.cnt)[type.cnt <= sum(type.cnt) * 0.01]
    if(length(rare_types) < 1){
      rare_types <- names(type.cnt)[which.min(type.cnt)]
    }

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
    res.df$para.symbol = parameter.symbol
    res.df$para.value = parameter.value
    return(res.df)
  }

  perf.df <- data.frame(Type=NA,
                        F1=NA,
                        Precision=NA,
                        Recall=NA,
                        Data=NA,
                        para.symbol=NA,
                        para.value=NA)
  for(j in 1:length(scenarios)){
    for(kk in 1:parameter.num){

      F1.OUR <- eval.rare.OUR(scenario = scenarios[j],
                              parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
      perf.df <- rbind(perf.df, F1.OUR)
    }
  }
  perf.df <- perf.df[2:dim(perf.df)[1],]
  perf.df$F1[is.na(perf.df$F1)] <- 0
  perf.df$Precision[is.na(perf.df$Precision)] <- 0
  perf.df$Recall[is.na(perf.df$Recall)] <- 0

  write.csv(perf.df, file = paste0('Figure/Test_robustness_20_datasets_', 'Reshuffle', '_F1.csv'))

  p.F1 <- ggplot(data=perf.df, aes(x=para.value, y=F1, group=factor(para.value)))+ geom_boxplot(width=0.5, outlier.shape = NA, fill='#7FD2FF') + #geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(color = 'black'),
          axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
          legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
    scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) + labs(title='F1 score',y=NULL) +
    guides(color=guide_legend(nrow=1, byrow=TRUE))

  p.Precision <- ggplot(data=perf.df, aes(x=para.value, y=Precision, group=factor(para.value))) + geom_boxplot(width=0.5, outlier.shape = NA, fill='#7FD2FF') + #geom_jitter(size=0.1, alpha=0.3, width = 0.1)+ #,outlier.shape=NA
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(color = 'black'),
          axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
          legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
    scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) +  labs(title='Precision',y=NULL) +
    guides(color=guide_legend(nrow=1, byrow=TRUE))

  p.Recall <- ggplot(data=perf.df, aes(x=para.value, y=Recall, group=factor(para.value))) +  geom_boxplot(width=0.5, outlier.shape = NA, fill='#7FD2FF') + #geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(color = 'black'),
          axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
          legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
    scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) +  labs(title='Recall',y=NULL) +
    guides(color=guide_legend(nrow=1, byrow=TRUE))

  p.20.datasets.comb <- ggarrange(p.F1, p.Precision, p.Recall, ncol=3, align='h', widths = c(0.8,0.8,0.8), common.legend = T, legend='none')
  return(p.20.datasets.comb)
}


ks <- c(0, 1:30) # 0 denote the unshuffled result


p.k <- Evaluate.parameter(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Figure/F1_score_20_datasets_shuffled_index.pdf'), width = 9, height = 3)
p.k
dev.off()







Evaluate.parameter.NMI <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                  "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                  "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                  "UUOkidney", "Zelsel", "Choroid_plexus")

  eval.rare.OUR <- function(scenario, rep_id, parameter.symbol, parameter.value){
    path <- paste0(scenario,'/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }
    labs <- readRDS(paste0(path, 'label.RDS'))
    type.cnt <- sort(table(labs))

    nmi.score = aricode::NMI(labs, res)
    res.df <- data.frame(NMI=nmi.score,
                         Data=scenario)
    res.df$para.symbol = parameter.symbol
    res.df$para.value = parameter.value

    return(res.df)
  }

  perf.df <- data.frame(NMI=NA,
                        Data=NA,
                        para.symbol=NA,
                        para.value=NA)
  for(j in 1:length(scenarios)){

    for(kk in 1:parameter.num){

      F1.OUR <- eval.rare.OUR(scenario = scenarios[j],
                              parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
      perf.df <- rbind(perf.df, F1.OUR)
    }
  }
  perf.df <- perf.df[2:dim(perf.df)[1],]

  write.csv(perf.df, file = paste0('Figure/Test_robustness_20_datasets_', 'Reshuffle', '_NMI.csv'))

  p.NMI <- ggplot(data=perf.df, aes(x=para.value, y=NMI, group=factor(para.value)))+ geom_boxplot(width=0.5, outlier.shape = NA, fill='#7FD2FF') + #geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(color = 'black'),
          axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
          legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
    scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) + labs(title='NMI',y=NULL) +
    guides(color=guide_legend(nrow=1, byrow=TRUE))

  return(p.NMI)
}

p.k.NMI <- Evaluate.parameter.NMI(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Figure/F1_score_20_datasets_shuffled_index_NMI.pdf'), width = 3, height = 3)
p.k.NMI
dev.off()









## Evaluate stability by Jaccard index between simulations
Evaluate.parameter.Jaccard <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                  "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                  "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                  "UUOkidney", "Zelsel", "Choroid_plexus")

  eval.rare.OUR <- function(scenario, parameter.symbol, parameter.value){

    path <- paste0(scenario,'/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }

    labs <- readRDS(paste0(path, 'label.RDS'))
    type.cnt <- sort(table(labs))

    res0 <- readRDS(paste0(path, 'OUR_result.RDS'))

    rare_types <- names(type.cnt)[type.cnt <= sum(type.cnt) * 0.01]
    if(length(rare_types) < 1){
      rare_types <- names(type.cnt)[which.min(type.cnt)]
    }

    Jaccard.vec <- c()

    for(rare_type in rare_types){
      rare.cluster <- tapply(labs, res, function(x){
        type.cnt <- table(x)
        if(names(type.cnt)[which.max(type.cnt)]==rare_type){
          return(T)
        }else{
          return(F)
        }
      })
      rare.cluster0 <- tapply(labs, res0, function(x){
        type.cnt <- table(x)
        if(names(type.cnt)[which.max(type.cnt)]==rare_type){
          return(T)
        }else{
          return(F)
        }
      })

      true.lab <- ifelse(labs %in% rare_type, 1, 0)
      pred0.lab <- ifelse(res0 %in% names(rare.cluster0)[rare.cluster0],1, 0)
      pred.lab <- ifelse(res %in% names(rare.cluster)[rare.cluster], 1, 0)

      id.pred = which(pred.lab==1)
      id.pred0 = which(pred0.lab==1)
      if(length(id.pred)==0 & length(id.pred0)==0){
        Jaccard.index = 1
      }else{
        Jaccard.index = length(intersect(id.pred, id.pred0))/length(union(id.pred, id.pred0))
      }

      Jaccard.vec <- c(Jaccard.vec, Jaccard.index)
    }

    res.df <- data.frame(Type=rare_types,
                         Jaccard=Jaccard.vec,
                         Data=scenario)
    res.df$para.symbol = parameter.symbol
    res.df$para.value = parameter.value
    return(res.df)
  }

  perf.df <- data.frame(Type=NA,
                        Jaccard=NA,
                        Data=NA,
                        para.symbol=NA,
                        para.value=NA)
  for(j in 1:length(scenarios)){
    for(kk in 1:parameter.num){

      Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j],
                              parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
      perf.df <- rbind(perf.df, Jaccard.OUR)
    }
  }
  perf.df <- perf.df[2:dim(perf.df)[1],]
  perf.df$Jaccard[is.na(perf.df$Jaccard)] <- 0


  write.csv(perf.df, file = paste0('Figure/Test_robustness_20_datasets_', 'Reshuffle', '_Jaccard.csv'))

  p.Jacccard <- ggplot(data=perf.df, aes(x=para.value, y=Jaccard, group=factor(para.value)))+ geom_boxplot(width=0.5, outlier.shape = NA, fill='#7FD2FF') + #geom_jitter(size=0.1, alpha=0.3, width = 0.1) + #,outlier.shape=NA
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(color = 'black'),
          axis.text.y = element_text(color = 'black'), plot.title = element_text(hjust=0.5),
          legend.title = element_blank(), axis.title.x = element_blank(), legend.position = 'top') +
    scale_y_continuous(limits = c(0,1), n.breaks = 3, expand=c(0.05,0)) + labs(title='Jaccard index',y=NULL) +
    guides(color=guide_legend(nrow=1, byrow=TRUE))

  return(p.Jacccard)
}


ks <- c(1:30)


p.k <- Evaluate.parameter.Jaccard(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Figure/Jaccard_20_datasets_shuffled_index.pdf'), width = 3, height = 3)
p.k
dev.off()













## ## Evaluate stability by Jaccard index between top 50 genes of rare cell types in 30 simulations
Evaluate.parameter.Jaccard.Top.Features <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c( "Airway", "Arc-ME", "B_lymphoma", "Cao", "Chen", "Cortex", "Heart",
                  "Kidney_ccRCC", "Kidney_normal", "Pediatric_gut", "MacParland", "Macosko",
                  "Mammary", "Pancreas", "Plasschaert", "Retina", "Shekhar",
                  "UUOkidney", "Zelsel", "Choroid_plexus")

  eval.rare.OUR <- function(scenario, parameter.symbol, parameter.value){

    path <- paste0(scenario,'/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'Top_marker_',parameter.symbol, '_', parameter.value,'.RDS'))
    }

    labs <- readRDS(paste0(path, 'label.RDS'))
    type.cnt <- sort(table(labs))

    res0 <- readRDS(paste0(path, 'Top_marker_ref.RDS'))

    rare_types <- unique(res0$type)

    Jaccard.vec <- c()

    for(rare_type in rare_types){

      mk0 = res0$marker[res0$type==rare_type]
      mk = res$marker[res$type==rare_type]

      Jaccard.index = length(intersect(mk0, mk))/length(union(mk0, mk))

      Jaccard.vec <- c(Jaccard.vec, Jaccard.index)
    }

    res.df <- data.frame(Type=rare_types,
                         Jaccard=Jaccard.vec,
                         Data=scenario)
    res.df$para.symbol = parameter.symbol
    res.df$para.value = parameter.value
    return(res.df)
  }

  perf.df <- data.frame(Type=NA,
                        Jaccard=NA,
                        Data=NA,
                        para.symbol=NA,
                        para.value=NA)
  for(j in 1:length(scenarios)){
    for(kk in 1:parameter.num){

      Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j],
                                   parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
      perf.df <- rbind(perf.df, Jaccard.OUR)
    }
  }
  perf.df <- perf.df[2:dim(perf.df)[1],]
  perf.df$Data_type = paste0(perf.df$Data,'_',perf.df$Type)

  perf.mat = dcast(perf.df, para.value~Data_type, value.var = 'Jaccard')
  perf.mat1 = as.matrix(perf.mat[,2:dim(perf.mat)[2]])

  rownames(perf.mat1) = parameter.list

  data.group = sapply(strsplit(colnames(perf.mat1), '_', fixed = T), '[', 1)

  perf.mat1.ord = do.call(cbind, tapply(1:dim(perf.mat1)[2], data.group,
                                        function(x){
                                          sub.mat = perf.mat1[,x,drop=F]
                                          sub.mat.ord = sub.mat[,order(colMeans(sub.mat), decreasing = T)]
                                          return(sub.mat.ord)
                                        }))

  col_fun = circlize::colorRamp2(seq(0,1, length.out=8), c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15"))
  ht = Heatmap(perf.mat1.ord, cluster_rows = F, cluster_columns = F, name = 'Jaccard index', show_column_names = F, col=col_fun, border = T,
          column_split = data.group)
  p.ht = ggplotify::as.ggplot(ht)
  return(p.ht)
}


ks <- c(1:30)


p.k <- Evaluate.parameter.Jaccard.Top.Features(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Figure/Jaccard_index_Top_Features_20_datasets_shuffled_index.pdf'), width = 6, height = 4)
p.k
dev.off()
