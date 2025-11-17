library(anndata)
library(reshape2)

setwd('/home/rstudio/Projects/Rare_cell/data/Jurkat/Sensitivity/')

sample.matrix <- readRDS('sample_matrix.RDS')

## For one rare cluster
repeats <- 1:dim(sample.matrix)[1]

eval.rare.OUR <- function(rep_id){

  res <- readRDS(paste0('OUR/','DE_',rep_id,'.RDS'))
  #lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
  #labs <- lab.df$V1
  labs <- readRDS('Cell_group.RDS')
  type.cnt <- table(labs)

  if(T){
    rare_type <- names(type.cnt)[which.min(type.cnt)]
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
    return(c(F1))

  }
}


eval.rare.OUR.Q.ablation <- function(rep_id){

  res <- readRDS(paste0('OUR/','DE_',rep_id,'_Q_Ablation.RDS'))
  #lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
  #labs <- lab.df$V1
  labs <- readRDS('Cell_group.RDS')
  type.cnt <- table(labs)

  if(T){
    rare_type <- names(type.cnt)[which.min(type.cnt)]
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
    return(c(F1))

  }
}


methods <- c('Standard','Ablation (Q)')

perf.mat1 <- matrix(NA, nrow=length(repeats), ncol=length(methods))

for(k in 1:length(repeats)){

  F1.OUR <- eval.rare.OUR(rep_id = repeats[k])
  F1.OUR.Q.ablation <- eval.rare.OUR.Q.ablation(rep_id = repeats[k])
  perf.mat1[k,] <- c(F1.OUR, F1.OUR.Q.ablation)
}

colnames(perf.mat1) <- methods


library(ComplexHeatmap)
method.order <- c('Standard','Ablation (Q)')
cols <- c('#FF0099', '#C3EF00')
names(cols) <- method.order

perf.df <- melt(perf.mat1)
names(perf.df) <- c('Var1','Method','F1')
perf.df$index <- rep(rep(1:(dim(perf.mat1)[1]/10), each=10),length(method.order))

perf.df.ave <- tapply(perf.df$F1, list(perf.df$Method, perf.df$index), mean)
perf.df.ave.mod <- melt(perf.df.ave)
names(perf.df.ave.mod) <- c('Method', 'index', 'F1')
perf.df.ave.mod$Method <- factor(perf.df.ave.mod$Method, levels = method.order, ordered = T)

p.F1 <- ggplot(data=perf.df.ave.mod, aes(x=index, y=F1, color=factor(Method))) + geom_line() + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'right', legend.title = element_blank()) +
  scale_color_manual(values = cols) + labs(y='F1 score',x='Number of DE genes')
p.F1
ggsave(p.F1, filename = 'Result/F1_Jurkat_sensitivity_Q_ablation.pdf', width=6, height = 2)

saveRDS(p.F1, file = 'Result/F1_Jurkat_sensitivity_Q_ablation.RDS')



