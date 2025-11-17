library(anndata)

setwd('/home/rstudio/Projects/Rare_cell/data/PBMCs/')


Evaluate.parameter.F1 <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10") # ,
  repeats <- c(0:49)

  eval.rare.OUR <- function(scenario, rep_id, parameter.symbol, parameter.value){
    path <- paste0('site1/',scenario,'/',rep_id,'/R_input/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }

    lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
    labs <- lab.df$V1
    type.cnt <- table(labs)

    if(scenario %in% c('Rare1_Ordinary4', 'Rare1_Ordinary9')){
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

    }else{
      F1.vec <- c()
      rare_types <- names(type.cnt)[order(type.cnt)[1:5]]
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
      }
      return(F1.vec)
    }
  }


  perf.mat1 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 1){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in c(1:parameter.num)){
        if(parameter.list[kk]==0){
          F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=0)
        }else{
          F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        }

        acc.vec <- c(acc.vec, F1.OUR)
      }
      perf.mat1[k,] <- acc.vec
    }
  }

  colnames(perf.mat1) <- parameter.list


  perf.mat2 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 2){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in 1:parameter.num){

        F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        acc.vec <- c(acc.vec, F1.OUR)
      }
      perf.mat2[k,] <- acc.vec
    }
  }

  colnames(perf.mat2) <- parameter.list

  perf.mat3 <- matrix(NA, nrow=length(repeats), ncol=parameter.num*5)
  for(j in 3){
    for(k in 1:length(repeats)){
      R1.vec <- c()
      R2.vec <- c()
      R3.vec <- c()
      R4.vec <- c()
      R5.vec <- c()
      for(kk in 1:parameter.num){

        F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        R1.vec <- c(R1.vec, F1.OUR[1])
        R2.vec <- c(R2.vec, F1.OUR[2])
        R3.vec <- c(R3.vec, F1.OUR[3])
        R4.vec <- c(R4.vec, F1.OUR[4])
        R5.vec <- c(R5.vec, F1.OUR[5])
      }
      perf.mat3[k,] <- c(R1.vec, R2.vec, R3.vec, R4.vec, R5.vec)
    }
  }

  #colnames(perf.mat3) <- c(methods)

  library(ComplexHeatmap)
  ord1 <- order(rowMeans(perf.mat1), decreasing = T)
  ord2 <- order(rowMeans(perf.mat2), decreasing = T)
  ord3 <- order(rowMeans(perf.mat3), decreasing = T)

  col_fun = circlize::colorRamp2(c(0,0.5,1), c("grey95",'#FF99A3', "red"))
  col_fun = circlize::colorRamp2(seq(0,1, length.out=8), c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15"))
  #col_fun = circlize::colorRamp2(c(0,0.5,1), c("#4DADCF",'yellow', "red"))

  method.order <- c('RareQ','scCAD','CellSIUS','RaceID','GiniClust2','FiRE','EDGE','GapClust')
  cols <- c('#FF0099', '#C3EF00', '#007ED3','#FF9D1E','#7FD2FF', '#00C19B', '#894FC6', '#D55E00')
  names(cols) <- method.order

  perf.mat1.ord <- perf.mat1[ord1,]
  perf.mat2.ord <- perf.mat2[ord2,]
  perf.mat3.ord <- perf.mat3[ord3,]

  perf.mat1.out <- perf.mat1.ord
  perf.mat2.out <- perf.mat2.ord
  perf.mat3.out <- perf.mat3.ord
  colnames(perf.mat3.out) <- rep(c('R1', 'R2', 'R3', 'R4', 'R5'), each=parameter.num)

  write.csv(perf.mat1.out, file = paste0('Result/Sim-PBMC-1_', 'Reshuffle', '.csv'))
  write.csv(perf.mat2.out, file = paste0('Result/Sim-PBMC-2_', 'Reshuffle', '.csv'))
  write.csv(perf.mat3.out, file = paste0('Result/Sim-PBMC-3_', 'Reshuffle', '.csv'))

  row_ha1 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(F1)' = anno_barplot(colMeans(perf.mat1.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha2 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(F1)' = anno_barplot(colMeans(perf.mat2.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha3 = rowAnnotation('Types' = rep(c('R1','R2','R3','R4','R5'), each=parameter.num),  'Mean(F1)' = anno_barplot(colMeans(perf.mat3.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90','R2'='grey70','R3'='grey50','R4'='grey30','R5'='grey10')),
                          annotation_width=c(0.25,1))

  p1 <- Heatmap(t(perf.mat1.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='F1 score',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha1, row_names_rot = 270)
  p2 <- Heatmap(t(perf.mat2.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='F1 score',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha2, row_names_rot = 270)
  p3 <- Heatmap(t(perf.mat3.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='F1 score',
                row_names_side = 'left', row_title_rot = 0, row_dend_width = unit(0, "mm"), right_annotation = row_ha3, show_row_names = T,
                row_split = factor(c(rep('R1', parameter.num),
                                     rep('R2', parameter.num),
                                     rep('R3', parameter.num),
                                     rep('R4', parameter.num),
                                     rep('R5', parameter.num)), levels=c('R1', 'R2', 'R3', 'R4', 'R5'), ordered = T))
  p.F1 <- (p1 %v% p2 %v% p3)
  return(p.F1)
}



ks <- c(0, 1:30) # 0 denote the unshuffled result


p.k <- Evaluate.parameter.F1(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Result/F1_score_PBMC_Sim123_shuffled_index.pdf'), width = 4, height = 6)
p.k
dev.off()







Evaluate.parameter.NMI <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10") # ,
  repeats <- c(0:49)

  eval.rare.OUR <- function(scenario, rep_id, parameter.symbol, parameter.value){
    path <- paste0('site1/',scenario,'/',rep_id,'/R_input/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }
    lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
    labs <- lab.df$V1
    type.cnt <- table(labs)

    nmi.score = aricode::NMI(lab.df$V1, res)
    return(nmi.score)
  }


  perf.mat1 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 1){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in 1:parameter.num){

        F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        acc.vec <- c(acc.vec, F1.OUR)
      }
      perf.mat1[k,] <- acc.vec
    }
  }

  colnames(perf.mat1) <- parameter.list


  perf.mat2 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 2){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in c(1:parameter.num)){
        if(parameter.list[kk]==0){
          F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=0)
        }else{
          F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        }

        acc.vec <- c(acc.vec, F1.OUR)
      }
      perf.mat2[k,] <- acc.vec
    }
  }

  colnames(perf.mat2) <- parameter.list

  perf.mat3 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 2){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in 1:parameter.num){

        F1.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        acc.vec <- c(acc.vec, F1.OUR)
      }
      perf.mat3[k,] <- acc.vec
    }
  }

  colnames(perf.mat3) <- parameter.list

  library(ComplexHeatmap)
  ord1 <- order(rowMeans(perf.mat1), decreasing = T)
  ord2 <- order(rowMeans(perf.mat2), decreasing = T)
  ord3 <- order(rowMeans(perf.mat3), decreasing = T)

  col_fun = circlize::colorRamp2(c(0,0.5,1), c("grey95",'#FF99A3', "red"))
  col_fun = circlize::colorRamp2(seq(0,1, length.out=8), c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15"))
  #col_fun = circlize::colorRamp2(c(0,0.5,1), c("#4DADCF",'yellow', "red"))

  method.order <- c('RareQ','scCAD','CellSIUS','RaceID','GiniClust2','FiRE','EDGE','GapClust')
  cols <- c('#FF0099', '#C3EF00', '#007ED3','#FF9D1E','#7FD2FF', '#00C19B', '#894FC6', '#D55E00')
  names(cols) <- method.order

  perf.mat1.ord <- perf.mat1[ord1,]
  perf.mat2.ord <- perf.mat2[ord2,]
  perf.mat3.ord <- perf.mat3[ord3,]

  perf.mat1.out <- perf.mat1.ord
  perf.mat2.out <- perf.mat2.ord
  perf.mat3.out <- perf.mat3.ord

  write.csv(perf.mat1.out, file = paste0('Result/Sim-PBMC-1_', 'Reshuffle', '_NMI.csv'))
  write.csv(perf.mat2.out, file = paste0('Result/Sim-PBMC-2_', 'Reshuffle', '_NMI.csv'))
  write.csv(perf.mat3.out, file = paste0('Result/Sim-PBMC-3_', 'Reshuffle', '_NMI.csv'))

  row_ha1 = rowAnnotation('Mean(F1)' = anno_barplot(colMeans(perf.mat1.ord), ylim=c(0,1)),
                          show_annotation_name = F)
  row_ha2 = rowAnnotation('Mean(F1)' = anno_barplot(colMeans(perf.mat2.ord), ylim=c(0,1)),
                          show_annotation_name = F)
  row_ha3 = rowAnnotation('Mean(F1)' = anno_barplot(colMeans(perf.mat2.ord), ylim=c(0,1)),
                          show_annotation_name = F)

  p1 <- Heatmap(t(perf.mat1.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='NMI',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha1, row_names_rot = 270)
  p2 <- Heatmap(t(perf.mat2.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='NMI',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha2, row_names_rot = 270)
  p3 <- Heatmap(t(perf.mat3.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='NMI',
                row_names_side = 'left', row_title_rot = 0, row_dend_width = unit(0, "mm"), right_annotation = row_ha3, show_row_names = T,
                row_names_rot = 270)
  p.F1 <- (p1 %v% p2 %v% p3)
  return(p.F1)

}

p.k.NMI <- Evaluate.parameter.NMI(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Result/F1_score_PBMC_Sim123_shuffled_index_NMI.pdf'), width = 4, height = 6)
p.k.NMI
dev.off()








## Evaluate Jaccard index
Evaluate.parameter.Jaccard <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10") # ,
  repeats <- c(0:49)

  eval.rare.OUR <- function(scenario, rep_id, parameter.symbol, parameter.value){
    path <- paste0('site1/',scenario,'/',rep_id,'/R_input/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'OUR_result.RDS'))
    }else{
      res <- readRDS(paste0(path, 'OUR_result_',parameter.symbol, '_', parameter.value,'.RDS'))
      ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      res <- res[order(ids)]
    }

    lab.df <- read.table(paste0(path, 'sub_label.txt'), sep='\t', header=F)
    labs <- lab.df$V1
    type.cnt <- table(labs)

    res0 <- readRDS(paste0(path, 'OUR_result.RDS'))

    if(scenario %in% c('Rare1_Ordinary4', 'Rare1_Ordinary9')){
      rare_type <- names(type.cnt)[which.min(type.cnt)]
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
      return(c(Jaccard.index))

    }else{
      Jaccard.vec <- c()
      rare_types <- names(type.cnt)[order(type.cnt)[1:5]]
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
      return(Jaccard.vec)
    }
  }


  perf.mat1 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 1){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in c(1:parameter.num)){
        if(parameter.list[kk]==0){
          Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=0)
        }else{
          Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                  parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        }

        acc.vec <- c(acc.vec, Jaccard.OUR)
      }
      perf.mat1[k,] <- acc.vec
    }
  }

  colnames(perf.mat1) <- parameter.list


  perf.mat2 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 2){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in 1:parameter.num){

        Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        acc.vec <- c(acc.vec, Jaccard.OUR)
      }
      perf.mat2[k,] <- acc.vec
    }
  }

  colnames(perf.mat2) <- parameter.list

  perf.mat3 <- matrix(NA, nrow=length(repeats), ncol=parameter.num*5)
  for(j in 3){
    for(k in 1:length(repeats)){
      R1.vec <- c()
      R2.vec <- c()
      R3.vec <- c()
      R4.vec <- c()
      R5.vec <- c()
      for(kk in 1:parameter.num){

        Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        R1.vec <- c(R1.vec, Jaccard.OUR[1])
        R2.vec <- c(R2.vec, Jaccard.OUR[2])
        R3.vec <- c(R3.vec, Jaccard.OUR[3])
        R4.vec <- c(R4.vec, Jaccard.OUR[4])
        R5.vec <- c(R5.vec, Jaccard.OUR[5])
      }
      perf.mat3[k,] <- c(R1.vec, R2.vec, R3.vec, R4.vec, R5.vec)
    }
  }

  #colnames(perf.mat3) <- c(methods)

  library(ComplexHeatmap)
  ord1 <- order(rowMeans(perf.mat1), decreasing = T)
  ord2 <- order(rowMeans(perf.mat2), decreasing = T)
  ord3 <- order(rowMeans(perf.mat3), decreasing = T)

  col_fun = circlize::colorRamp2(c(0,0.5,1), c("grey95",'#FF99A3', "red"))
  col_fun = circlize::colorRamp2(seq(0,1, length.out=8), c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15"))
  #col_fun = circlize::colorRamp2(c(0,0.5,1), c("#4DADCF",'yellow', "red"))

  method.order <- c('RareQ','scCAD','CellSIUS','RaceID','GiniClust2','FiRE','EDGE','GapClust')
  cols <- c('#FF0099', '#C3EF00', '#007ED3','#FF9D1E','#7FD2FF', '#00C19B', '#894FC6', '#D55E00')
  names(cols) <- method.order

  perf.mat1.ord <- perf.mat1[ord1,]
  perf.mat2.ord <- perf.mat2[ord2,]
  perf.mat3.ord <- perf.mat3[ord3,]

  perf.mat1.out <- perf.mat1.ord
  perf.mat2.out <- perf.mat2.ord
  perf.mat3.out <- perf.mat3.ord
  colnames(perf.mat3.out) <- rep(c('R1', 'R2', 'R3', 'R4', 'R5'), each=parameter.num)

  write.csv(perf.mat1.out, file = paste0('Result/Sim-PBMC-1_', 'Reshuffle', '_Jaccard.csv'))
  write.csv(perf.mat2.out, file = paste0('Result/Sim-PBMC-2_', 'Reshuffle', '_Jaccard.csv'))
  write.csv(perf.mat3.out, file = paste0('Result/Sim-PBMC-3_', 'Reshuffle', '_Jaccard.csv'))

  row_ha1 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat1.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha2 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat2.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha3 = rowAnnotation('Types' = rep(c('R1','R2','R3','R4','R5'), each=parameter.num),  'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat3.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90','R2'='grey70','R3'='grey50','R4'='grey30','R5'='grey10')),
                          annotation_width=c(0.25,1))

  p1 <- Heatmap(t(perf.mat1.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha1, row_names_rot = 270)
  p2 <- Heatmap(t(perf.mat2.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha2, row_names_rot = 270)
  p3 <- Heatmap(t(perf.mat3.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', row_title_rot = 0, row_dend_width = unit(0, "mm"), right_annotation = row_ha3, show_row_names = T,
                row_split = factor(c(rep('R1', parameter.num),
                                     rep('R2', parameter.num),
                                     rep('R3', parameter.num),
                                     rep('R4', parameter.num),
                                     rep('R5', parameter.num)), levels=c('R1', 'R2', 'R3', 'R4', 'R5'), ordered = T))
  p.F1 <- (p1 %v% p2 %v% p3)
  return(p.F1)
}



ks <- c(1:30)


p.k <- Evaluate.parameter.Jaccard(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Result/Jaccard_index_PBMC_Sim123_shuffled_index.pdf'), width = 4, height = 6)
p.k
dev.off()










## Evaluate Jaccard index of top features
Evaluate.parameter.Jaccard.Top.Features <- function(parameter.list, parameter.symbol){

  parameter.num = length(parameter.list)

  ## For one rare cluster
  scenarios <- c("Rare1_Ordinary4","Rare1_Ordinary9","Rare5_Ordinary10") # ,
  repeats <- c(0:49)

  eval.rare.OUR <- function(scenario, rep_id, parameter.symbol, parameter.value){
    path <- paste0('site1/',scenario,'/',rep_id,'/R_input/')
    if(parameter.value==0){
      res <- readRDS(paste0(path, 'Top_marker_ref.RDS'))
    }else{
      res <- readRDS(paste0(path, 'Top_marker_',parameter.symbol, '_', parameter.value,'.RDS'))
    }

    res0 <- readRDS(paste0(path, 'Top_marker_ref.RDS'))

    if(scenario %in% c('Rare1_Ordinary4', 'Rare1_Ordinary9')){
      rare_type <- unique(res0$type)

      Jaccard.index = length(intersect(res$marker, res0$marker))/length(union(res$marker, res0$marker))

      return(c(Jaccard.index))

    }else{
      Jaccard.vec <- c()
      rare_types <- unique(res0$type)
      for(rare_type in rare_types){

        Jaccard.index = length(intersect(res$marker[res$type==rare_type], res0$marker[res0$type==rare_type]))/length(union(res$marker[res$type==rare_type], res0$marker[res0$type==rare_type]))
        Jaccard.vec <- c(Jaccard.vec, Jaccard.index)
      }
      return(Jaccard.vec)
    }
  }


  perf.mat1 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 1){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in c(1:parameter.num)){
        if(parameter.list[kk]==0){
          Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                       parameter.symbol=parameter.symbol, parameter.value=0)
        }else{
          Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                       parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        }

        acc.vec <- c(acc.vec, Jaccard.OUR)
      }
      perf.mat1[k,] <- acc.vec
    }
  }

  colnames(perf.mat1) <- parameter.list


  perf.mat2 <- matrix(NA, nrow=length(repeats), ncol=parameter.num)
  for(j in 2){
    for(k in 1:length(repeats)){
      acc.vec <- c()
      for(kk in 1:parameter.num){

        Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                     parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        acc.vec <- c(acc.vec, Jaccard.OUR)
      }
      perf.mat2[k,] <- acc.vec
    }
  }

  colnames(perf.mat2) <- parameter.list

  perf.mat3 <- matrix(NA, nrow=length(repeats), ncol=parameter.num*5)
  for(j in 3){
    for(k in 1:length(repeats)){
      R1.vec <- c()
      R2.vec <- c()
      R3.vec <- c()
      R4.vec <- c()
      R5.vec <- c()
      for(kk in 1:parameter.num){

        Jaccard.OUR <- eval.rare.OUR(scenario = scenarios[j], rep_id = repeats[k],
                                     parameter.symbol=parameter.symbol, parameter.value=parameter.list[kk])
        R1.vec <- c(R1.vec, Jaccard.OUR[1])
        R2.vec <- c(R2.vec, Jaccard.OUR[2])
        R3.vec <- c(R3.vec, Jaccard.OUR[3])
        R4.vec <- c(R4.vec, Jaccard.OUR[4])
        R5.vec <- c(R5.vec, Jaccard.OUR[5])
      }
      perf.mat3[k,] <- c(R1.vec, R2.vec, R3.vec, R4.vec, R5.vec)
    }
  }

  #colnames(perf.mat3) <- c(methods)

  library(ComplexHeatmap)
  ord1 <- order(rowMeans(perf.mat1), decreasing = T)
  ord2 <- order(rowMeans(perf.mat2), decreasing = T)
  ord3 <- order(rowMeans(perf.mat3), decreasing = T)

  col_fun = circlize::colorRamp2(c(0,0.5,1), c("grey95",'#FF99A3', "red"))
  col_fun = circlize::colorRamp2(seq(0,1, length.out=8), c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15"))
  #col_fun = circlize::colorRamp2(c(0,0.5,1), c("#4DADCF",'yellow', "red"))

  method.order <- c('RareQ','scCAD','CellSIUS','RaceID','GiniClust2','FiRE','EDGE','GapClust')
  cols <- c('#FF0099', '#C3EF00', '#007ED3','#FF9D1E','#7FD2FF', '#00C19B', '#894FC6', '#D55E00')
  names(cols) <- method.order

  perf.mat1.ord <- perf.mat1[ord1,]
  perf.mat2.ord <- perf.mat2[ord2,]
  perf.mat3.ord <- perf.mat3[ord3,]

  perf.mat1.out <- perf.mat1.ord
  perf.mat2.out <- perf.mat2.ord
  perf.mat3.out <- perf.mat3.ord
  colnames(perf.mat3.out) <- rep(c('R1', 'R2', 'R3', 'R4', 'R5'), each=parameter.num)


  row_ha1 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat1.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha2 = rowAnnotation('Types' = rep('R1', parameter.num), 'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat2.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90')), annotation_width=c(0.25,1), show_annotation_name = F)
  row_ha3 = rowAnnotation('Types' = rep(c('R1','R2','R3','R4','R5'), each=parameter.num),  'Mean(Jaccard index)' = anno_barplot(colMeans(perf.mat3.ord), ylim=c(0,1)),
                          col=list('Types' = c('R1'='grey90','R2'='grey70','R3'='grey50','R4'='grey30','R5'='grey10')),
                          annotation_width=c(0.25,1))

  p1 <- Heatmap(t(perf.mat1.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha1, row_names_rot = 270)
  p2 <- Heatmap(t(perf.mat2.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', show_row_dend = F, show_row_names = T, right_annotation = row_ha2, row_names_rot = 270)
  p3 <- Heatmap(t(perf.mat3.ord), cluster_columns = F, cluster_rows = F, border = T, col=col_fun, name='Jaccard index',
                row_names_side = 'left', row_title_rot = 0, row_dend_width = unit(0, "mm"), right_annotation = row_ha3, show_row_names = T,
                row_split = factor(c(rep('R1', parameter.num),
                                     rep('R2', parameter.num),
                                     rep('R3', parameter.num),
                                     rep('R4', parameter.num),
                                     rep('R5', parameter.num)), levels=c('R1', 'R2', 'R3', 'R4', 'R5'), ordered = T))
  p.F1 <- (p1 %v% p2 %v% p3)
  return(p.F1)
}



ks <- c(1:30)


p.k <- Evaluate.parameter.Jaccard.Top.Features(parameter.list = ks, parameter.symbol = 'shuffled_cols_seed')
pdf(file=paste0('Result/Jaccard_index_Top_Features_PBMC_Sim123_shuffled_index.pdf'), width = 4, height = 6)
p.k
dev.off()

