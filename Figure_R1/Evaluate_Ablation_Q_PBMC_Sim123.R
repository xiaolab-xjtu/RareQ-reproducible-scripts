library(anndata)

setwd('/home/rstudio/Projects/Rare_cell/data/PBMCs/')


# Evaluate rare cell detection accuracy via F1 score
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
      # ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      # res <- res[order(ids)]
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

  saveRDS(perf.mat1.out, file = paste0('Result/Sim-PBMC-1_', 'Ablation_Q', '_F1.RDS'))
  saveRDS(perf.mat2.out, file = paste0('Result/Sim-PBMC-2_', 'Ablation_Q', '_F1.RDS'))
  saveRDS(perf.mat3.out, file = paste0('Result/Sim-PBMC-3_', 'Ablation_Q', '_F1.RDS'))

}



ks <- c(0,1) # 0 denote the unshuffled result
p.k <- Evaluate.parameter.F1(parameter.list = ks, parameter.symbol = 'Ablation_Q')







## Evaluate Global clustering via NMI metric
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
      #ids <- readRDS(paste0(path, parameter.symbol, '_', parameter.value,'.RDS'))
      #res <- res[order(ids)]
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

  saveRDS(perf.mat1.out, file = paste0('Result/Sim-PBMC-1_', 'Ablation_Q', '_NMI.RDS'))
  saveRDS(perf.mat2.out, file = paste0('Result/Sim-PBMC-2_', 'Ablation_Q', '_NMI.RDS'))
  saveRDS(perf.mat3.out, file = paste0('Result/Sim-PBMC-3_', 'Ablation_Q', '_NMI.RDS'))

}

ks <- c(0,1) # 0 denote the unshuffled result
p.k.NMI <- Evaluate.parameter.NMI(parameter.list = ks, parameter.symbol = 'Ablation_Q')









