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
  perf.df$para.value = factor(perf.df$para.value, levels=c(0,1), labels = c('Standard', 'Merging_ablation'), ordered = T)

  saveRDS(perf.df, file = 'Figure/20_datasets_Ablation_Merging_F1.RDS')

}


ks <- c(0, 1) # 0 denote the unshuffled result


p.k <- Evaluate.parameter(parameter.list = ks, parameter.symbol = 'Ablation_Merging')







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
  perf.df$para.value = factor(perf.df$para.value, levels=c(0,1), labels = c('Standard', 'Merging_ablation'), ordered = T)

  saveRDS(perf.df, file = 'Figure/20_datasets_Ablation_Merging_NMI.RDS')
}

ks <- c(0, 1) # 0 denote the unshuffled result

p.k.NMI <- Evaluate.parameter.NMI(parameter.list = ks, parameter.symbol = 'Ablation_Merging')









