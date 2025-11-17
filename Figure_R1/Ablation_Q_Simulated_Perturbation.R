setwd('/home/rstudio/Projects/Rare_cell/data/Jurkat/')

library(RareQ)
library(Seurat)
library(ggplot)



source('FindRare_Q_Ablation.R', encoding='utf-8')



sc_object = readRDS('seu_obj.RDS')
DimPlot(sc_object, group.by = 'cluster')
sc_object$cell_type = factor(sc_object$cluster, levels=c(33,1125,1549), labels=c('293T_1', '293T_2', 'Jurkat'))

set.seed(2025)
sample.id <- sample(which(sc_object$cell_type == '293T_1'), size=15)
top.genes = sort(rowMeans(sc_object@assays$RNA@data[,sample.id]), decreasing = T)[1:200]
count = sc_object@assays$RNA@counts

labels = as.character(sc_object$cell_type)
labels_perturb = labels
labels_perturb[sample.id] = '293T_perturb'


DEG.nums = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
for(DEG.num in DEG.nums){
  for(rep in 1:10){
    set.seed(2*rep - rep^2)
    gene.sample <- sample(names(top.genes), DEG.num)
    count.change = count
    count.change[gene.sample,sample.id] = count[gene.sample,sample.id] * 3

    sc_object_perturb <- CreateSeuratObject(count=count.change, project = "sc_object_perturb", min.cells = 3)
    sc_object_perturb <- NormalizeData(sc_object_perturb) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData()
    sc_object_perturb <- RunPCA(sc_object_perturb, features = VariableFeatures(object = sc_object_perturb))
    sc_object_perturb <- FindNeighbors(object = sc_object_perturb,
                                        k.param = 20,
                                        compute.SNN = F,
                                        prune.SNN = 0,
                                        reduction = "pca",
                                        dims = 1:50,
                                        force.recalc = F, return.neighbor = T)
    cluster_ablation_Q = FindRare_Ablation_Q(sc_object = sc_object_perturb)
    cluster_standard = FindRare(sc_object = sc_object_perturb)

    saveRDS(cluster_standard, file = paste0('Perturb/OUR_result_', DEG.num,'_DEGs_seed_', rep,'.RDS'))
    saveRDS(cluster_ablation_Q, file = paste0('Perturb/OUR_result_', DEG.num,'_DEGs_seed_',rep, '_Q_Ablation.RDS'))
  }
}





eval.rare <- function(labs, res, method){
  type.cnt =  table(labs)
  rare_types <- c('293T_perturb','293T_2','Jurkat')

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
                       Recall=Recall.vec)
  res.df$Method = method
  return(res.df)
}


res.df <- data.frame(Type=NA,
                     F1=NA,
                     Precision=NA,
                     Recall=NA,
                     Method=NA,
                     Num=NA,
                     Rep=NA)
res.df.perturb <- data.frame(Type=NA,
                     F1=NA,
                     Precision=NA,
                     Recall=NA,
                     Method=NA,
                     Num=NA,
                     Rep=NA)

for(DEG.num in DEG.nums){
  for(rep in 1:10){

    cluster_standard = readRDS(paste0('Perturb/OUR_result_', DEG.num,'_DEGs_seed_', rep,'.RDS'))
    cluster_ablation_Q = readRDS(paste0('Perturb/OUR_result_', DEG.num,'_DEGs_seed_',rep, '_Q_Ablation.RDS'))

    acc1 = eval.rare(labs = labels_perturb, res = cluster_standard, method = 'Standard')
    acc2 = eval.rare(labs = labels_perturb, res = cluster_ablation_Q, method = 'Q_Ablation')
    acc1$Num = DEG.num
    acc1$Rep = rep
    res.df = rbind(res.df, acc1)

    acc2$Num = DEG.num
    acc2$Rep = rep
    res.df.perturb = rbind(res.df.perturb, acc2)

  }
}

res.df.perturb = res.df.perturb[!is.na(res.df.perturb$Type),]
res.df.perturb$Precision[is.na(res.df.perturb$Precision)] = 0
res.df = res.df[!is.na(res.df$Type),]
res.df$Precision[is.na(res.df$Precision)] = 0

# saveRDS(res.df, file = 'Perturb/Standard_res.RDS')
# saveRDS(res.df.perturb, file = 'Perturb/Q_Ablation_res.RDS')


res.df.ave = data.frame(Num=as.integer(names(tapply(res.df$F1, res.df$Num, mean))),
                        F1=tapply(res.df$F1, res.df$Num, mean),
                        Precision=tapply(res.df$Precision, res.df$Num, mean),
                        Recall=tapply(res.df$Recall, res.df$Num, mean),
                        Method='Standard')

res.df.perturb.ave = data.frame(Num=as.integer(names(tapply(res.df.perturb$F1, res.df.perturb$Num, mean))),
                                F1=tapply(res.df.perturb$F1, res.df.perturb$Num, mean),
                                Precision=tapply(res.df.perturb$Precision, res.df.perturb$Num, mean),
                                Recall=tapply(res.df.perturb$Recall, res.df.perturb$Num, mean),
                                Method='Q_Ablation')

res.ave.comb = rbind(res.df.ave, res.df.perturb.ave)
res.ave.comb$Method = factor(res.ave.comb$Method, levels = c('Standard', 'Q_Ablation'), ordered = T)

p.F1 = ggplot(data=res.ave.comb, aes(x=Num, y=F1, color=factor(Method), shape=factor(Method))) + geom_point() + geom_line() + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#FF0099', '#C3EF00')) + labs(x='Num of perturbed genes') +
  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))
p.Precision = ggplot(data=res.ave.comb, aes(x=Num, y=Precision, color=factor(Method), shape=factor(Method))) + geom_point() + geom_line() + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#FF0099', '#C3EF00')) + labs(x='Num of perturbed genes') +
  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))
p.Recall = ggplot(data=res.ave.comb, aes(x=Num, y=Recall, color=factor(Method), shape=factor(Method))) + geom_point() + geom_line() + theme_bw() +
  theme(panel.border = element_rect(fill=NA), axis.text = element_text(colour = 'black'), legend.title = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c('#FF0099', '#C3EF00')) + labs(x='Num of perturbed genes') +
  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100))

p.acc = ggpubr::ggarrange(p.F1, p.Precision, p.Recall, ncol=3, align='hv', common.legend = T, legend = 'right')
ggsave(p.acc, filename = 'Perturb/Acuracy_comparsion_perturbation.pdf', width = 8, height = 2)


