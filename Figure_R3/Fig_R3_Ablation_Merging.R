library(ggplot2)
library(dplyr)
library(tidyr)

# F1 score----------------------------
## sim1&2-------------------------
df <- readRDS("Sim-PBMC-1_Ablation_Q_F1.RDS") %>% as.data.frame()
colnames(df) <- c("Standard", "Ablation (Merging)")

df_long <- df %>%
  select(Standard, `Ablation (Merging)`) %>% 
  pivot_longer(
    cols = everything(),                 
    names_to = "Condition",              
    values_to = "F1"                    
  )

df_long <- df_long %>%
  mutate(Condition = factor(Condition, levels = c("Standard", "Ablation (Merging)")))

df_cdf <- df_long %>%
  group_by(Condition) %>%
  arrange(F1) %>%
  mutate(cumulative = seq_along(F1) / n())


theme1 <- theme(axis.title=element_text(size=17,color = 'black'), 
                   axis.text=element_text(size=14,color = 'black',angle = 0, lineheight = 0.6), 
                   legend.title = element_text(size=16),   
                   legend.text = element_text(size=14),   
                   legend.position="right",  
                   axis.line = element_line(linewidth=0, colour = 'black'), 
                   panel.background = element_rect(linewidth=0.6,color='black'),  
                   legend.key = element_blank(),
                   panel.grid.major = element_line(color = "gray80", linewidth = 0.3)
)

p1 <- ggplot(df_cdf, aes(x = cumulative, y = F1, color = Condition)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = expression(F[1]~"score"),
    title = "Sim-PBMC-1",
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p1


## sim3--------------------------
df <- readRDS("Sim-PBMC-3_Ablation_Q_F1.RDS") %>% as.data.frame()
colnames(df) <- paste0(colnames(df), "_", c("1","2","1","2","1","2","1","2","1","2"))
groups <- unique(gsub("_.*", "", colnames(df)))

for (Ri in groups){
  df_i <- df[, grep(paste0("^", Ri, "_"), colnames(df))]
  colnames(df_i) <- c("Standard", "Ablation (Merging)")
  
  df_long <- df_i %>%
    select(Standard, `Ablation (Merging)`) %>% 
    pivot_longer(
      cols = everything(),                 
      names_to = "Condition",              
      values_to = "F1"                     
    )
  df_long <- df_long %>%
    mutate(Condition = factor(Condition, levels = c("Standard", "Ablation (Merging)")))
  
  df_cdf <- df_long %>%
    group_by(Condition) %>%
    arrange(F1) %>%
    mutate(cumulative = seq_along(F1) / n())
  
  p <- ggplot(df_cdf, aes(x = cumulative, y = F1, color = Condition)) +
    geom_step(size = 0.7) +
    theme_classic(base_size = 14) +
    labs(
      x = "Cumulative proportion",
      y = expression(F[1]~"score"),
      title = paste0("Sim-PBMC-3 ", Ri),
      color = ""
    ) +
    scale_y_continuous(
      limits = c(0, 1),          
      breaks = seq(0, 1, 0.25)   
    ) +
    scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
    theme1
  
  out_path <- paste0("~/Rare_Q/visualisasion/fig/sim3_", Ri, "_F1.pdf")
  ggsave(out_path, p, width = 5, height = 3)
}


# 20 real datasets--------------------
df <- readRDS("F1_Precision_Recall_20_datasets_Ablation_Q.RDS") %>% as.data.frame()
df <- df %>% 
  mutate(para.value = gsub("Merging_ablation","Ablation (Merging)",para.value))
df <- df %>%
  mutate(para.value = factor(para.value, levels = c("Standard", "Ablation (Merging)")))

## F1------
df_cdf <- df %>%
  select(F1,para.value) %>% 
  group_by(para.value) %>%
  arrange(F1) %>%
  mutate(cumulative = seq_along(F1) / n())

p_f1 <- ggplot(df_cdf, aes(x = cumulative, y = F1, color = para.value)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = "",
    title = expression(F[1]~"score"),
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p_f1


## Recall--------------
df_cdf <- df %>%
  select(Recall,para.value) %>% 
  group_by(para.value) %>%
  arrange(Recall) %>%
  mutate(cumulative = seq_along(Recall) / n())

p_Recall <- ggplot(df_cdf, aes(x = cumulative, y = Recall, color = para.value)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = "",
    title = "Recall",
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p_Recall

## Precision-----------------
df_cdf <- df %>%
  select(Precision,para.value) %>% 
  group_by(para.value) %>%
  arrange(Precision) %>%
  mutate(cumulative = seq_along(Precision) / n())

p_Precision <- ggplot(df_cdf, aes(x = cumulative, y = Precision, color = para.value)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = "",
    title = "Precision",
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p_Precision



#NMI-----------------------------------------------
## sim1&2&3-------------------------
df <- readRDS("Sim-PBMC-1_Ablation_Q_NMI.RDS") %>% as.data.frame()
colnames(df) <- c("Standard", "Ablation (Merging)")

df_long <- df %>%
  select(Standard, `Ablation (Merging)`) %>% 
  pivot_longer(
    cols = everything(),                
    names_to = "Condition",             
    values_to = "NMI"                   
  )

df_long <- df_long %>%
  mutate(Condition = factor(Condition, levels = c("Standard", "Ablation (Merging)")))

df_cdf <- df_long %>%
  group_by(Condition) %>%
  arrange(NMI) %>%
  mutate(cumulative = seq_along(NMI) / n())

p2 <- ggplot(df_cdf, aes(x = cumulative, y = NMI, color = Condition)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = "NMI",
    title = "Sim-PBMC-1",
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p2


## 20 real datasets--------------------
df <- readRDS("NMI_20_datasets_Ablation_Q.RDS") %>% as.data.frame()
df <- df %>% 
  mutate(para.value = gsub("Merging_ablation","Ablation (Merging)",para.value))
df <- df %>%
  mutate(para.value = factor(para.value, levels = c("Standard", "Ablation (Merging)")))

df_cdf <- df %>%
  select(NMI,para.value) %>% 
  group_by(para.value) %>%
  arrange(NMI) %>%
  mutate(cumulative = seq_along(NMI) / n())

p_NMI <- ggplot(df_cdf, aes(x = cumulative, y = NMI, color = para.value)) +
  geom_step(size = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cumulative proportion",
    y = "",
    title = "NMI",
    color = ""
  ) +
  scale_color_manual(values = c("Standard" = "#C72228", "Ablation (Merging)" = "#3F84AA")) +
  theme1
p_NMI