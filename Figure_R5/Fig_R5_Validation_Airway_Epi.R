# object preprocessing--------------------------------
object <- NormalizeData(object)
all.genes <- rownames(object@assays$Spatial@data)
object <- ScaleData(object, features = all.genes, assay = "Spatial")

### C3/C6-------
thresholds <- c('KRT5' = 0.1, "TP63" = 0.1, 'NGFR' = 0.1, "TOP2A" = 0.1, "MKI67" = 0.1, "H2AFZ" = 0.1)

for (gene in names(thresholds)) {
  status_colname <- paste0(gene, "_status")
  object@meta.data[[status_colname]] <- ifelse(object@assays$Spatial@scale.data[gene, ] > thresholds[gene],
                                               "High", "Low")
}

object$C3_C6 <- "Other"

object$C3_C6[object$TOP2A_status == "Low" & object$MKI67_status == "High" & 
               object$H2AFZ_status == "High" & object$KRT5_status == "High" & 
               object$TP63_status == "High" & object$NGFR_status == "High"] <- "C3"

object$C3_C6[object$TOP2A_status == "High" & object$MKI67_status == "High" & 
               object$H2AFZ_status == "High" & object$KRT5_status == "High" & 
               object$TP63_status == "High" & object$NGFR_status == "High"] <- "C6"

SpatialDimPlot(object,
               group.by = "C3_C6",
               pt.size.factor = 2) +
  scale_fill_manual(values = c("Other" = "#532C8A", "C3" = "#00F5FF", "C6" = "#FACB12")) +
  theme(legend.position = "right")+
  ggtitle("")


# C4------------
thresholds <- c('CCNO' = 0.1, 'FOXJ1' = 0.1, "LRRC23" = 0.1, "MEIG1" = 0.1, "TPPP3" = 0.1)

for (gene in names(thresholds)) {
  status_colname <- paste0(gene, "_status")
  object@meta.data[[status_colname]] <- ifelse(object@assays$Spatial@scale.data[gene, ] > thresholds[gene],
                                               "High", "Low")
}

object$C4 <- "Other"
object$C4[object$CCNO_status == "High" & object$FOXJ1_status == "High" & 
            object$LRRC23_status == "High" & object$MEIG1_status == "High" &
            object$TPPP3_status == "High"] <- "C4"

SpatialDimPlot(object,
               group.by = "C4",
               pt.size.factor = 2) +
  scale_fill_manual(values = c("Other" = "#532C8A", "C4" = "#39FF14")) +
  theme(legend.position = "right")+
  ggtitle("")


# C10-----------
thresholds <- c('ISG15' = 0.1, 'SCGB1A1' = 0.1,"IFIT3" = 0.1,"FOXJ1" = 0.1,'LRRC23' = 0.1,'TUBA1A' = 0.1)

for (gene in names(thresholds)) {
  status_colname <- paste0(gene, "_status")
  object@meta.data[[status_colname]] <- ifelse(object@assays$Spatial@scale.data[gene, ] > thresholds[gene],
                                               "High", "Low")
}

object$C10 <- "Other"

object$C10[object$FOXJ1_status == "High" & object$LRRC23_status == "High" &
             object$TUBA1A_status == "High" & object$SCGB1A1_status == "High" &
             object$ISG15_status == "Low" & object$IFIT3_status == "Low"] <- "C10"

SpatialDimPlot(object,
               group.by = "C10",
               pt.size.factor = 2) +
  scale_fill_manual(values = c("Other" = "#532C8A", "C10" = "#f77b59")) +
  theme(legend.position = "right")+
  ggtitle("")