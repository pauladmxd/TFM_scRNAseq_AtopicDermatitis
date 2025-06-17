# Databricks notebook source
# MAGIC %md
# MAGIC #Find markers with seurat
# MAGIC ##AD/PN – Alkon et al, 2023
# MAGIC - Does not include non-lesional samples.
# MAGIC - 5 AD, 7 PN, 1 AP (atopic prurigo) and 4 healthy control​​
# MAGIC
# MAGIC Here will be performed the analysis to find all markers of each relevant cell type in Alkon et al, 2023 dataset.
# MAGIC
# MAGIC ###Most relevant cell types: 
# MAGIC T-cells (TC), Fibroblasts, Keratinocytes (KC),  Monocytes, Macrophages, Dendritic cells, Natural killers, Treg and MastC
# MAGIC
# MAGIC ###Constrast: 
# MAGIC   - Lesional vs Healthy control (LvsHC)

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(SingleCellExperiment)

# COMMAND ----------

.libPaths(c("/dbfs/home/pdelgadom@almirall.com/my_r_packages/tfm_paula_4", .libPaths()))
library(MAST)

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat_v2", .libPaths()))
library(Seurat)

# COMMAND ----------

#Read seurat object
alkon <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")

# COMMAND ----------

head(alkon@meta.data, 5)

# COMMAND ----------

unique(alkon$h_celltype_v4)

# COMMAND ----------

# Check if all conditions have at least 50 cells for all cell types
# Print result for each cell type
cell_types <- unique(alkon$h_celltype_v4)
for(cell_type in cell_types) {
  sub <- subset(alkon, subset = h_celltype_v4 == cell_type)
  if(all((table(sub@meta.data$Condition)) >= 50)) {
    print(paste("All conditions for", cell_type, "have at least 50 cells."))
  } else {
    print(paste("Not all conditions for", cell_type, "have at least 50 cells."))
  }
}

# COMMAND ----------

# MAGIC %md
# MAGIC To check that the analysis is significant to be performed, I first check that each condition has at least 50 cells (decided as a minimum)

# COMMAND ----------

# MAGIC %md
# MAGIC ## T cell markers

# COMMAND ----------

#First it is performed a subset with only T cells
subtcell <- subset(alkon, subset = h_celltype_v4 == "TC")
unique(subtcell$h_celltype_v4) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers tcell

# COMMAND ----------

Idents(subtcell) <- subtcell$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

?FindAllMarkers

# COMMAND ----------

subtcell.markers <- FindAllMarkers(subtcell, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

alkon_LvsHC_tcell_markers <- subtcell.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_LvsHC_tcell_markers)

# COMMAND ----------

# %sh
# mkdir /dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method

# COMMAND ----------

#Save results
#write.xlsx(alkon_LvsHC_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_LvsHC_tcell_markers.xlsx")
write.xlsx(subtcell.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

#First it is performed a subset with only Fibroblasts
subfb <- subset(alkon, subset = h_celltype_v4 == "Fibroblasts")
unique(subfb$h_celltype_v4) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers fibroblasts

# COMMAND ----------

Idents(subfb) <- subfb$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

subfb.markers <- FindAllMarkers(subfb, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

alkon_fb_LvsHC_markers <- subfb.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_fb_LvsHC_markers)

# COMMAND ----------

# Save results
#write.xlsx(alkon_fb_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_markers.xlsx")
write.xlsx(subfb.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

#First it is performed a subset with only KC
subkc <- subset(alkon, subset = h_celltype_v4 == "KC")
unique(subkc$h_celltype_v4) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers KC

# COMMAND ----------

Idents(subkc) <- subkc$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

subkc.markers <- FindAllMarkers(subkc, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

alkon_kc_LvsHC_markers <- subkc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_kc_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(alkon_kc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_markers.xlsx")
write.xlsx(subkc.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only macrophages
submacro <- subset(alkon, subset = h_celltype_v4 == "Macro")
unique(submacro$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(submacro@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}

# COMMAND ----------

table(submacro@meta.data$Condition, submacro@meta.data$h_celltype_v4)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers macrophages

# COMMAND ----------

Idents(submacro) <- submacro$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

macro_LvsHC.markers <- FindAllMarkers(submacro, only.pos=TRUE,  test.use="MAST")

# COMMAND ----------

alkon_macro_LvsHC_markers <- macro_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_macro_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_macro_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_markers.xlsx")
write.xlsx( macro_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

#First it is performed a subset with only T reg
subtreg <- subset(alkon, subset = h_celltype_v4 == "Treg")
unique(subtreg$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(subtreg@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(subtreg@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC As it is just 5 below I continue performing the analysis but taking it into account.

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers treg

# COMMAND ----------

Idents(subtreg) <- subtreg$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

subtreg.markers <- FindAllMarkers(subtreg, only.pos=TRUE,  test.use="MAST")

# COMMAND ----------

alkon_LvsHC_treg_markers <- subtreg.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_LvsHC_treg_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_LvsHC_treg_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_treg_LvsHC_markers.xlsx")
write.xlsx(subtreg.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC #Volcano Plots

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(dplyr)
library(openxlsx)

# COMMAND ----------

.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(EnhancedVolcano)

# COMMAND ----------

volcano_generator<- function(resultsDE){
resultsDE<- as.data.frame(resultsDE)

#Create anotations for volcano plot
resultsDE0 <- resultsDE
rownames(resultsDE0) <- resultsDE0$gene

top10_genes <- resultsDE0 %>% 
filter(avg_log2FC > 2 & p_val_adj < 0.05 ) %>%
arrange(p_val_adj) %>% top_n(10, -p_val_adj)

bottom10_genes <-  resultsDE0 %>% 
filter(avg_log2FC < -1 & p_val_adj < 0.05 ) %>%
arrange(p_val_adj) %>% top_n(10, -p_val_adj)

#Plot Volcano
volcano <- (EnhancedVolcano(resultsDE0,
lab = rownames(resultsDE0),
x = 'avg_log2FC',
y = 'p_val_adj',
pCutoff = 0.05,
selectLab = c(top10_genes$gene, bottom10_genes$gene),
labSize = 5,
drawConnectors = TRUE,
widthConnectors = 0.5,
colConnectors = 'black'))
volcano

}

# COMMAND ----------

subtreg.markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

  subtreg.markers <- subtreg.markers %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(subtreg.markers)

# COMMAND ----------

subkc.markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

subkc.markers.seurat <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

  subkc.markers <- subkc.markers %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

  subkc.markers.seurat <- subkc.markers.seurat %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(subkc.markers.seurat)

# COMMAND ----------

volcano_generator(subkc.markers)

# COMMAND ----------

subtcell.markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/alkon_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

  subtcell.markers <- subtcell.markers %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(subtcell.markers)
