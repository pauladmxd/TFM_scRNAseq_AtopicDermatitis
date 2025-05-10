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
library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

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
# MAGIC ## T cell markers

# COMMAND ----------

#First it is performed a subset with only T cells
subtcell <- subset(alkon, subset = h_celltype_v4 == "TC")
unique(subtcell$h_celltype_v4) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC To check that the analysis is significant to be performed, I first check that each condition has at least 50 cells (decided as a minimum)

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers tcell

# COMMAND ----------

Idents(subtcell) <- subtcell$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

subtcell.markers <- FindAllMarkers(subtcell, only.pos=TRUE)

# COMMAND ----------

alkon_LvsHC_tcell_markers <- subtcell.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_LvsHC_tcell_markers)

# COMMAND ----------

#Save results
#write.xlsx(alkon_LvsHC_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_LvsHC_tcell_markers.xlsx")
write.xlsx(subtcell.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_LvsHC_tcell_allmarkers.xlsx")

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

subfb.markers <- FindAllMarkers(subfb, only.pos=TRUE)

# COMMAND ----------

alkon_fb_LvsHC_markers <- subfb.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_fb_LvsHC_markers)

# COMMAND ----------

# Save results
#write.xlsx(alkon_fb_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_markers.xlsx")
write.xlsx(subfb.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_allmarkers.xlsx")

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

subkc.markers <- FindAllMarkers(subkc, only.pos=TRUE)

# COMMAND ----------

alkon_kc_LvsHC_markers <- subkc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_kc_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(alkon_kc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_markers.xlsx")
write.xlsx(subkc.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only monocytes
submono <- subset(alkon, subset = h_celltype_v4 == "Mono")
unique(submono$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(submono@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(submono@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC In this case, HC condition of Monocytes has only 18 cells.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers monocytes

# COMMAND ----------

# Idents(submono) <- submono$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

# submono_LvsHC.markers <- FindAllMarkers(submono, only.pos=TRUE)

# COMMAND ----------

# alkon_mono_LvsHC_markers <- submono_LvsHC.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# display(alkon_mono_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_mono_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_mono_LvsHC_markers.xlsx")

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

macro_LvsHC.markers <- FindAllMarkers(submacro, only.pos=TRUE)

# COMMAND ----------

alkon_macro_LvsHC_markers <- macro_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_macro_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_macro_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_markers.xlsx")
write.xlsx( macro_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only DC
subdc <- subset(alkon, subset = h_celltype_v4 == "DC")
unique(subdc$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(subdc@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(subdc@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC In this case we only have 8 Dendritic cells in HC condition

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers DC

# COMMAND ----------

# Idents(subdc) <- subdc$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

# dc_LvsHC.markers <- FindAllMarkers(subdc, only.pos=TRUE)

# COMMAND ----------

# alkon_dc_LvsHC_markers <- dc_LvsHC.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# display(alkon_dc_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_dc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_dc_LvsHC_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Natural killers

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only NK
subnk <- subset(alkon, subset = h_celltype_v4 == "NK")
unique(subnk$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(subnk@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(subnk@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC There are only 2 cells and only in the condition AD.

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers NK

# COMMAND ----------

# Idents(subnk) <- subnk$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

# nk_LvsHC.markers <- FindAllMarkers(subnk, only.pos=TRUE)

# COMMAND ----------

# alkon_nk_LvsHC_markers <- nk_LvsHC.markers %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
# display(alkon_nk_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_nk_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_nk_LvsHC_markers.xlsx")

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

subtreg.markers <- FindAllMarkers(subtreg, only.pos=TRUE)

# COMMAND ----------

alkon_LvsHC_treg_markers <- subtreg.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(alkon_LvsHC_treg_markers)

# COMMAND ----------

#Save results
# write.xlsx(alkon_LvsHC_treg_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_treg_LvsHC_markers.xlsx")
write.xlsx(subtreg.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

#First it is performed a subset with only monocytes
submast <- subset(alkon, subset = h_celltype_v4 == "MastC")
unique(submast$h_celltype_v4) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(submast@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(submast@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Majority in AD patients, but not enough cells to do DEA
