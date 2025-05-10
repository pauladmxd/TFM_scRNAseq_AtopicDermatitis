# Databricks notebook source
# MAGIC %md
# MAGIC #Find markers
# MAGIC ## Liu et al, 2021
# MAGIC - Includes lesional and non-lesional samples
# MAGIC - It only has immune cells
# MAGIC - 4 AD patients, 3 psoriasis and 5 healthy controls.
# MAGIC
# MAGIC Here will be performed the analysis to find all markers of each relevant cell type in Liu et al, 2021 dataset.
# MAGIC
# MAGIC ###Most relevant cell types: 
# MAGIC T-cells (TC), T-reg, Innate linfoid cells (ILC), Monocytes (Mono), Macrophagues (Macro), Dendritic cells (DC), Natural killers (NK), MastC
# MAGIC
# MAGIC ###Constrasts: 
# MAGIC   - Lesional vs Healthy control (LvsHC)
# MAGIC   - Non lesional vs Healthy control (NLvsHC)
# MAGIC   - Lesional vs Non-lesional (LvsNL)

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

#Read seurat object
liu <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/LIU_PROCESSED_TFM.rds")

# COMMAND ----------

head(liu@meta.data, 5)

# COMMAND ----------

# Check if all conditions have at least 50 cells for all cell types
# Print result for each cell type
cell_types <- unique(liu$h_celltype)
for(cell_type in cell_types) {
  sub <- subset(liu, subset = h_celltype == cell_type)
  if(all((table(sub@meta.data$Condition)) >= 50)) {
    print(paste("All conditions for", cell_type, "have at least 50 cells."))
  } else {
    print(paste("Not all conditions for", cell_type, "have at least 50 cells."))
  }
}

# COMMAND ----------

# MAGIC %md
# MAGIC I create a new column called Condition to keep the information of healthy or AD and lesional or non lesional all in one

# COMMAND ----------

# Create a new column 'Condition' with default value 'healthy'
liu$Condition <- "healthy"

# Update 'Condition' based on the 'Status' and 'Site' columns
liu$Condition[liu$Status == "Eczema" & liu$Site == "lesion"] <- "lesional"
liu$Condition[liu$Status == "Eczema" & liu$Site == "non_lesion"] <- "non lesional"

unique(liu$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC ## T cell markers

# COMMAND ----------

#First it is performed a subset with only T cells
subtcell <- subset(liu, subset = h_celltype == "TC")
unique(subtcell$h_celltype) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC

# COMMAND ----------

#Then, it is performed a subset with only conditions for this constrast
tcell_LvsHC <- subset(subtcell, Condition %in% c("lesional", "healthy"))
unique(tcell_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers tcell

# COMMAND ----------

Idents(tcell_LvsHC) <- tcell_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

tcell_LvsHC.markers <- FindAllMarkers(tcell_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_LvsHC_tcell_markers <- tcell_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_LvsHC_tcell_markers)

# COMMAND ----------

#Save results
# write.xlsx(liu_LvsHC_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_LvsHC_tcell_markers.xlsx")
write.xlsx(tcell_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only ILC
subILC <- subset(liu, subset = h_celltype == "ILC")
unique(subILC$h_celltype) #to check

# COMMAND ----------

ILC_LvsHC <- subset(subILC, Condition %in% c("lesional", "healthy"))
unique(ILC_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers ILC

# COMMAND ----------

Idents(ILC_LvsHC) <- ILC_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

ILC_LvsHC.markers <- FindAllMarkers(ILC_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_ILC_LvsHC_markers <- ILC_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5  & p_val_adj < 0.05)
display(liu_ILC_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(liu_ILC_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_ILC_LvsHC_markers.xlsx")
write.xlsx(ILC_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_ILC_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##T reg

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only ILC
subtreg <- subset(liu, subset = h_celltype == "Treg")
unique(subtreg$h_celltype) #to check

# COMMAND ----------

treg_LvsHC <- subset(subtreg, Condition %in% c("lesional", "healthy"))
unique(treg_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers T reg

# COMMAND ----------

Idents(treg_LvsHC) <- treg_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

treg_LvsHC.markers <- FindAllMarkers(treg_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_treg_LvsHC_markers <- treg_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_treg_LvsHC_markers)

# COMMAND ----------

# Save results
#write.xlsx(liu_treg_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_treg_LvsHC_markers.xlsx")
write.xlsx(treg_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only monocytes
submono <- subset(liu, subset = h_celltype == "Mono")
unique(submono$h_celltype) #to check

# COMMAND ----------

mono_LvsHC <- subset(submono, Condition %in% c("lesional", "healthy"))
unique(mono_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers monocytes

# COMMAND ----------

Idents(mono_LvsHC) <- mono_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

mono_LvsHC.markers <- FindAllMarkers(mono_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_mono_LvsHC_markers <- mono_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_mono_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(liu_mono_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mono_LvsHC_markers.xlsx")
write.xlsx(mono_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

# First it is performed a subset with only macrophages
submacro <- subset(liu, subset = h_celltype == "Macro")
unique(submacro$h_celltype) #to check

# COMMAND ----------

macro_LvsHC <- subset(submacro, Condition %in% c("lesional", "healthy"))
unique(macro_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers macrophages

# COMMAND ----------

Idents(macro_LvsHC) <- macro_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

macro_LvsHC.markers <- FindAllMarkers(macro_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_macro_LvsHC_markers <- macro_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_macro_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(liu_macro_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_macro_LvsHC_markers.xlsx")
write.xlsx(macro_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only DC
subdc <- subset(liu, subset = h_celltype == "DC")
unique(subdc$h_celltype) #to check

# COMMAND ----------

dc_LvsHC <- subset(subdc, Condition %in% c("lesional", "healthy"))
unique(dc_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers DC

# COMMAND ----------

Idents(dc_LvsHC) <- dc_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

dc_LvsHC.markers <- FindAllMarkers(dc_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_dc_LvsHC_markers <- dc_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_dc_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(liu_dc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_dc_LvsHC_markers.xlsx")
write.xlsx(dc_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Natural killers

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only NK
subnk <- subset(liu, subset = h_celltype == "NK")
unique(subnk$h_celltype) #to check

# COMMAND ----------

nk_LvsHC <- subset(subnk, Condition %in% c("lesional", "healthy"))
unique(nk_LvsHC$Condition) #to check

# COMMAND ----------

table(nk_LvsHC@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough cells but I did the analysis

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers NK

# COMMAND ----------

Idents(nk_LvsHC) <- nk_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

nk_LvsHC.markers <- FindAllMarkers(nk_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_nk_LvsHC_markers <- nk_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_nk_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(liu_nk_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_nk_LvsHC_markers.xlsx")
write.xlsx(nk_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_nk_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

# MAGIC %md
# MAGIC ###LvsHC

# COMMAND ----------

#First it is performed a subset with only mast cells
submast <- subset(liu, subset = h_celltype == "MastC")
unique(submast$h_celltype) #to check

# COMMAND ----------

mast_LvsHC <- subset(submast, Condition %in% c("lesional", "healthy"))
unique(mast_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC ### Find markers mastC

# COMMAND ----------

Idents(mast_LvsHC) <- mast_LvsHC$Condition

# COMMAND ----------

mast_LvsHC.markers <- FindAllMarkers(mast_LvsHC, only.pos=TRUE)

# COMMAND ----------

liu_mast_LvsHC_markers <- mast_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_mast_LvsHC_markers)

# COMMAND ----------

# Save results
# write.xlsx(liu_mast_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mast_LvsHC_markers.xlsx")
write.xlsx(mast_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mast_LvsHC_allmarkers.xlsx")
