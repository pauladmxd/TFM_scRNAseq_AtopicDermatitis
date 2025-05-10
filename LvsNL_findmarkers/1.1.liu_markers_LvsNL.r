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
# MAGIC   - Lesional vs Non-lesional (LvsNL)  <--- Here this one will be done

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
# MAGIC ## T cell markers

# COMMAND ----------

#First it is performed a subset with only T cells
subtcell <- subset(liu, subset = h_celltype == "TC")
unique(subtcell$h_celltype) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL

# COMMAND ----------

#Then, it is performed a subset with only conditions for this constrast
tcell_LvsNL <- subset(subtcell, Condition %in% c("lesional", "non lesional"))
unique(tcell_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers tcell

# COMMAND ----------

Idents(tcell_LvsNL) <- tcell_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

tcell_LvsNL.markers <- FindAllMarkers(tcell_LvsNL, only.pos=TRUE)

# COMMAND ----------

tcell_LvsNL.markers <- tcell_LvsHC.markers

# COMMAND ----------

liu_LvsNL_tcell_markers <- tcell_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_LvsNL_tcell_markers)

# COMMAND ----------

#Save results
write.xlsx(liu_LvsNL_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_LvsNL_tcell_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only ILC
subILC <- subset(liu, subset = h_celltype == "ILC")
unique(subILC$h_celltype) #to check

# COMMAND ----------

ILC_LvsNL <- subset(subILC, Condition %in% c("lesional", "non lesional"))
unique(ILC_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers ILC

# COMMAND ----------

Idents(ILC_LvsNL) <- ILC_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

ILC_LvsNL.markers <- FindAllMarkers(ILC_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_ILC_LvsNL_markers <- ILC_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5  & p_val_adj < 0.05)
display(liu_ILC_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(liu_ILC_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_ILC_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##T reg

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only ILC
subtreg <- subset(liu, subset = h_celltype == "Treg")
unique(subtreg$h_celltype) #to check

# COMMAND ----------

treg_LvsNL <- subset(subtreg, Condition %in% c("lesional", "non lesional"))
unique(treg_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers T reg

# COMMAND ----------

Idents(treg_LvsNL) <- treg_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

treg_LvsNL.markers <- FindAllMarkers(treg_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_treg_LvsNL_markers <- treg_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_treg_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(liu_treg_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_treg_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only monocytes
submono <- subset(liu, subset = h_celltype == "Mono")
unique(submono$h_celltype) #to check

# COMMAND ----------

mono_LvsNL <- subset(submono, Condition %in% c("lesional", "non lesional"))
unique(mono_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers monocytes

# COMMAND ----------

Idents(mono_LvsNL) <- mono_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

mono_LvsNL.markers <- FindAllMarkers(mono_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_mono_LvsNL_markers <- mono_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_mono_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(liu_mono_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_mono_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

# First it is performed a subset with only macrophages
submacro <- subset(liu, subset = h_celltype == "Macro")
unique(submacro$h_celltype) #to check

# COMMAND ----------

macro_LvsNL <- subset(submacro, Condition %in% c("lesional", "non lesional"))
unique(macro_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers macrophages

# COMMAND ----------

Idents(macro_LvsNL) <- macro_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

macro_LvsNL.markers <- FindAllMarkers(macro_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_macro_LvsNL_markers <- macro_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_macro_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(liu_macro_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_macro_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only DC
subdc <- subset(liu, subset = h_celltype == "DC")
unique(subdc$h_celltype) #to check

# COMMAND ----------

dc_LvsNL <- subset(subdc, Condition %in% c("lesional", "non lesional"))
unique(dc_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers DC

# COMMAND ----------

Idents(dc_LvsNL) <- dc_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

dc_LvsNL.markers <- FindAllMarkers(dc_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_dc_LvsNL_markers <- dc_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_dc_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(liu_dc_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_dc_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Natural killers

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only NK
subnk <- subset(liu, subset = h_celltype == "NK")
unique(subnk$h_celltype) #to check

# COMMAND ----------

nk_LvsNL <- subset(subnk, Condition %in% c("lesional", "non lesional"))
unique(nk_LvsNL$Condition) #to check

# COMMAND ----------

table(nk_LvsNL@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough cells but I did the analysis

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers NK

# COMMAND ----------

Idents(nk_LvsNL) <- nk_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

nk_LvsNL.markers <- FindAllMarkers(nk_LvsNL, only.pos=TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC Not significant results

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

# MAGIC %md
# MAGIC ###LvsNL

# COMMAND ----------

#First it is performed a subset with only mast cells
submast <- subset(liu, subset = h_celltype == "MastC")
unique(submast$h_celltype) #to check

# COMMAND ----------

mast_LvsNL <- subset(submast, Condition %in% c("lesional", "non lesional"))
unique(mast_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC
# MAGIC ### Find markers mastC

# COMMAND ----------

Idents(mast_LvsNL) <- mast_LvsNL$Condition

# COMMAND ----------

mast_LvsNL.markers <- FindAllMarkers(mast_LvsNL, only.pos=TRUE)

# COMMAND ----------

liu_mast_LvsNL_markers <- mast_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(liu_mast_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(liu_mast_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_mast_LvsNL_markers.xlsx")
