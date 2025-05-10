# Databricks notebook source
# MAGIC %md
# MAGIC #Find markers
# MAGIC ## AD/PS Reynolds et al, 2021
# MAGIC - Includes lesional and non-lesional samples
# MAGIC
# MAGIC Here will be performed the analysis to find all markers of each relevant cell type in Reynolds et al, 2021 dataset.
# MAGIC
# MAGIC ###Most relevant cell types: 
# MAGIC T-cells (TC), Treg, Innate linfoid cells (ILC), Fibroblasts, Keratinocytes (KC), Monocytes (Mono), Macrophagues (Macro), Dendritic cells (DC), Natural killers (NK) and MastC
# MAGIC
# MAGIC ###Constrasts: 
# MAGIC   - Lesional vs Healthy control (LvsHC)
# MAGIC   - Non lesional vs Healthy control (NLvsHC)
# MAGIC   - Lesional vs Non-lesional (LvsNL) <-- Here this contrast is done

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

#Read seurat object
reynolds <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

head(reynolds@meta.data, 5)

# COMMAND ----------

# MAGIC %md
# MAGIC I create a new column called Condition to keep the information of healthy or AD and lesional or non lesional all in one

# COMMAND ----------

# Create a new column 'Condition' with default value 'healthy'
reynolds$Condition <- "healthy"

# Update 'Condition' based on the 'Status' and 'Site' columns
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "lesion"] <- "lesional"
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "non_lesion"] <- "non lesional"

unique(reynolds$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC First check if all conditions in each cell type have at least 50 cells to perform the DEA analysis with enough power

# COMMAND ----------

# Check if all conditions have at least 50 cells for all cell types
# Print result for each cell type
cell_types <- unique(reynolds$h_celltype)
for(cell_type in cell_types) {
  sub <- subset(reynolds, subset = h_celltype == cell_type)
  if(all((table(sub$Condition)) >= 50)) {
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
subtcell <- subset(reynolds, subset = h_celltype == "TC")
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

reynolds_LvsNL_tcell_markers <- tcell_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_LvsNL_tcell_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_LvsNL_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_LvsNL_tcell_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##T reg

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only Treg
subtreg <- subset(reynolds, subset = h_celltype == "Treg")
unique(subtreg$h_celltype) #to check

# COMMAND ----------

treg_LvsNL <- subset(subtreg, Condition %in% c("lesional", "non lesional"))
unique(treg_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers Treg

# COMMAND ----------

Idents(treg_LvsNL) <- treg_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

treg_LvsNL.markers <- FindAllMarkers(treg_LvsNL, only.pos=TRUE)

# COMMAND ----------

treg_LvsNL_markers <- treg_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(treg_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(treg_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/treg_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only ILC
subILC <- subset(reynolds, subset = h_celltype == "ILC")
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

reynolds_ILC_LvsNL_markers <- ILC_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_ILC_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_ILC_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_ILC_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only ILC
subfb <- subset(reynolds, subset = h_celltype == "Fibroblasts")
unique(subfb$h_celltype) #to check

# COMMAND ----------

fb_LvsNL <- subset(subfb, Condition %in% c("lesional", "non lesional"))
unique(fb_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers fibroblasts

# COMMAND ----------

Idents(fb_LvsNL) <- fb_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

fb_LvsNL.markers <- FindAllMarkers(fb_LvsNL, only.pos=TRUE)

# COMMAND ----------

reynolds_fb_LvsNL_markers <- fb_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_fb_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_fb_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_fb_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only KC
subkc <- subset(reynolds, subset = h_celltype == "KC")
unique(subkc$h_celltype) #to check

# COMMAND ----------

kc_LvsNL <- subset(subkc, Condition %in% c("lesional", "non lesional"))
unique(kc_LvsNL$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers KC

# COMMAND ----------

Idents(kc_LvsNL) <- kc_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

kc_LvsNL.markers <- FindAllMarkers(kc_LvsNL, only.pos=TRUE)

# COMMAND ----------

reynolds_kc_LvsNL_markers <- kc_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_kc_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_kc_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_kc_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only monocytes
submono <- subset(reynolds, subset = h_celltype == "Mono")
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

reynolds_mono_LvsNL_markers <- mono_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_mono_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_mono_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_mono_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only macrophages
submacro <- subset(reynolds, subset = h_celltype == "Macro")
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

reynolds_macro_LvsNL_markers <- macro_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_macro_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_macro_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_macro_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only DC
subdc <- subset(reynolds, subset = h_celltype == "DC")
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

reynolds_dc_LvsNL_markers <- dc_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_dc_LvsNL_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_dc_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_dc_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Natural killers

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs NL
# MAGIC

# COMMAND ----------

#First it is performed a subset with only NK
subnk <- subset(reynolds, subset = h_celltype == "NK")
unique(subnk$h_celltype) #to check

# COMMAND ----------

nk_LvsNL <- subset(subnk, Condition %in% c("lesional", "non lesional"))
unique(nk_LvsNL$Condition) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(nk_LvsNL@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(nk_LvsNL@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough cells but I did it just to have

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers NK

# COMMAND ----------

Idents(nk_LvsNL) <- nk_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

nk_LvsNL.markers <- FindAllMarkers(nk_LvsNL, only.pos=TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs NL

# COMMAND ----------

#First it is performed a subset with only MAST
submast <- subset(reynolds, subset = h_celltype == "MastC")
unique(submast$h_celltype) #to check

# COMMAND ----------

mast_LvsNL <- subset(submast, Condition %in% c("lesional", "non lesional"))
unique(mast_LvsNL$Condition) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(mast_LvsNL$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(mast_LvsNL$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Just 6 cells below 50, I will perform the analysis but taking it into account

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers mastC

# COMMAND ----------

Idents(mast_LvsNL) <- mast_LvsNL$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

mast_LvsNL.markers <- FindAllMarkers(mast_LvsNL, only.pos=TRUE)

# COMMAND ----------

reynolds_mast_LvsNL_markers <- mast_LvsNL.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_mast_LvsNL_markers)

# COMMAND ----------

# Save results
write.xlsx(reynolds_mast_LvsNL_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_mast_LvsNL_markers.xlsx")
