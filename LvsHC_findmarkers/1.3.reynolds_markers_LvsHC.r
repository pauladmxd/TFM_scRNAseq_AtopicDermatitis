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

reynolds_LvsHC_tcell_markers <- tcell_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_LvsHC_tcell_markers)

# COMMAND ----------

#Save results
write.xlsx(reynolds_LvsHC_tcell_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_LvsHC_tcell_markers.xlsx")
write.xlsx(tcell_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##T reg

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only Treg
subtreg <- subset(reynolds, subset = h_celltype == "Treg")
unique(subtreg$h_celltype) #to check

# COMMAND ----------

treg_LvsHC <- subset(subtreg, Condition %in% c("lesional", "healthy"))
unique(treg_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers Treg

# COMMAND ----------

Idents(treg_LvsHC) <- treg_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

treg_LvsHC.markers <- FindAllMarkers(treg_LvsHC, only.pos=TRUE)

# COMMAND ----------

reynolds_treg_LvsHC_markers <- treg_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_treg_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(treg_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_markers.xlsx")
write.xlsx(treg_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only ILC
subILC <- subset(reynolds, subset = h_celltype == "ILC")
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

reynolds_ILC_LvsHC_markers <- ILC_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_ILC_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_ILC_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_ILC_LvsHC_markers.xlsx")
write.xlsx(ILC_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_ILC_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only ILC
subfb <- subset(reynolds, subset = h_celltype == "Fibroblasts")
unique(subfb$h_celltype) #to check

# COMMAND ----------

fb_LvsHC <- subset(subfb, Condition %in% c("lesional", "healthy"))
unique(fb_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ###Find markers fibroblasts

# COMMAND ----------

Idents(fb_LvsHC) <- fb_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

fb_LvsHC.markers <- FindAllMarkers(fb_LvsHC, only.pos=TRUE)

# COMMAND ----------

reynolds_fb_LvsHC_markers <- fb_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_fb_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_fb_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_fb_LvsHC_markers.xlsx")
write.xlsx(fb_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only KC
subkc <- subset(reynolds, subset = h_celltype == "KC")
unique(subkc$h_celltype) #to check

# COMMAND ----------

kc_LvsHC <- subset(subkc, Condition %in% c("lesional", "healthy"))
unique(kc_LvsHC$Condition) #to check

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers KC

# COMMAND ----------

Idents(kc_LvsHC) <- kc_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

kc_LvsHC.markers <- FindAllMarkers(kc_LvsHC, only.pos=TRUE)

# COMMAND ----------

reynolds_kc_LvsHC_markers <- kc_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_kc_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_kc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_kc_LvsHC_markers.xlsx")
write.xlsx(kc_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only monocytes
submono <- subset(reynolds, subset = h_celltype == "Mono")
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

reynolds_mono_LvsHC_markers <- mono_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_mono_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_mono_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mono_LvsHC_markers.xlsx")
write.xlsx(mono_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only macrophages
submacro <- subset(reynolds, subset = h_celltype == "Macro")
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

reynolds_macro_LvsHC_markers <- macro_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_macro_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_macro_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_markers.xlsx")
write.xlsx(macro_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only DC
subdc <- subset(reynolds, subset = h_celltype == "DC")
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

reynolds_dc_LvsHC_markers <- dc_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_dc_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_dc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_dc_LvsHC_markers.xlsx")
write.xlsx(dc_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Natural killers

# COMMAND ----------

# MAGIC %md
# MAGIC ### L vs HC
# MAGIC

# COMMAND ----------

#First it is performed a subset with only NK
subnk <- subset(reynolds, subset = h_celltype == "NK")
unique(subnk$h_celltype) #to check

# COMMAND ----------

nk_LvsHC <- subset(subnk, Condition %in% c("lesional", "healthy"))
unique(nk_LvsHC$Condition) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(nk_LvsHC@meta.data$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(nk_LvsHC@meta.data$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough cells but I did it just to have

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers NK

# COMMAND ----------

Idents(nk_LvsHC) <- nk_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

nk_LvsHC.markers <- FindAllMarkers(nk_LvsHC, only.pos=TRUE)

# COMMAND ----------

display(nk_LvsHC.markers)

# COMMAND ----------

reynolds_nk_LvsHC_markers <- nk_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_nk_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_nk_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_nk_LvsHC_markers.xlsx")
write.xlsx(nk_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_nk_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

# MAGIC %md
# MAGIC ###L vs HC

# COMMAND ----------

#First it is performed a subset with only MAST
submast <- subset(reynolds, subset = h_celltype == "MastC")
unique(submast$h_celltype) #to check

# COMMAND ----------

mast_LvsHC <- subset(submast, Condition %in% c("lesional", "healthy"))
unique(mast_LvsHC$Condition) #to check

# COMMAND ----------

# Check if all conditions have at least 50 cells
# Print result
if(all((table(mast_LvsHC$Condition)) >= 50)) {
  print("All conditions have at least 50 cells.")
} else {
  print("Not all conditions have at least 50 cells.")
}
table(mast_LvsHC$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC Just 6 cells below 50, I will perform the analysis but taking it into account

# COMMAND ----------

# MAGIC %md
# MAGIC ### Find markers mastC

# COMMAND ----------

Idents(mast_LvsHC) <- mast_LvsHC$Condition # We want condition to be the contrast to find markers

# COMMAND ----------

mast_LvsHC.markers <- FindAllMarkers(mast_LvsHC, only.pos=TRUE)

# COMMAND ----------

reynolds_mast_LvsHC_markers <- mast_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_mast_LvsHC_markers)

# COMMAND ----------

# Save results
#write.xlsx(reynolds_mast_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mast_LvsHC_markers.xlsx")
write.xlsx(mast_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mast_LvsHC_allmarkers.xlsx")
