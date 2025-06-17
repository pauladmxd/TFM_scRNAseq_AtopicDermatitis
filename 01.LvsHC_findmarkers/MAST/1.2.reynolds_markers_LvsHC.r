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

tcell_LvsHC.markers <- FindAllMarkers(tcell_LvsHC, only.pos=TRUE, test.use="MAST" )

# COMMAND ----------

reynolds_LvsHC_tcell_markers <- tcell_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_LvsHC_tcell_markers)

# COMMAND ----------

# MAGIC %sh
# MAGIC mkdir /dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method

# COMMAND ----------

#Save results
write.xlsx(tcell_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_LvsHC_tcell_allmarkers.xlsx")

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

treg_LvsHC.markers <- FindAllMarkers(treg_LvsHC, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

reynolds_treg_LvsHC_markers <- treg_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_treg_LvsHC_markers)

# COMMAND ----------

#Save results
# write.xlsx(treg_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_markers.xlsx")
write.xlsx(treg_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_treg_LvsHC_allmarkers.xlsx")

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

fb_LvsHC.markers <- FindAllMarkers(fb_LvsHC, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

reynolds_fb_LvsHC_markers <- fb_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_fb_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_fb_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_fb_LvsHC_markers.xlsx")
write.xlsx(fb_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_fb_LvsHC_allmarkers.xlsx")

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

kc_LvsHC.markers <- FindAllMarkers(kc_LvsHC, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

reynolds_kc_LvsHC_markers <- kc_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_kc_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_kc_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_kc_LvsHC_markers.xlsx")
write.xlsx(kc_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_kc_LvsHC_allmarkers.xlsx")

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

macro_LvsHC.markers <- FindAllMarkers(macro_LvsHC, only.pos=TRUE, test.use="MAST")

# COMMAND ----------

reynolds_macro_LvsHC_markers <- macro_LvsHC.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5 & p_val_adj < 0.05)
display(reynolds_macro_LvsHC_markers)

# COMMAND ----------

#Save results
#write.xlsx(reynolds_macro_LvsHC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_markers.xlsx")
write.xlsx(macro_LvsHC.markers, file="/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_macro_LvsHC_allmarkers.xlsx")

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

subkc.markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

  subkc.markers <- subkc.markers %>%
    mutate(avg_log2FC = ifelse(cluster == "healthy", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(subkc.markers)

# COMMAND ----------

subtcell.markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/reynolds_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

  subtcell.markers <- subtcell.markers %>%
    mutate(avg_log2FC = ifelse(cluster == "healthy", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(subtcell.markers)
