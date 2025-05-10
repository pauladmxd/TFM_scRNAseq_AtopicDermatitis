# Databricks notebook source
# MAGIC %md
# MAGIC # Common DEGs from the different datasets
# MAGIC
# MAGIC Here I will merge the DEGs I found for each cell type on each dataset using the **pseudobulk** approach, for the **LvsHC** contrast

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)
library(dplyr)
library(openxlsx)
.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(VennDiagram)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Lesional vs Healthy Control

# COMMAND ----------

# MAGIC %md
# MAGIC ### T-cells

# COMMAND ----------

#Read each dataset DEGs saved previously
reynolds_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_tcell_LvsHC_bulk.xlsx")

alkon_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_tcell_LvsHC_bulk.xlsx")

liu_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_tcell_LvsHC_bulk.xlsx")

# COMMAND ----------

# Merge the first two data frames with suffixes
tcell_markersRL <- merge(reynolds_LvsHC_tcell_markers, liu_LvsHC_tcell_markers, by = "gene", suffixes = c(".reynolds", ".liu"))

alkon_LvsHC_tcell_markers_suffix <- alkon_LvsHC_tcell_markers %>%
rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  
# Merge the result with the third data frame
tcell_markersRLA <- merge(tcell_markersRL, alkon_LvsHC_tcell_markers_suffix, by = "gene")

# Display the merged data frame sorted by p_val_adj
display(tcell_markersRLA %>% arrange(p_val_adj.alkon))

# COMMAND ----------

colnames(tcell_markersRLA)

# COMMAND ----------

#To save all the pairwise common markers with custom column names
tcell_markersRA <- merge(reynolds_LvsHC_tcell_markers, alkon_LvsHC_tcell_markers, by = "gene", suffixes = c(".reynolds", ".alkon"))
tcell_markersLA <- merge(liu_LvsHC_tcell_markers, alkon_LvsHC_tcell_markers, by = "gene", suffixes = c(".liu", ".alkon"))

# COMMAND ----------

display(tcell_markersRA)

# COMMAND ----------

filtered_tcell_markersRLA <- tcell_markersRLA %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu) & sign(avg_log2FC.liu) == sign(avg_log2FC.alkon)))  %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3)

display(filtered_tcell_markersRLA)

# COMMAND ----------

#Filter all pairwise common markers
# For reynolds - liu
filtered_tcell_markersRL <- tcell_markersRL %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu))) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) 

# For liu - alkon
filtered_tcell_markersLA <- tcell_markersLA %>%
  filter((sign(avg_log2FC.liu) == sign(avg_log2FC.alkon))) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2)
# For reynolds- alkon
filtered_tcell_markersRA <- tcell_markersRA %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.alkon))) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) 


# COMMAND ----------

display(filtered_tcell_markersRL)

# COMMAND ----------

write.xlsx(filtered_tcell_markersRLA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Tcell/ARL_Tcell_LvsHC_bulk.xlsx")
write.xlsx(filtered_tcell_markersRL, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Tcell/RL_Tcell_LvsHC_bulk.xlsx")
write.xlsx(filtered_tcell_markersRA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Tcell/AR_Tcell_LvsHC_bulk.xlsx")
write.xlsx(filtered_tcell_markersLA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Tcell/AL_Tcell_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### Fibroblast

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_fb_LvsHC_bulk.xlsx")
alkon_LvsHC_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_fb_LvsHC_bulk.xlsx")

# COMMAND ----------

fb_markers <- merge(alkon_LvsHC_fb_markers, reynolds_LvsHC_fb_markers, by ="gene",suffixes = c(".alkon", ".reynolds"))
display(fb_markers %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where the sign of avg_log2FC is the same for both datasets.
filtered_fb_markers <- fb_markers %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.alkon)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2)
display(filtered_fb_markers)

# COMMAND ----------

write.xlsx(filtered_fb_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Fibroblast/AR_fb_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### Keratinocytes

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_kc_LvsHC_bulk.xlsx")
alkon_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_kc_LvsHC_bulk.xlsx")

# COMMAND ----------

kc_markers <- merge(reynolds_LvsHC_KC_markers, alkon_LvsHC_KC_markers, by ="gene", suffixes= c(".reynolds", ".alkon"))
display(kc_markers %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where avg_log2FC has the same sign.
filtered_kc_markers <- kc_markers %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.alkon)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2)
display(filtered_kc_markers)

# COMMAND ----------

write.xlsx(filtered_kc_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Keratinocytes/AR_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ####Common DEGs between the different methods and the different datasets

# COMMAND ----------

filtered_kc_markers_sc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/AR_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

filtered_kc_markers_bulk <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Keratinocytes/AR_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

common_kc_markers <- merge(filtered_kc_markers_sc, filtered_kc_markers_bulk, by ="gene", suffixes= c(".sc", ".bulk"))
display(common_kc_markers)

# COMMAND ----------

# Filter rows where avg_log2FC has the same sign.
filtered_common_kc_markers <- common_kc_markers %>%
  filter(sign(avg_log2FC.sc) == sign(avg_log2FC.bulk)) %>%
  mutate(avg_log2FC = (avg_log2FC.sc + avg_log2FC.bulk) / 2) %>%
  mutate(avg_pvalue = (avg_pvalue.sc + avg_pvalue.bulk) / 2) %>%
  mutate(avg_pvalue_adj = (avg_pvalue_adj.sc + avg_pvalue_adj.bulk) / 2)
display(filtered_common_kc_markers)

# COMMAND ----------

filtered_common_kc_markers <- filtered_common_kc_markers %>%
  filter((avg_pvalue_adj.sc < 0.05) & (avg_pvalue_adj.bulk < 0.05) & (avg_log2FC.sc > 0.5) & (avg_log2FC.bulk > 0.5))
display(filtered_common_kc_markers)

# COMMAND ----------

display(filtered_kc_markers_sc)

# COMMAND ----------

kc.de.sig <- filtered_kc_markers_sc[filtered_kc_markers_sc$avg_pvalue_adj < 0.05 & filtered_kc_markers_sc$avg_log2FC> 0.5,]
bulk.kc.de.sig <- filtered_kc_markers_bulk[filtered_kc_markers_bulk$avg_pvalue_adj < 0.05 & filtered_kc_markers_bulk$avg_log2FC > 0.5,]

common <- merge(kc.de.sig, bulk.kc.de.sig, by= "gene")

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = na.omit(kc.de.sig$gene),
  "Bulk" = na.omit(bulk.kc.de.sig$gene)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs AR Keratinocytes", # Add title,
  cat.dist = c(0.04, 0.04), # Adjust the distance of the category names from the circles
  main.cex = 2, # Increase title size
  cat.cex = 1.5, # Increase label size
  cat.pos = c(-20, 20), # Position labels more on the top
  cex = 1.5 # Increase numbers size
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ### DC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_dc_LvsHC_bulk.xlsx")

liu_LvsHC_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_dc_LvsHC_bulk.xlsx")

# COMMAND ----------

dc_markers <- merge(reynolds_LvsHC_dc_markers, liu_LvsHC_dc_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(dc_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where avg_log2FC has the same sign.
filtered_dc_markers <- dc_markers %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) 
display(filtered_dc_markers)

# COMMAND ----------

write.xlsx(filtered_dc_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/DC/RL_dc_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###ILC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_ILC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_ilc_LvsHC_bulk.xlsx")
liu_ILC_LvsHC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_ilc_LvsHC_bulk.xlsx")

# COMMAND ----------

ILC_markers <- merge(liu_ILC_LvsHC_markers, reynolds_LvsHC_ILC_markers, by ="gene", suffixes= c(".liu", ".reynolds"))
display(ILC_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where sign in log2FC is the same for both
filtered_ILC_markers <- ILC_markers %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2)
display(filtered_ILC_markers)

# COMMAND ----------

write.xlsx(filtered_ILC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/ILC/RL_ilc_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Macrophages

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_macro_LvsHC_bulk.xlsx")

alkon_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_macro_LvsHC_bulk.xlsx")

liu_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_macro_LvsHC_bulk.xlsx")

# COMMAND ----------

alkon_LvsHC_macro_markers_suffix <- alkon_LvsHC_macro_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  

macro_markersRL <- merge(reynolds_LvsHC_macro_markers, liu_LvsHC_macro_markers, by = "gene", suffixes=c(".reynolds", ".liu"))
macro_markersAR <- merge(reynolds_LvsHC_macro_markers, alkon_LvsHC_macro_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
macro_markersAL <- merge(liu_LvsHC_macro_markers, alkon_LvsHC_macro_markers, by = "gene", suffixes=c(".alkon", ".liu"))
macro_markersARL <- merge(macro_markersRL, alkon_LvsHC_macro_markers_suffix, by = "gene") 
display(macro_markersARL %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where the sign of avglog2FC is the same
filtered_macro_markersARL <- macro_markersARL %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu) & sign(avg_log2FC.liu) == sign(avg_log2FC.alkon)))  %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3) 
filtered_macro_markersARL

# COMMAND ----------

# Filtering for all pairwise common markers

# For alkon-reynolds
filtered_macro_markersAR <- macro_markersAR %>%
  filter(sign(avg_log2FC.alkon) == sign(avg_log2FC.reynolds)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2)

# For alkon-liu
filtered_macro_markersAL <- macro_markersAL %>%
  filter(sign(avg_log2FC.alkon) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2) 
# For reynolds-liu
filtered_macro_markersRL <- macro_markersRL %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2)

# COMMAND ----------

display(filtered_macro_markersRL)

# COMMAND ----------

write.xlsx(filtered_macro_markersARL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Macrophage/ARL_macro_LvsHC_bulk.xlsx")
write.xlsx(filtered_macro_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Macrophage/RL_macro_LvsHC_bulk.xlsx")
write.xlsx(filtered_macro_markersAL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Macrophage/AL_macro_LvsHC_bulk.xlsx")
write.xlsx(filtered_macro_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Macrophage/AR_macro_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Monocytes

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_mono_LvsHC_bulk.xlsx")
liu_LvsHC_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_mono_LvsHC_bulk.xlsx")

# COMMAND ----------

mono_markers <- merge(reynolds_LvsHC_mono_markers, liu_LvsHC_mono_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mono_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where sign of avglog2FC is the same
filtered_mono_markers <- mono_markers %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2)
display(filtered_mono_markers)

# COMMAND ----------

write.xlsx(filtered_mono_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Monocyte/RL_mono_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Treg

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_treg_LvsHC_bulk.xlsx")

alkon_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_treg_LvsHC_bulk.xlsx")

liu_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_treg_LvsHC_bulk.xlsx")

# COMMAND ----------

alkon_LvsHC_treg_markers_suffix <- alkon_LvsHC_treg_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  

treg_markersRL <- merge(reynolds_LvsHC_treg_markers, liu_LvsHC_treg_markers, by = "gene", suffixes=c(".reynolds", ".liu"))
treg_markersAR <- merge(reynolds_LvsHC_treg_markers, alkon_LvsHC_treg_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
treg_markersAL <- merge(liu_LvsHC_treg_markers, alkon_LvsHC_treg_markers, by = "gene", suffixes=c(".alkon", ".liu"))
treg_markersARL <- merge(treg_markersRL, alkon_LvsHC_treg_markers_suffix, by = "gene") 
display(treg_markersARL %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where sign of avglog2fc is the same for all of them
filtered_treg_markersARL <- treg_markersARL %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu) & sign(avg_log2FC.liu) == sign(avg_log2FC.alkon)))  %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3)
display(filtered_treg_markersARL)

# COMMAND ----------

# Filtering for all pairwise common markers

# For alkon-reynolds
filtered_treg_markersAR <- treg_markersAR %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.alkon)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2)

# For alkon-liu
filtered_treg_markersAL <- treg_markersAL %>%
  filter(sign(avg_log2FC.liu) == sign(avg_log2FC.alkon)) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2) 

# For reynolds-liu
filtered_treg_markersRL <- treg_markersRL %>%
  filter(sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu)) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) 

display(filtered_treg_markersRL %>% arrange(avg_pvalue_adj))

# COMMAND ----------

write.xlsx(filtered_treg_markersARL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Treg/ARL_treg_LvsHC_bulk.xlsx")
write.xlsx(filtered_treg_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Treg/RL_treg_LvsHC_bulk.xlsx")
write.xlsx(filtered_treg_markersAL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Treg/AL_treg_LvsHC_bulk.xlsx")
write.xlsx(filtered_treg_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/Treg/AR_treg_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###NK

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_nk_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_nk_LvsHC_bulk.xlsx")
liu_LvsHC_nk_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_nk_LvsHC_bulk.xlsx")

# COMMAND ----------

nk_markers <- merge(reynolds_LvsHC_nk_markers, liu_LvsHC_nk_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(nk_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_nk_markers <- nk_markers %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu))) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2)
display(filtered_nk_markers)

# COMMAND ----------

write.xlsx(filtered_nk_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/NK/RL_nk_LvsHC_bulk.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###MastC

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_mast_LvsHC_bulk.xlsx")
liu_LvsHC_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_mast_LvsHC_bulk.xlsx")

# COMMAND ----------

mast_markers <- merge(reynolds_LvsHC_mast_markers, liu_LvsHC_mast_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mast_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_mast_markers <- mast_markers %>%
  filter((sign(avg_log2FC.reynolds) == sign(avg_log2FC.liu))) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) 
display(filtered_mast_markers)

# COMMAND ----------

write.xlsx(filtered_mast_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/pseudobulk/MastC/RL_mast_LvsHC_allmarkers.xlsx")
