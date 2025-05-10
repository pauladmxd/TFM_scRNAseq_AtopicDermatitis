# Databricks notebook source
# MAGIC %md
# MAGIC # Common markers from the different datasets
# MAGIC
# MAGIC Here I will merge the markers I found for each cell type on each dataset, for the **LvsHC** contrast

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Lesional vs Healthy Control

# COMMAND ----------

# MAGIC %md
# MAGIC ### T-cells

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_LvsHC_tcell_allmarkers.xlsx")

alkon_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_LvsHC_tcell_allmarkers.xlsx")

liu_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_LvsHC_tcell_allmarkers.xlsx")

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
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.liu == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3) %>%
  select(-cluster.liu, -cluster.reynolds, -cluster.alkon)

display(filtered_tcell_markersRLA)

# COMMAND ----------

#Filter all pairwise common markers
# For reynolds - liu
filtered_tcell_markersRL <- tcell_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "healthy" & cluster.reynolds == "healthy")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)

# For liu - alkon
filtered_tcell_markersLA <- tcell_markersLA %>%
  filter((cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.liu == "healthy" & cluster.alkon == "healthy")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.alkon)

# For reynolds- alkon
filtered_tcell_markersRA <- tcell_markersRA %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)


# COMMAND ----------

display(filtered_tcell_markersRL)

# COMMAND ----------

# MAGIC %md
# MAGIC As I see, alkon dataset is the most restringent

# COMMAND ----------

write.xlsx(filtered_tcell_markersRLA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/ARL_Tcell_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_tcell_markersRL, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_Tcell_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_tcell_markersRA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/AR_Tcell_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_tcell_markersLA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/AL_Tcell_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### Fibroblast

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_fb_LvsHC_allmarkers.xlsx")
alkon_LvsHC_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

fb_markers <- merge(alkon_LvsHC_fb_markers, reynolds_LvsHC_fb_markers, by ="gene",suffixes = c(".alkon", ".reynolds"))
display(fb_markers %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster_x is "AD" and cluster_y is "lesional" at the same time or healthy for both.
filtered_fb_markers <- fb_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)

display(filtered_fb_markers)

# COMMAND ----------

write.xlsx(filtered_fb_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/AR_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### Keratinocytes

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_kc_LvsHC_allmarkers.xlsx")
alkon_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

kc_markers <- merge(reynolds_LvsHC_KC_markers, alkon_LvsHC_KC_markers, by ="gene", suffixes= c(".reynolds", ".alkon"))
display(kc_markers %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster is the same for both datasets.
filtered_kc_markers <- kc_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)
display(filtered_kc_markers)

# COMMAND ----------

write.xlsx(filtered_kc_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/AR_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### DC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_dc_LvsHC_allmarkers.xlsx")

liu_LvsHC_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

dc_markers <- merge(reynolds_LvsHC_dc_markers, liu_LvsHC_dc_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(dc_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is the same for all the datasets.
filtered_dc_markers <- dc_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "healthy" & cluster.liu == "healthy")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_dc_markers)

# COMMAND ----------

write.xlsx(filtered_dc_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/RL_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###ILC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_ILC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_ILC_LvsHC_allmarkers.xlsx")
liu_ILC_LvsHC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_ILC_LvsHC_allmarkers.xlsx")

# COMMAND ----------

ILC_markers <- merge(liu_ILC_LvsHC_markers, reynolds_LvsHC_ILC_markers, by ="gene", suffixes= c(".liu", ".reynolds"))
display(ILC_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is the same for both
filtered_ILC_markers <- ILC_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "healthy" & cluster.liu == "healthy")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_ILC_markers)

# COMMAND ----------

write.xlsx(filtered_ILC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/RL_ilc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Macrophages

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_allmarkers.xlsx")

alkon_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_allmarkers.xlsx")

liu_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

alkon_LvsHC_macro_markers_suffix <- alkon_LvsHC_macro_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  

macro_markersRL <- merge(reynolds_LvsHC_macro_markers, liu_LvsHC_macro_markers, by = "gene", suffixes=c(".reynolds", ".liu"))
macro_markersAR <- merge(reynolds_LvsHC_macro_markers, alkon_LvsHC_macro_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
macro_markersAL <- merge(liu_LvsHC_macro_markers, alkon_LvsHC_macro_markers, by = "gene", suffixes=c(".alkon", ".liu"))
macro_markersARL <- merge(macro_markersRL, alkon_LvsHC_macro_markers_suffix, by = "gene") 
display(macro_markersARL %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_macro_markersARL <- macro_markersARL %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.liu == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3) %>%
  select(-cluster.liu, -cluster.reynolds, -cluster.alkon)

display(filtered_macro_markersARL)

# COMMAND ----------

# Filtering for all pairwise common markers

#For alkon-reynolds
filtered_macro_markersAR <- macro_markersAR %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)

#For alkon-liu
filtered_macro_markersAL <- macro_markersAL %>%
  filter((cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.liu == "healthy" & cluster.alkon == "healthy")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.alkon)

#For reynolds-liu
filtered_macro_markersRL <- macro_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "healthy" & cluster.reynolds == "healthy")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)


# COMMAND ----------

display(filtered_macro_markersRL)

# COMMAND ----------

write.xlsx(filtered_macro_markersARL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/ARL_macro_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_macro_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_macro_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_macro_markersAL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/AL_macro_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_macro_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/AR_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Monocytes

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mono_LvsHC_allmarkers.xlsx")
liu_LvsHC_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

mono_markers <- merge(reynolds_LvsHC_mono_markers, liu_LvsHC_mono_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mono_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_mono_markers <- mono_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "healthy" & cluster.liu == "healthy")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_mono_markers)

# COMMAND ----------

write.xlsx(filtered_mono_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/RL_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Treg

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_allmarkers.xlsx")

alkon_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_treg_LvsHC_allmarkers.xlsx")

liu_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

alkon_LvsHC_treg_markers_suffix <- alkon_LvsHC_treg_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  

treg_markersRL <- merge(reynolds_LvsHC_treg_markers, liu_LvsHC_treg_markers, by = "gene", suffixes=c(".reynolds", ".liu"))
treg_markersAR <- merge(reynolds_LvsHC_treg_markers, alkon_LvsHC_treg_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
treg_markersAL <- merge(liu_LvsHC_treg_markers, alkon_LvsHC_treg_markers, by = "gene", suffixes=c(".alkon", ".liu"))
treg_markersARL <- merge(treg_markersRL, alkon_LvsHC_treg_markers_suffix, by = "gene") 
display(treg_markersARL %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_treg_markersARL <- treg_markersARL %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.liu == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu + avg_log2FC.alkon) / 3) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu + p_val.alkon) / 3) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu + p_val_adj.alkon) / 3) %>%
  select(-cluster.liu, -cluster.reynolds, -cluster.alkon)

display(filtered_treg_markersARL)

# COMMAND ----------

# Filtering for all pairwise common markers

#For alkon-reynolds
filtered_treg_markersAR <- treg_markersAR %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)

#For alkon-liu
filtered_treg_markersAL <- treg_markersAL %>%
  filter((cluster.liu == "lesional" & cluster.alkon == "AD") | 
         (cluster.liu == "healthy" & cluster.alkon == "healthy")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.alkon)

#For reynolds-liu
filtered_treg_markersRL <- treg_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "healthy" & cluster.reynolds == "healthy")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)

display(filtered_treg_markersRL %>% arrange(avg_pvalue_adj))

# COMMAND ----------

# MAGIC %md
# MAGIC In this case Alkon markers are very restrictive, also taking into account that we found only 45 cells o Treg in HC condition.

# COMMAND ----------

write.xlsx(filtered_treg_markersARL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/ARL_treg_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_treg_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_treg_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_treg_markersAL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/AL_treg_LvsHC_allmarkers.xlsx")
write.xlsx(filtered_treg_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/AR_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###NK

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_nk_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_nk_LvsHC_allmarkers.xlsx")
liu_LvsHC_nk_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_nk_LvsHC_allmarkers.xlsx")

# COMMAND ----------

nk_markers <- merge(reynolds_LvsHC_nk_markers, liu_LvsHC_nk_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(nk_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_nk_markers <- nk_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "healthy" & cluster.liu == "healthy")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)

display(filtered_nk_markers)

# COMMAND ----------

write.xlsx(filtered_nk_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/NK/RL_nk_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###MastC

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mast_LvsHC_allmarkers.xlsx")
liu_LvsHC_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mast_LvsHC_allmarkers.xlsx")

# COMMAND ----------

mast_markers <- merge(reynolds_LvsHC_mast_markers, liu_LvsHC_mast_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mast_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_mast_markers <- mast_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "healthy" & cluster.liu == "healthy")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_mast_markers)

# COMMAND ----------

write.xlsx(filtered_mast_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/MastC/RL_mast_LvsHC_allmarkers.xlsx")
