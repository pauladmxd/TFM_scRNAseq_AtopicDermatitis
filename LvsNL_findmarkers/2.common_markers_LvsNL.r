# Databricks notebook source
# MAGIC %md
# MAGIC # Common markers from the different datasets
# MAGIC
# MAGIC Here I will merge the markers I found for each cell type on each dataset, for the **LvsNL** contrast

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Lesional vs Non lesional

# COMMAND ----------

# MAGIC %md
# MAGIC ### T-cells

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsNL_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_LvsNL_tcell_markers.xlsx")

liu_LvsNL_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_LvsNL_tcell_markers.xlsx")

# COMMAND ----------

# Merge the two data frames with suffixes
tcell_markersRL <- merge(reynolds_LvsNL_tcell_markers, liu_LvsNL_tcell_markers, by = "gene", suffixes = c(".reynolds", ".liu"))

# Display the merged data frame sorted by p_val_adj
display(tcell_markersRL %>% arrange(p_val_adj.liu))

# COMMAND ----------

#Filter common markers and add average pvalues and FC
filtered_tcell_markersRL <- tcell_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "non lesional" & cluster.reynolds == "non lesional")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)


# COMMAND ----------

display(filtered_tcell_markersRL)

# COMMAND ----------

write.xlsx(filtered_tcell_markersRL, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_Tcell_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ### DC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsNL_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_dc_LvsNL_markers.xlsx")

liu_LvsNL_dc_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_dc_LvsNL_markers.xlsx")

# COMMAND ----------

dc_markers <- merge(reynolds_LvsNL_dc_markers, liu_LvsNL_dc_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(dc_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is the same for all the datasets.
filtered_dc_markers <- dc_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "non lesional" & cluster.liu == "non lesional")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_dc_markers)

# COMMAND ----------

write.xlsx(filtered_dc_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/RL_dc_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###ILC

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsNL_ILC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_ILC_LvsNL_markers.xlsx")
liu_ILC_LvsNL_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_ILC_LvsNL_markers.xlsx")

# COMMAND ----------

ILC_markers <- merge(liu_ILC_LvsNL_markers, reynolds_LvsNL_ILC_markers, by ="gene", suffixes= c(".liu", ".reynolds"))
display(ILC_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is the same for both
filtered_ILC_markers <- ILC_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "non lesional" & cluster.liu == "non lesional")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_ILC_markers)

# COMMAND ----------

write.xlsx(filtered_ILC_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/RL_ilc_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Macrophages

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsNL_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_macro_LvsNL_markers.xlsx")

liu_LvsNL_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_macro_LvsNL_markers.xlsx")

# COMMAND ----------


macro_markersRL <- merge(reynolds_LvsNL_macro_markers, liu_LvsNL_macro_markers, by = "gene", suffixes=c(".reynolds", ".liu"))
display(macro_markersRL)

# COMMAND ----------

# Filtering common markers and add avg FC and pvlue
filtered_macro_markersRL <- macro_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "non lesional" & cluster.reynolds == "non lesional")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)


# COMMAND ----------

display(filtered_macro_markersRL)

# COMMAND ----------

write.xlsx(filtered_macro_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_macro_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Monocytes

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsNL_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_mono_LvsNL_markers.xlsx")
liu_LvsNL_mono_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_mono_LvsNL_markers.xlsx")

# COMMAND ----------

mono_markers <- merge(reynolds_LvsNL_mono_markers, liu_LvsNL_mono_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mono_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "AD" or "lesional" for all of them or healthy for all
filtered_mono_markers <- mono_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "non lesional" & cluster.liu == "non lesional")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_mono_markers)

# COMMAND ----------

write.xlsx(filtered_mono_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/RL_mono_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Treg

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsNL_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/treg_LvsNL_markers.xlsx")
liu_LvsNL_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_treg_LvsNL_markers.xlsx")

# COMMAND ----------

treg_markersRL <- merge(reynolds_LvsNL_treg_markers, liu_LvsNL_treg_markers, by = "gene", suffixes=c(".reynolds", ".liu"))

# COMMAND ----------

#Filtering common markers and add avregae FC and pvalues
filtered_treg_markersRL <- treg_markersRL %>%
  filter((cluster.liu == "lesional" & cluster.reynolds == "lesional") | 
         (cluster.liu == "non lesional" & cluster.reynolds == "non lesional")) %>%
  mutate(Condition = cluster.reynolds) %>%
  mutate(avg_log2FC = (avg_log2FC.reynolds + avg_log2FC.liu) / 2) %>%
  mutate(avg_pvalue = (p_val.reynolds + p_val.liu) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.reynolds + p_val_adj.liu) / 2) %>%
  select(-cluster.liu, -cluster.reynolds)

display(filtered_treg_markersRL %>% arrange(avg_pvalue_adj))

# COMMAND ----------


write.xlsx(filtered_treg_markersRL, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_treg_LvsNL_markers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###MastC

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsNL_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsNL/reynolds_mast_LvsNL_markers.xlsx")
liu_LvsNL_mast_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsNL/liu_mast_LvsNL_markers.xlsx")

# COMMAND ----------

mast_markers <- merge(reynolds_LvsNL_mast_markers, liu_LvsNL_mast_markers, by = "gene", suffixes= c(".reynolds", ".liu"))
display(mast_markers %>% arrange(p_val_adj.liu))

# COMMAND ----------

# Filter rows where cluster is "lesional" for all of them or non lesional for all
filtered_mast_markers <- mast_markers %>%
  filter((cluster.reynolds == "lesional" & cluster.liu == "lesional") | (cluster.reynolds == "non lesional" & cluster.liu == "non lesional")) %>%
  mutate(Condition = cluster.liu) %>%
  mutate(avg_log2FC = (avg_log2FC.liu + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.liu + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.liu + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.liu)
display(filtered_mast_markers)

# COMMAND ----------

write.xlsx(filtered_mast_markers, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/MastC/RL_mast_LvsNL_markers.xlsx")
