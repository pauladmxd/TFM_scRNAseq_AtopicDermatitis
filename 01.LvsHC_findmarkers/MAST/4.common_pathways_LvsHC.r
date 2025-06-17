# Databricks notebook source
# MAGIC %md
# MAGIC #Common pathways across datasets 
# MAGIC Check each cell type if the pathways enriched with a significant pvalue are similar and with the same sign in the NES

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(openxlsx)
library(tidyverse)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

# MAGIC %md
# MAGIC ###Reactome

# COMMAND ----------

reynolds_kc_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/GSEA_individually/res_kc_0.05_reynolds.xlsx", sheet =2)
alkon_kc_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/GSEA_individually/res_kc_0.05_alkon.xlsx", sheet=2)

# COMMAND ----------

common_kc_reactome <- reynolds_kc_react[reynolds_kc_react$ID %in% alkon_kc_react$ID & sign(reynolds_kc_react$NES) == sign(alkon_kc_react$NES[match(reynolds_kc_react$ID, alkon_kc_react$ID)]), ]

# COMMAND ----------

common_kc_reactome_contrary <- reynolds_kc_react[reynolds_kc_react$ID %in% alkon_kc_react$ID & sign(reynolds_kc_react$NES) != sign(alkon_kc_react$NES[match(reynolds_kc_react$ID, alkon_kc_react$ID)]), ]
common_kc_reactome_contrary

# COMMAND ----------

rows <- min(nrow(reynolds_kc_react), nrow(alkon_kc_react))

# COMMAND ----------

display(reynolds_kc_react)

# COMMAND ----------

display(select(common_kc_reactome, -ncol(common_kc_reactome)))

# COMMAND ----------

percent_similarity_kc <- nrow(common_kc_reactome) / rows *100
percent_similarity_kc

# COMMAND ----------

# MAGIC %md
# MAGIC 34% similarity

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tcell

# COMMAND ----------

# MAGIC %md
# MAGIC ###reactome

# COMMAND ----------

reynolds_tcell_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/GSEA_individually/res_tcell_0.05_reynolds.xlsx", sheet =2)
alkon_tcell_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/MAST_method/GSEA_individually/res_tcell_0.05_alkon.xlsx", sheet=2)

# COMMAND ----------

rows <- min(nrow(alkon_tcell_react), nrow(reynolds_tcell_react))

# COMMAND ----------

common_tcell_react <- reynolds_tcell_react[reynolds_tcell_react$ID %in% alkon_tcell_react$ID & sign(reynolds_tcell_react$NES) == sign(alkon_tcell_react$NES[match(reynolds_tcell_react$ID, alkon_tcell_react$ID)]), ]
common_tcell_react

# COMMAND ----------

common_tcell_react_contrary <- reynolds_tcell_react[reynolds_tcell_react$ID %in% alkon_tcell_react$ID & sign(reynolds_tcell_react$NES) != sign(alkon_tcell_react$NES[match(reynolds_tcell_react$ID, alkon_tcell_react$ID)]), ]
common_tcell_react_contrary

# COMMAND ----------

percent_similarity_tcell_react <- nrow(common_tcell_react)/rows * 100
percent_similarity_tcell_react

# COMMAND ----------

# MAGIC %md
# MAGIC 0% similarity

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

# MAGIC %md
# MAGIC ###reactome

# COMMAND ----------

reynolds_fb_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/MAST_method/GSEA_individually/res_fb_0.05_reynolds.xlsx", sheet =2)
alkon_fb_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_fb_0.05_alkon.xlsx", sheet=2)

# COMMAND ----------

common_fb_reactome <- reynolds_fb_react[reynolds_fb_react$ID %in% alkon_fb_react$ID & sign(reynolds_fb_react$NES) == sign(alkon_fb_react$NES[match(reynolds_fb_react$ID, alkon_fb_react$ID)]), ]
common_fb_reactome

# COMMAND ----------

common_fb_reactome_contrary <- reynolds_fb_react[reynolds_fb_react$ID %in% alkon_fb_react$ID & sign(reynolds_fb_react$NES) != sign(alkon_fb_react$NES[match(reynolds_fb_react$ID, alkon_fb_react$ID)]), ]
common_fb_reactome_contrary

# COMMAND ----------

rows <- min(nrow(reynolds_fb_react), nrow(alkon_fb_react))

# COMMAND ----------

percent_similarity_fb <- nrow(common_fb_reactome) / rows *100
percent_similarity_fb

# COMMAND ----------

# MAGIC %md
# MAGIC 30% of the terms are common with same sign

# COMMAND ----------

percent_similarity_contrary_fb <- nrow(common_fb_reactome_contrary) / rows *100
percent_similarity_contrary_fb

# COMMAND ----------

# MAGIC %md
# MAGIC 50% are common terms with contrary NES sign
