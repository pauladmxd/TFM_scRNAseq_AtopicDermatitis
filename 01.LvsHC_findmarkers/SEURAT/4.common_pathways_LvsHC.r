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
# MAGIC ###KEGG

# COMMAND ----------

alkon_kc_kegg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_kc_0.05_alkon.xlsx", sheet=1)

# COMMAND ----------

# MAGIC %md
# MAGIC no terms in alkon

# COMMAND ----------

# MAGIC %md
# MAGIC ###Reactome

# COMMAND ----------

reynolds_kc_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_kc_0.05_reynolds.xlsx", sheet =2)
alkon_kc_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_kc_0.05_alkon.xlsx", sheet=2)

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
# MAGIC 18% similarity

# COMMAND ----------

# MAGIC %md
# MAGIC ###Hallmarks

# COMMAND ----------

reynolds_kc_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_kc_0.05_reynolds.xlsx", sheet =3)
alkon_kc_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_kc_0.05_alkon.xlsx", sheet=3)

# COMMAND ----------

common_kc_hallmarks <- reynolds_kc_hall[reynolds_kc_hall$ID %in% alkon_kc_hall$ID & sign(reynolds_kc_hall$NES) == sign(alkon_kc_hall$NES[match(reynolds_kc_hall$ID, alkon_kc_hall$ID)]), ]

# COMMAND ----------

common_kc_hallmarks_contrary <- reynolds_kc_hall[reynolds_kc_hall$ID %in% alkon_kc_hall$ID & sign(reynolds_kc_hall$NES) != sign(alkon_kc_hall$NES[match(reynolds_kc_hall$ID, alkon_kc_hall$ID)]), ]
common_kc_hallmarks_contrary

# COMMAND ----------

display(common_kc_hallmarks)

# COMMAND ----------

rows <- min(nrow(reynolds_kc_hall), nrow(alkon_kc_hall))

# COMMAND ----------

percent_similarity_kc <- nrow(common_kc_hallmarks) / rows *100
percent_similarity_kc

# COMMAND ----------

# MAGIC %md
# MAGIC Hallmarks C2

# COMMAND ----------

reynolds_kc_hall2 <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_kc_0.05_reynolds.xlsx", sheet =4)
alkon_kc_hall2 <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_kc_0.05_alkon.xlsx", sheet=4)

# COMMAND ----------

common_kc_hallmarks2 <- reynolds_kc_hall2[reynolds_kc_hall2$ID %in% alkon_kc_hall2$ID & sign(reynolds_kc_hall2$NES) == sign(alkon_kc_hall2$NES[match(reynolds_kc_hall2$ID, alkon_kc_hall2$ID)]), ]
display(common_kc_hallmarks2)

# COMMAND ----------

common_kc_hallmarks2_contrary <- reynolds_kc_hall2[reynolds_kc_hall2$ID %in% alkon_kc_hall2$ID & sign(reynolds_kc_hall2$NES) != sign(alkon_kc_hall2$NES[match(reynolds_kc_hall2$ID, alkon_kc_hall2$ID)]), ]
common_kc_hallmarks2_contrary

# COMMAND ----------

common_kc <- list(
  kegg = NULL,
  Reactome = common_kc_reactome,
  Hallmarks1 = common_kc_hallmarks,
  Hallmarks2 = common_kc_hallmarks2
)

# COMMAND ----------

write.xlsx(common_kc, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/common_GSEA_individually/common_pathways_kc.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tcell

# COMMAND ----------

# MAGIC %md
# MAGIC ###kegg

# COMMAND ----------

reynolds_tcell_kegg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_tcell_0.05_reynolds.xlsx", sheet =1)
alkon_tcell_kegg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_tcell_0.05_alkon.xlsx", sheet=1)

# COMMAND ----------

# MAGIC %md
# MAGIC no terms in alkon

# COMMAND ----------

# MAGIC %md
# MAGIC ###reactome

# COMMAND ----------

reynolds_tcell_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_tcell_0.05_reynolds.xlsx", sheet =2)
alkon_tcell_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_tcell_0.05_alkon.xlsx", sheet=2)

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
# MAGIC ###Hallmark H2 

# COMMAND ----------

reynolds_tcell_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_tcell_0.05_reynolds.xlsx", sheet =3)
alkon_tcell_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_tcell_0.05_alkon.xlsx", sheet=3)

# COMMAND ----------

reynolds_tcell_hall

# COMMAND ----------

alkon_tcell_hall

# COMMAND ----------

# MAGIC %md
# MAGIC No terms are common with same sign. 

# COMMAND ----------

common_tcell_hallmarks_contrary <- reynolds_tcell_hall[reynolds_tcell_hall$ID %in% alkon_tcell_hall$ID & sign(reynolds_tcell_hall$NES) != sign(alkon_tcell_hall$NES[match(reynolds_tcell_hall$ID, alkon_tcell_hall$ID)]), ]
common_tcell_hallmarks_contrary

# COMMAND ----------

# MAGIC %md
# MAGIC HALLMARK_TNFA_SIGNALING_VIA_NFKB is positive in Reynolds and negative enriched in ALkon

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

# MAGIC %md
# MAGIC ###Kegg

# COMMAND ----------

reynolds_fb_kegg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_0.05_reynolds.xlsx", sheet =1)
alkon_fb_kegg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_fb_0.05_alkon.xlsx", sheet=1)

# COMMAND ----------

common_fb_kegg <- reynolds_fb_kegg[reynolds_fb_kegg$ID %in% alkon_fb_kegg$ID & sign(reynolds_fb_kegg$NES) == sign(alkon_fb_kegg$NES[match(reynolds_fb_kegg$ID, alkon_fb_kegg$ID)]), ]
common_fb_kegg

# COMMAND ----------

# MAGIC %md
# MAGIC ###reactome

# COMMAND ----------

reynolds_fb_react <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_0.05_reynolds.xlsx", sheet =2)
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
# MAGIC 53% are common terms with contrary NES sign

# COMMAND ----------

# MAGIC %md
# MAGIC ###Hallmarks H2

# COMMAND ----------

reynolds_fb_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_0.05_reynolds.xlsx", sheet =3)
alkon_fb_hall <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_fb_0.05_alkon.xlsx", sheet=3)

# COMMAND ----------

common_fb_hallmarks <- reynolds_fb_hall[reynolds_fb_hall$ID %in% alkon_fb_hall$ID & sign(reynolds_fb_hall$NES) == sign(alkon_fb_hall$NES[match(reynolds_fb_hall$ID, alkon_fb_hall$ID)]), ]
common_fb_hallmarks

# COMMAND ----------

common_fb_hallmarks_contrary <- reynolds_fb_hall[reynolds_fb_hall$ID %in% alkon_fb_hall$ID & sign(reynolds_fb_hall$NES) != sign(alkon_fb_hall$NES[match(reynolds_fb_hall$ID, alkon_fb_hall$ID)]), ]
common_fb_hallmarks_contrary

# COMMAND ----------

# MAGIC %md
# MAGIC No common but 4 are common with contrary NES sign:
# MAGIC HALLMARK_INTERFERON_GAMMA_RESPONSE Negative in Reynolds
# MAGIC HALLMARK_INTERFERON_ALPHA_RESPONSE Negative in Reynolds
# MAGIC HALLMARK_MYC_TARGETS_V1   Negative in Reynolds         
# MAGIC HALLMARK_MTORC1_SIGNALING   Negative in Reynolds      
# MAGIC HALLMARK_APICAL_JUNCTION Negative in Alkon

# COMMAND ----------

# MAGIC %md
# MAGIC ###Hallmarks C2

# COMMAND ----------

reynolds_fb_hall2 <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_0.05_reynolds.xlsx", sheet =4)
alkon_fb_hall2 <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_fb_0.05_alkon.xlsx", sheet=4)

# COMMAND ----------

common_fb_hallmarks2 <- reynolds_fb_hall2[reynolds_fb_hall2$ID %in% alkon_fb_hall2$ID & sign(reynolds_fb_hall2$NES) == sign(alkon_fb_hall2$NES[match(reynolds_fb_hall2$ID, alkon_fb_hall2$ID)]), ]
common_fb_hallmarks2

# COMMAND ----------

common_fb_hallmarks2_contrary <- reynolds_fb_hall2[reynolds_fb_hall2$ID %in% alkon_fb_hall2$ID & sign(reynolds_fb_hall2$NES) != sign(alkon_fb_hall2$NES[match(reynolds_fb_hall2$ID, alkon_fb_hall2$ID)]), ]
common_fb_hallmarks2_contrary

# COMMAND ----------

rows <- min(nrow(reynolds_fb_hall2), nrow(alkon_fb_hall2))

# COMMAND ----------

percent_similarity_fb <- nrow(common_fb_hallmarks2) / rows *100
percent_similarity_fb

# COMMAND ----------

# MAGIC %md
# MAGIC 29% are common with the same NES sign in hallmarks c2

# COMMAND ----------

percent_similarity_contrary_fb <- nrow(common_fb_hallmarks2_contrary) / rows *100
percent_similarity_contrary_fb

# COMMAND ----------

# MAGIC %md
# MAGIC 33,3% are common with contrary NES sign in hallmarks C2

# COMMAND ----------

common_fb <- list(
  Kegg = common_fb_kegg,
  Reactome = common_fb_reactome,
  Hallmarks1 = common_fb_hallmarks,
  Hallmarks2 = common_fb_hallmarks2
)

# COMMAND ----------

write.xlsx(common_fb, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/common_GSEA_individually/common_pathways_fb.xlsx")
