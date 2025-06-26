# Databricks notebook source
# MAGIC %md
# MAGIC # Common markers from the different datasets
# MAGIC
# MAGIC Here I will merge the markers I found for each cell type on each dataset, for the **LvsHC** contrast to see which are common and have the same sign so represent a consistent result.

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(dittoSeq)
library(dplyr)
library(openxlsx)
.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(VennDiagram)

# COMMAND ----------

# MAGIC %md
# MAGIC ## Lesional vs Healthy Control

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

# Read each dataset markers saved previously
reynolds_LvsHC_KC_markers <- reynolds_LvsHC_KC_markers %>%
  filter(avg_log2FC  > 1, p_val_adj < 0.05)

alkon_LvsHC_KC_markers <- alkon_LvsHC_KC_markers %>%
  filter(avg_log2FC  > 1, p_val_adj < 0.05)

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "Alkon" = na.omit(reynolds_LvsHC_KC_markers$gene),
  "Reynolds" = na.omit(alkon_LvsHC_KC_markers$gene)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs btw datasets - KC (FindAllMarkers)", # Add title,
  main.fontface = "bold", # Make title bold
  cat.dist = c(0.03, 0.03), # Adjust the distance of the category names from the circles
  main.cex = 1.3, # Increase title size
  cat.cex = 1.3, # Increase label size
  cat.pos = c(-17, 17), # Position labels more on the top
  cex = 1.3, # Increase numbers size
  cat.col = c("red", "blue") # Set label colors to match circles
)


# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# Filter rows where cluster_x is "AD" and cluster_y is "lesional" at the same time or healthy for both.
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
# MAGIC ### T-cells

# COMMAND ----------

#Read each dataset markers saved previously
reynolds_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_LvsHC_tcell_allmarkers.xlsx")

alkon_LvsHC_tcell_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

#Common markers
tcell_markersRA <- merge(reynolds_LvsHC_tcell_markers, alkon_LvsHC_tcell_markers, by = "gene", suffixes = c(".reynolds", ".alkon"))

# COMMAND ----------

display(tcell_markersRA)

# COMMAND ----------

#Filter common markers

filtered_tcell_markersRA <- tcell_markersRA %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)


# COMMAND ----------

display(filtered_tcell_markersRA)

# COMMAND ----------

write.xlsx(filtered_tcell_markersRA, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/AR_Tcell_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Macrophages

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_allmarkers.xlsx")

alkon_LvsHC_macro_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

alkon_LvsHC_macro_markers_suffix <- alkon_LvsHC_macro_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  

macro_markersAR <- merge(reynolds_LvsHC_macro_markers, alkon_LvsHC_macro_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
display(macro_markersAR %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster_x is "AD" and cluster_y is "lesional" at the same time or healthy for both.
filtered_macro_markersAR <- macro_markersAR %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)


# COMMAND ----------

display(filtered_macro_markersAR)

# COMMAND ----------

write.xlsx(filtered_macro_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/AR_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Treg

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_allmarkers.xlsx")

alkon_LvsHC_treg_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

alkon_LvsHC_treg_markers_suffix <- alkon_LvsHC_treg_markers %>%
  rename_with(~ paste0(., ".alkon"), -gene) #add also alkon name, not possible with suffix because they don't share the same name then
  
treg_markersAR <- merge(reynolds_LvsHC_treg_markers, alkon_LvsHC_treg_markers, by = "gene", suffixes=c(".alkon", ".reynolds"))
display(treg_markersAR %>% arrange(p_val_adj.alkon))

# COMMAND ----------

# Filter rows where cluster_x is "AD" and cluster_y is "lesional" at the same time or healthy for both.
filtered_treg_markersAR <- treg_markersAR %>%
  filter((cluster.reynolds == "lesional" & cluster.alkon == "AD") | 
         (cluster.reynolds == "healthy" & cluster.alkon == "HC")) %>%
  mutate(Condition = ifelse(cluster.alkon == "AD", "lesional", "healthy")) %>%
  mutate(avg_log2FC = (avg_log2FC.alkon + avg_log2FC.reynolds) / 2) %>%
  mutate(avg_pvalue = (p_val.alkon + p_val.reynolds) / 2) %>%
  mutate(avg_pvalue_adj = (p_val_adj.alkon + p_val_adj.reynolds) / 2) %>%
  select(-cluster.reynolds, -cluster.alkon)


display(filtered_treg_markersAR %>% arrange(avg_pvalue_adj))

# COMMAND ----------

write.xlsx(filtered_treg_markersAR, file="/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/AR_treg_LvsHC_allmarkers.xlsx")
