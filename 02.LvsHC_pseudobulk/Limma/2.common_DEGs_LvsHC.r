# Databricks notebook source
# MAGIC %md
# MAGIC # Common DEGs from the different datasets
# MAGIC
# MAGIC Here I will merge the DEGs I found for each cell type on each dataset using the **pseudobulk** approach, for the **LvsHC** contrast

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
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
reynolds_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/pseudobulk/reynolds_limma_results_kc.xlsx")
alkon_LvsHC_KC_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_kc.xlsx")

rownames(reynolds_LvsHC_KC_markers) <- reynolds_LvsHC_KC_markers$gene
rownames(alkon_LvsHC_KC_markers) <- alkon_LvsHC_KC_markers$gene

# COMMAND ----------

# Read each dataset markers saved previously
reynolds_LvsHC_KC_markers <- reynolds_LvsHC_KC_markers %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

alkon_LvsHC_KC_markers <- alkon_LvsHC_KC_markers %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

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
  main = "Common DEGs btw datasets - KC (Limma)", # Add title,
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
