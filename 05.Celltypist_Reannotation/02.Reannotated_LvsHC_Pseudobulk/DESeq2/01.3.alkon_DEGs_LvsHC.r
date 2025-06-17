# Databricks notebook source
# MAGIC %md
# MAGIC #Find DEGs with pseudobulk
# MAGIC ##AD/PN – Alkon et al, 2023
# MAGIC - Does not include non-lesional samples.
# MAGIC - 5 AD, 7 PN, 1 AP (atopic prurigo) and 4 healthy control​​
# MAGIC
# MAGIC Here will be performed the analysis to find DEGs of each relevant cell type in Alkon et al, 2023 dataset following: https://satijalab.org/seurat/articles/de_vignette#perform-de-analysis-after-pseudobulking).
# MAGIC
# MAGIC ###Most relevant cell types: 
# MAGIC T-cells (TC), Fibroblasts, Keratinocytes (KC),  Monocytes, Macrophages, Dendritic cells, Natural killers, Treg and MastC
# MAGIC
# MAGIC ###Constrast: 
# MAGIC   - Lesional vs Healthy control (LvsHC)

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(dplyr)
library(openxlsx)

# COMMAND ----------

my_library <- "/dbfs/home/pdelgadom@almirall.com/my_r_packages/tfm_paula_4"
dir.create(my_library, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(my_library, .libPaths()))
if (!requireNamespace("remotes")) install.packages("remotes")
 
install_from_github <- function(pkg_name, my_library=NULL) {
  if (is.null(my_library)) {
    my_library <- .libPaths()[1]
    message("Installing ", pkg_name, " to ", my_library)
  }
 
  temp_library <- tempfile()
  dir.create(temp_library)
  remotes::install_cran(pkg_name, lib = temp_library, upgrade=FALSE)
  #remotes::install_bioc(pkg_name, lib=temp_library, upgrade=FALSE)
  #remotes::install_github(pkg_name, lib = temp_library, upgrade=FALSE)
  for (x in list.files(temp_library)) {
    file.copy(
      file.path(temp_library, x),
      my_library,
      recursive=TRUE
    )
  }
}

# COMMAND ----------

if (!requireNamespace("DESeq2"))install_from_github("DESeq2")

# COMMAND ----------

.libPaths(c("/dbfs/home/pdelgadom@almirall.com/my_r_packages/tfm_paula_4", .libPaths()))
library(DESeq2)

# COMMAND ----------

.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(VennDiagram)
library(EnhancedVolcano)

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat_v2", .libPaths()))
library(Seurat)

# COMMAND ----------

volcano_generator <- function(resultsDE, given_title) {
  library(dplyr)

  resultsDE <- as.data.frame(resultsDE)
  
  # Create annotations for volcano plot
  resultsDE0 <- resultsDE
  resultsDE0$gene_id <- rownames(resultsDE0)

  # Ensure unique row names and remove rows with missing gene id
  resultsDE0 <- resultsDE0 %>%
    distinct(gene_id, .keep_all = TRUE)
  rownames(resultsDE0) <- resultsDE0$gene_id
  
  # Determine column names for p-value and log2 fold change
  p_val_col <- if ("p_val_adj" %in% colnames(resultsDE0)) "p_val_adj" else "padj"
  log2fc_col <- if ("avg_log2FC" %in% colnames(resultsDE0)) "avg_log2FC" else "log2FoldChange"
  
  top10_genes <- resultsDE0 %>%
    filter(!!sym(log2fc_col) > 1 & !!sym(p_val_col) < 0.05) %>%
    arrange(!!sym(p_val_col)) %>% top_n(10, -!!sym(p_val_col))
  
  bottom10_genes <- resultsDE0 %>%
    filter(!!sym(log2fc_col) < -1 & !!sym(p_val_col) < 0.05) %>%
    arrange(!!sym(p_val_col)) %>% top_n(10, -!!sym(p_val_col))
  
  # Plot Volcano
  volcano <- EnhancedVolcano(resultsDE0,
    lab = rownames(resultsDE0),
    x = log2fc_col,
    y = p_val_col,
    pCutoff = 0.05,
    selectLab = c(top10_genes$gene_id, bottom10_genes$gene_id),
    labSize = 5,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'black',
    title = given_title)
  volcano
}

# COMMAND ----------

# MAGIC %md
# MAGIC ##Read data

# COMMAND ----------

#Read seurat object
alkon <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/ALKON_CELLTYPIST_TFM.rds") #new annotation

# COMMAND ----------

# alkon$Condition <- ifelse(alkon$Condition == "AD", "Lesional", alkon$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Filtering variables that have at least 3 counts

# COMMAND ----------

counts_matrix <- alkon[["RNA"]]$counts
dim(counts_matrix)

# COMMAND ----------

# Keep only rows that have a count of at least 3 counts in 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts_matrix >= 3) >= smallestGroupSize
counts_keep <- counts_matrix[keep,]

# Subset the Seurat object to keep only the features in counts_keep
alkon_f <- subset(alkon, features = rownames(counts_keep))

# Assign the filtered counts to the new Seurat object
alkon_f[["RNA"]]$counts <- counts_keep

# Check dimensions
dim(alkon_f[["RNA"]]$counts)

# COMMAND ----------

head(alkon_f[["RNA"]]$counts)

# COMMAND ----------

unique(alkon_f$celltypist)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Pseudobulk the counts based on the donor id

# COMMAND ----------

# pseudobulk the counts based on donor-condition-celltype
pseudo_alkon <- AggregateExpression(alkon_f, assays = "RNA", return.seurat = T, group.by = c("Condition", "Sample_id", "celltypist"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_alkon))

# COMMAND ----------

pseudo_alkon

# COMMAND ----------

pseudo_alkon$celltype.cond <- paste(pseudo_alkon$celltypist, pseudo_alkon$Condition, sep = "_")

# COMMAND ----------

Idents(pseudo_alkon) <- "celltype.cond"

# COMMAND ----------

unique(pseudo_alkon$celltype.cond)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (Undifferentiated)

# COMMAND ----------

bulk.kc.de <- FindMarkers(object = pseudo_alkon, 
                         ident.1 = "Undifferentiated-KC_AD", 
                         ident.2 = "Undifferentiated-KC_HC",
                         test.use = "DESeq2")
head(bulk.kc.de, n = 15)

# COMMAND ----------

library(ggplot2)
bulk.kc.de <- as.data.frame(bulk.kc.de)
ggplot(bulk.kc.de, aes(x = p_val)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of p-value", x = "p-value", y = "Frequency")

# COMMAND ----------

volcano_generator(bulk.kc.de, "Pseudobulk - Undifferentiated KC - Alkon")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (Differentiated)

# COMMAND ----------

bulk.dif.kc.de <- FindMarkers(object = pseudo_alkon, 
                         ident.1 = "Differentiated-KC_AD", 
                         ident.2 = "Differentiated-KC_HC",
                         test.use = "DESeq2")
head(bulk.dif.kc.de, n = 15)

# COMMAND ----------

library(ggplot2)
bulk.kc.de <- as.data.frame(bulk.dif.kc.de)
ggplot(bulk.kc.de, aes(x = p_val)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of p-value", x = "p-value", y = "Frequency")

# COMMAND ----------

volcano_generator(bulk.dif.kc.de, "Pseudobulk - Differentiated KC - Alkon")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Write and save

# COMMAND ----------

bulk.kc.de$gene <- rownames(bulk.kc.de)
bulk.dif.kc.de$gene <- rownames(bulk.dif.kc.de)

write.xlsx(bulk.kc.de, file="/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_alkon_bulk_undif_kc_de.xlsx")
write.xlsx(bulk.dif.kc.de, file="/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_alkon_bulk_dif_kc_de.xlsx")
