# Databricks notebook source
# MAGIC %md
# MAGIC #Find DEGs with pseudobulk
# MAGIC ## AD/PS Reynolds et al, 2021
# MAGIC
# MAGIC Here will be performed the analysis to find DEGs in the subannotation of KC only with one dataset to compare with the merged analysis results.

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

#Load required libraries
.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(VennDiagram)
library(EnhancedVolcano)

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat_v2", .libPaths()))
library(Seurat)

# COMMAND ----------

# MAGIC %md
# MAGIC The packages have to be loaded in this specific order, if not an error arises

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
reynolds <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

head(reynolds@meta.data, 5)

# COMMAND ----------

# Create a new column 'Condition' with default value 'healthy'
reynolds$Condition <- "healthy"

# Update 'Condition' based on the 'Status' and 'Site' columns
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "lesion"] <- "lesional"
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "non_lesion"] <- "non lesional"

unique(reynolds$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC ###Filter to remove what is not relevant
# MAGIC - non lesional samples 
# MAGIC - all celltypes except to Tc, Th, Treg, Undifferentiated_KC, Differentiated_KC
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC Filter cell types, in that way I can run DESEQ2 faster

# COMMAND ----------

desired_cell_types <- c("Tc", "Th", "Treg", "Undifferentiated_KC*", "Differentiated_KC")  # Relevant celltypes
reynolds_f <- subset(reynolds, final_clustering %in% desired_cell_types)

# COMMAND ----------

unique(reynolds_f$final_clustering)

# COMMAND ----------

reynolds_f$final_clustering <- ifelse(reynolds_f$final_clustering == "Undifferentiated_KC*", "Undifferentiated_KC", reynolds_f$final_clustering)

# COMMAND ----------

unique(reynolds_f$final_clustering)

# COMMAND ----------

table(reynolds_f$Condition, reynolds_f$final_clustering)

# COMMAND ----------

unique(reynolds_f$Condition)

# COMMAND ----------

reynolds_f <- subset(reynolds_f, Condition %in% c("lesional", "healthy"))

# COMMAND ----------

unique(reynolds_f$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Filtering variables that have at least 3 counts
# MAGIC

# COMMAND ----------

counts_matrix <- reynolds_f[["RNA"]]$counts
dim(counts_matrix)

# COMMAND ----------

# Keep only rows that have a count of at least 3 counts in 3 samples
smallestGroupSize <- 3
keep <- rowSums(counts_matrix >= 3) >= smallestGroupSize
counts_keep <- counts_matrix[keep,]

# Subset the Seurat object to keep only the features in counts_keep
reynolds_f <- subset(reynolds_f, features = rownames(counts_keep))

# Assign the filtered counts to the new Seurat object
reynolds_f[["RNA"]]$counts <- counts_keep

# Check dimensions
dim(reynolds_f[["RNA"]]$counts)

# COMMAND ----------

head(reynolds_f[["RNA"]]$counts)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Pseudobulk the counts based on the donor id

# COMMAND ----------

# pseudobulk the counts based on donor-condition-celltype
pseudo_reynolds <- AggregateExpression(reynolds_f, assays = "RNA", return.seurat = T, group.by = c("Condition", "donor_id", "final_clustering"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
unique(Cells(pseudo_reynolds))

# COMMAND ----------

pseudo_reynolds

# COMMAND ----------

pseudo_reynolds$celltype.cond <- paste(pseudo_reynolds$final_clustering, pseudo_reynolds$Condition, sep = "_")

# COMMAND ----------

Idents(pseudo_reynolds) <- "celltype.cond"

# COMMAND ----------

unique(pseudo_reynolds$celltype.cond)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (Undifferentiated)

# COMMAND ----------

bulk.kc.de <- FindMarkers(object = pseudo_reynolds, 
                         ident.1 = "Undifferentiated-KC_lesional", 
                         ident.2 = "Undifferentiated-KC_healthy",
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

volcano_generator(bulk.kc.de, "Pseudobulk - Undifferentiated KC - Reynolds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (Differentiated)

# COMMAND ----------

bulk.dif.kc.de <- FindMarkers(object = pseudo_reynolds, 
                         ident.1 = "Differentiated-KC_lesional", 
                         ident.2 = "Differentiated-KC_healthy",
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

volcano_generator(bulk.dif.kc.de, "Pseudobulk - Differentiated KC - Reynolds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Write and save

# COMMAND ----------

bulk.kc.de$gene <- rownames(bulk.kc.de)
bulk.dif.kc.de$gene <- rownames(bulk.dif.kc.de)

write.xlsx(bulk.kc.de, file="/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_reynolds_bulk_undif_kc_de.xlsx")
write.xlsx(bulk.dif.kc.de, file="/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_reynolds_bulk_dif_kc_de.xlsx")
