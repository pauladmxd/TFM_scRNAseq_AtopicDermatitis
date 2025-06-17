# Databricks notebook source
# MAGIC %md
# MAGIC #Find DEGs with pseudobulk (limma)
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
  #remotes::install_cran(pkg_name, lib = temp_library, upgrade=FALSE)
  remotes::install_bioc(pkg_name, lib=temp_library, upgrade=FALSE)
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

if (!requireNamespace("limma"))install_from_github("limma")

# COMMAND ----------

if (!requireNamespace("edgeR"))install_from_github("edgeR")

# COMMAND ----------

.libPaths(c("/dbfs/home/pdelgadom@almirall.com/my_r_packages/tfm_paula_4", .libPaths()))
library(limma)
library(edgeR)

# COMMAND ----------

.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(EnhancedVolcano)
library(VennDiagram)

# COMMAND ----------

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat_v2/", .libPaths()))
library(Seurat)

# COMMAND ----------

volcano_generator<- function(resultsDE){
resultsDE<- as.data.frame(resultsDE)

#Create anotations for volcano plot
resultsDE0 <- resultsDE
rownames(resultsDE0) <- resultsDE0$gene

top10_genes <- resultsDE0 %>% 
filter(logFC > 1 & adj.P.Val < 0.05 ) %>%
arrange(adj.P.Val) %>% top_n(10, -adj.P.Val)

bottom10_genes <-  resultsDE0 %>% 
filter(logFC < -1 & adj.P.Val < 0.05 ) %>%
arrange(adj.P.Val) %>% top_n(10, -adj.P.Val)

#Plot Volcano
volcano <- (EnhancedVolcano(resultsDE0,
lab = rownames(resultsDE0),
x = 'logFC',
y = 'adj.P.Val',
pCutoff = 0.05,
selectLab = c(top10_genes$gene, bottom10_genes$gene),
labSize = 5,
drawConnectors = TRUE,
widthConnectors = 0.5,
colConnectors = 'black'))
volcano

}

# COMMAND ----------

# MAGIC %md
# MAGIC ##Read data

# COMMAND ----------

#Read seurat object
alkon <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")

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

# MAGIC %md
# MAGIC ##Pseudobulk the counts based on the donor id

# COMMAND ----------

# pseudobulk the counts based on donor-condition-celltype
pseudo_alkon <- AggregateExpression(alkon_f, assays = "RNA", return.seurat = T, group.by = c("Condition", "Sample_id", "h_celltype_v4"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_alkon))

# COMMAND ----------

pseudo_alkon

# COMMAND ----------

pseudo_alkon$celltype.cond <- paste(pseudo_alkon$h_celltype_v4, pseudo_alkon$Condition, sep = "_")

# COMMAND ----------

Idents(pseudo_alkon) <- "celltype.cond"

# COMMAND ----------

# MAGIC %md
# MAGIC ##PCA

# COMMAND ----------


# Step 2: Normalize the data
pseudo_alkon <- NormalizeData(pseudo_alkon)

# Step 3: Find variable features
pseudo_alkon <- FindVariableFeatures(pseudo_alkon)

# Step 4: Scale the data
pseudo_alkon <- ScaleData(pseudo_alkon)

# Step 5: Run PCA
pseudo_alkon <- RunPCA(pseudo_alkon, features = VariableFeatures(object = pseudo_alkon))

# COMMAND ----------

# Visualize PCA results with sample ID labels and color by celltype
plot1 <- DimPlot(pseudo_alkon, reduction = "pca", group.by = "h_celltype_v4", label = TRUE, repel = TRUE)

# Visualize PCA results with sample ID labels and color by sample id
plot2 <- DimPlot(pseudo_alkon, reduction = "pca", group.by = "Condition", label = TRUE, repel = TRUE)

options(repr.plot.width=1500, repr.plot.height=1200)

plot1 + plot2

# COMMAND ----------

# MAGIC %md
# MAGIC #Limma for pseudobulk

# COMMAND ----------

expr_matrix <- pseudo_alkon$RNA$counts
metadata <- data.frame(sample_id= colnames(expr_matrix))

# COMMAND ----------

# Split names into components by underscore
metadata$condition <- sapply(strsplit(as.character(metadata$sample_id), "_"), `[`, 1)
metadata$sample <- sapply(strsplit(as.character(metadata$sample_id), "_"), `[`, 2)
metadata$celltype <- sapply(strsplit(as.character(metadata$sample_id), "_"), `[`, 3)
rownames(metadata) <- metadata$sample_id

all(rownames(metadata) == colnames(expr_matrix))  # should be TRUE
 

# COMMAND ----------

# Extract sample names
samples <- colnames(expr_matrix)
 
# Parse to get cell types
celltypes <- sapply(strsplit(samples, "_"), function(x) x[3])
 
# Get unique cell types
unique_celltypes <- unique(celltypes)

# COMMAND ----------

expr_matrix <- as.matrix(expr_matrix)

# COMMAND ----------

wanted_celltypes <- c("NK", "Macro", "TC", "KC", "Fibroblasts", "Treg")
# 'samples' assumed to be: colnames of expr_matrix
samples <- colnames(expr_matrix)
alkon_limma_res <- list()
for (ct in wanted_celltypes) {
  # Subset sample names that end with the current cell type
  cols_ct <- samples[grepl(paste0("_", ct, "$"), samples)]
  if (length(cols_ct) < 2) {
    cat("Skipping cell type:", ct, "- not enough samples\n")
    next
  }
  # Extract metadata from column names
  split_list <- strsplit(cols_ct, "_")
  sample_info_ct <- do.call(rbind, split_list)
  colnames(sample_info_ct) <- c("condition", "sample", "celltype")
  sample_info_ct <- as.data.frame(sample_info_ct, stringsAsFactors = TRUE)
  rownames(sample_info_ct) <- cols_ct
  sample_info_ct$condition <- factor(sample_info_ct$condition, levels = c("HC", "AD"))
  # Subset expression matrix
  expr_ct <- expr_matrix[, cols_ct, drop = FALSE]
  expr_ct <- as.matrix(expr_ct)
  mode(expr_ct) <- "numeric"
  # Create DGEList and filter
  dge_ct <- DGEList(counts = expr_ct)
  keep <- filterByExpr(dge_ct, group = sample_info_ct$condition)
  dge_ct <- dge_ct[keep, , keep.lib.sizes = FALSE]
  dge_ct <- calcNormFactors(dge_ct)
  # Design matrix
  design_ct <- model.matrix(~ condition, data = sample_info_ct)
  # voom and limma
  v_ct <- voom(dge_ct, design_ct, plot = FALSE)
  fit_ct <- lmFit(v_ct, design_ct)
  fit_ct <- eBayes(fit_ct)
  # Extract DE results
  res_ct <- topTable(fit_ct, coef = "conditionAD", number = Inf)
  res_ct$celltype <- ct
  # Store results
  alkon_limma_res[[ct]] <- res_ct
  cat("Completed DE for cell type:", ct, "\n")
}
# Combine all cell type DE results
alkon_limma_res <- bind_rows(alkon_limma_res)
head(alkon_limma_res)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tcells

# COMMAND ----------

TC_results <- alkon_limma_res %>% filter(celltype=="TC")

# COMMAND ----------

TC_results$gene <- rownames(TC_results)

# COMMAND ----------

TC_results$logFC

# COMMAND ----------

class(TC_results$adj.P.Val)

# COMMAND ----------

sum(is.na(TC_results$logFC))

# COMMAND ----------

volcano_generator(TC_results)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

KC_results <-alkon_limma_res %>% filter(celltype=="KC")

# COMMAND ----------

KC_results$gene <- rownames(KC_results)

# COMMAND ----------

KC_results$logFC

# COMMAND ----------

class(KC_results$adj.P.Val)

# COMMAND ----------

sum(is.na(KC_results$logFC))

# COMMAND ----------

volcano_generator(KC_results)

# COMMAND ----------

display(FB_results)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblast

# COMMAND ----------

FB_results <- alkon_limma_res %>% filter(celltype=="Fibroblasts")

# COMMAND ----------

FB_results$gene <- rownames(FB_results)

# COMMAND ----------

volcano_generator(FB_results)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

Treg_results <- alkon_limma_res %>% filter(celltype=="Treg")

# COMMAND ----------

Treg_results$gene <- rownames(Treg_results)

# COMMAND ----------

volcano_generator(Treg_results)

# COMMAND ----------



# COMMAND ----------

# MAGIC %md
# MAGIC ##Macro

# COMMAND ----------

Macro_results <- alkon_limma_res %>% filter(celltype=="Macro")

# COMMAND ----------

Macro_results$gene <- rownames(Macro_results)

# COMMAND ----------

volcano_generator(Macro_results)

# COMMAND ----------

# MAGIC %md
# MAGIC #Save

# COMMAND ----------

# Delete existing Excel files
# file.remove("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_tcell.xlsx")
# file.remove("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_kc.xlsx")
# file.remove("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_FB.xlsx")
# file.remove("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_macro.xlsx")
# file.remove("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_treg.xlsx")

# Write the data frame to new Excel files
write.xlsx(TC_results, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_tcell.xlsx")
write.xlsx(KC_results, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_kc.xlsx")
write.xlsx(FB_results, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_FB.xlsx")
write.xlsx(Macro_results, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_macro.xlsx")
write.xlsx(Treg_results, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_limma_results_treg.xlsx")

# COMMAND ----------

display(TC_results)

# COMMAND ----------

display(KC_results)
