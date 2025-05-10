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

#Load required libraries
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)

.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))
library(VennDiagram)
library(EnhancedVolcano)

# COMMAND ----------

volcano_generator<- function(resultsDE){
resultsDE<- as.data.frame(resultsDE)

#Create anotations for volcano plot
resultsDE0 <- resultsDE
rownames(resultsDE0) <- resultsDE0$gene

top10_genes <- resultsDE0 %>% 
filter(avg_log2FC > 2 & p_val_adj < 0.05 ) %>%
arrange(p_val_adj) %>% top_n(10, -p_val_adj)

bottom10_genes <-  resultsDE0 %>% 
filter(avg_log2FC < -1 & p_val_adj < 0.05 ) %>%
arrange(p_val_adj) %>% top_n(10, -p_val_adj)

#Plot Volcano
volcano <- (EnhancedVolcano(resultsDE0,
lab = rownames(resultsDE0),
x = 'avg_log2FC',
y = 'p_val_adj',
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
# MAGIC saveRDS for Aitor

# COMMAND ----------

# saveRDS(pseudo_alkon, file="/dbfs/mnt/sandbox/TFM_PAULA/ALKON_aggregated_expression_TFM.rds")

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
# MAGIC ##Tcells

# COMMAND ----------

bulk.tcell.de <- FindMarkers(object = pseudo_alkon, 
                         ident.1 = "TC_AD", 
                         ident.2 = "TC_HC",
                         test.use = "DESeq2")
head(bulk.tcell.de, n = 15)


# COMMAND ----------

bulk.tcell.de.orig <- bulk.tcell.de
bulk.tcell.de.orig$gene <- rownames(bulk.tcell.de.orig)

# COMMAND ----------

tcell.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_tcell_LvsHC_allmarkers.xlsx")

# COMMAND ----------

rownames(tcell.de) <- tcell.de$gene

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.tcell.de) <- paste0(names(bulk.tcell.de), ".bulk")
bulk.tcell.de$gene <- rownames(bulk.tcell.de)

names(tcell.de) <- paste0(names(tcell.de), ".sc")
tcell.de$gene <- rownames(tcell.de)

merge_dat <- merge(tcell.de, bulk.tcell.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val_adj.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05 & 
                                merge_dat$avg_log2FC.bulk > 0.5 & 
                                merge_dat$avg_log2FC.sc > 0.5)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05 & 
                                  merge_dat$avg_log2FC.sc > 0.5)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05 & 
                                    merge_dat$avg_log2FC.bulk > 0.5)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))

# COMMAND ----------

tcell.de$p_val_adj.sc

# COMMAND ----------

tcell.de.sig <- tcell.de[tcell.de$p_val_adj.sc < 0.05 & tcell.de$avg_log2FC.sc > 0.5,]
bulk.tcell.de.sig <- bulk.tcell.de[bulk.tcell.de$p_val_adj.bulk < 0.05 & bulk.tcell.de$avg_log2FC.bulk > 0.5,]

# COMMAND ----------

common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05 & 
                                merge_dat$avg_log2FC.bulk > 0.5 & 
                                merge_dat$avg_log2FC.sc > 0.5)]
common <- common[order(merge_dat$p_val_adj.bulk[match(common, merge_dat$gene)])]
print(common[1:20])

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(tcell.de.sig),
  "Bulk" = rownames(bulk.tcell.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Alkon Tcells", # Add title
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
# MAGIC ##Fibroblasts

# COMMAND ----------

unique(alkon$h_celltype_v4)

# COMMAND ----------

bulk.fb.de <- FindMarkers(object = pseudo_alkon, 
                         ident.1 = "Fibroblasts_AD", 
                         ident.2 = "Fibroblasts_HC",
                         test.use = "DESeq2")
head(bulk.fb.de, n = 15)

# COMMAND ----------

bulk.fb.de.orig <- bulk.fb.de
bulk.fb.de.orig$gene <- rownames(bulk.fb.de.orig)

# COMMAND ----------

fb.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_allmarkers.xlsx")
rownames(fb.de) <- fb.de$gene
fb.de.orig <- fb.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.fb.de) <- paste0(names(bulk.fb.de), ".bulk")
bulk.fb.de$gene <- rownames(bulk.fb.de)

names(fb.de) <- paste0(names(fb.de), ".sc")
fb.de$gene <- rownames(fb.de)

merge_dat <- merge(fb.de, bulk.fb.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val_adj.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))

# COMMAND ----------

fb.de.sig <- fb.de[fb.de$p_val_adj.sc < 0.05 & fb.de$avg_log2FC.sc > 1,]
bulk.fb.de.sig <- bulk.fb.de[bulk.fb.de$p_val_adj.bulk < 0.05 & bulk.fb.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(fb.de.sig),
  "Bulk" = rownames(bulk.fb.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Alkon Fibroblasts", # Add title
  cat.pos = c(0, 0), # Move the name of each circle higher
  cat.dist = c(0.02, 0.02) # Adjust the distance of the category names from the circles
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

volcano_generator(bulk.fb.de.orig)

# COMMAND ----------

fb.de.orig <- fb.de.orig %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(fb.de.orig)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

bulk.kc.de <- FindMarkers(object = pseudo_alkon, 
                         ident.1 = "KC_AD", 
                         ident.2 = "KC_HC",
                         test.use = "DESeq2")
head(bulk.kc.de, n = 15)

# COMMAND ----------

bulk.kc.de.orig <- bulk.kc.de
bulk.kc.de.orig$gene <- rownames(bulk.kc.de.orig)

# COMMAND ----------

kc.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_kc_LvsHC_allmarkers.xlsx")
rownames(kc.de) <- kc.de$gene
kc.de.orig <- kc.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.kc.de) <- paste0(names(bulk.kc.de), ".bulk")
bulk.kc.de$gene <- rownames(bulk.kc.de)

names(kc.de) <- paste0(names(kc.de), ".sc")
kc.de$gene <- rownames(kc.de)

merge_dat <- merge(kc.de, bulk.kc.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val_adj.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))

# COMMAND ----------

common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05 & 
                                merge_dat$avg_log2FC.bulk > 0.5 & 
                                merge_dat$avg_log2FC.sc > 0.5)]
common <- common[order(merge_dat$p_val_adj.bulk[match(common, merge_dat$gene)])]
print(common[1:20])

# COMMAND ----------

kc.de.sig <- kc.de[kc.de$p_val_adj.sc < 0.05 & kc.de$avg_log2FC.sc > 0.5,]
bulk.kc.de.sig <- bulk.kc.de[bulk.kc.de$p_val_adj.bulk < 0.05 & bulk.kc.de$avg_log2FC.bulk > 0.5,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(kc.de.sig),
  "Bulk" = rownames(bulk.kc.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Alkon Keratinocytes",,
  cat.dist = c(0.02, 0.02), # Adjust the distance of the category names from the circles
  main.cex = 2, # Increase title size
  cat.cex = 1.5, # Increase label size
  cat.pos = c(-20, 20), # Position labels more on the top
  cex = 1.5 # Increase numbers size
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

volcano_generator(bulk.kc.de.orig)

# COMMAND ----------

kc.de.orig <- kc.de.orig %>%
    mutate(avg_log2FC = ifelse(cluster == "HC", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign

# COMMAND ----------

volcano_generator(kc.de.orig)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

bulk.treg.de <- FindMarkers(object = pseudo_alkon, 
                            ident.1 = "Treg_AD", 
                            ident.2 = "Treg_HC",
                            test.use = "DESeq2")
head(bulk.treg.de, n = 15)

# COMMAND ----------

bulk.treg.de.orig <- bulk.treg.de
bulk.treg.de.orig$gene <- rownames(bulk.treg.de.orig)

# COMMAND ----------

treg.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_allmarkers.xlsx")
rownames(treg.de) <- treg.de$gene
treg.de.orig <- treg.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.treg.de) <- paste0(names(bulk.treg.de), ".bulk")
bulk.treg.de$gene <- rownames(bulk.treg.de)

names(treg.de) <- paste0(names(treg.de), ".sc")
treg.de$gene <- rownames(treg.de)

merge_dat <- merge(treg.de, bulk.treg.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val_adj.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))

# COMMAND ----------

treg.de.sig <- treg.de[treg.de$p_val_adj.sc < 0.05 & treg.de$avg_log2FC.sc > 1,]
bulk.treg.de.sig <- bulk.treg.de[bulk.treg.de$p_val_adj.bulk < 0.05 & bulk.treg.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(treg.de.sig),
  "Bulk" = rownames(bulk.treg.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Alkon Tregs" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

bulk.macro.de <- FindMarkers(object = pseudo_alkon, 
                             ident.1 = "Macro_AD", 
                             ident.2 = "Macro_HC",
                             test.use = "DESeq2")
head(bulk.macro.de, n = 15)

# COMMAND ----------

bulk.macro.de.orig <- bulk.macro.de
bulk.macro.de.orig$gene <- rownames(bulk.macro.de.orig)

# COMMAND ----------

macro.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_macro_LvsHC_allmarkers.xlsx")
rownames(macro.de) <- macro.de$gene
macro.de.orig <- macro.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.macro.de) <- paste0(names(bulk.macro.de), ".bulk")
bulk.macro.de$gene <- rownames(bulk.macro.de)

names(macro.de) <- paste0(names(macro.de), ".sc")
macro.de$gene <- rownames(macro.de)

merge_dat <- merge(macro.de, bulk.macro.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val_adj.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                merge_dat$p_val_adj.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))

# COMMAND ----------

macro.de.sig <- macro.de[macro.de$p_val_adj.sc < 0.05 & macro.de$avg_log2FC.sc > 1,]
bulk.macro.de.sig <- bulk.macro.de[bulk.macro.de$p_val_adj.bulk < 0.05 & bulk.macro.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(macro.de.sig),
  "Bulk" = rownames(bulk.macro.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Alkon Macrophages", # Add title
  cat.pos = c(0, 0), 
  cat.dist = c(0.02, 0.02) # Adjust the distance of the category names from the circles
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC The rest of the cell types are not interesting or have very low number of cells

# COMMAND ----------



# COMMAND ----------

# MAGIC %md
# MAGIC ##Write and save

# COMMAND ----------

write.xlsx(bulk.fb.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_fb_LvsHC_bulk_v2.xlsx")
write.xlsx(bulk.treg.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_treg_LvsHC_bulk_v2.xlsx")
write.xlsx(bulk.kc.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_kc_LvsHC_bulk_v2.xlsx")
write.xlsx(bulk.tcell.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_tcell_LvsHC_bulk_v2.xlsx")
write.xlsx(bulk.macro.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_macro_LvsHC_bulk_v2.xlsx")
