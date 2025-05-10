# Databricks notebook source
# MAGIC %md
# MAGIC #Find DEGs with pseudobulk
# MAGIC ## Liu et al, 2021
# MAGIC - Includes lesional and non-lesional samples
# MAGIC - It only has immune cells
# MAGIC - 4 AD patients, 3 psoriasis and 5 healthy controls.
# MAGIC
# MAGIC Here will be performed the analysis to find all markers of each relevant cell type in Liu et al, 2021 dataset following: https://satijalab.org/seurat/articles/de_vignette#perform-de-analysis-after-pseudobulking).
# MAGIC
# MAGIC ###Most relevant cell types: 
# MAGIC T-cells (TC), T-reg, Innate linfoid cells (ILC), Monocytes (Mono), Macrophagues (Macro), Dendritic cells (DC), Natural killers (NK), MastC
# MAGIC
# MAGIC ###Constrasts: 
# MAGIC   - Lesional vs Healthy control (LvsHC)
# MAGIC   - Non lesional vs Healthy control (NLvsHC)
# MAGIC   - Lesional vs Non-lesional (LvsNL)

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

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(Seurat)

# COMMAND ----------

# MAGIC %md
# MAGIC The packages have to be loaded in this specific order, if not an error arises

# COMMAND ----------

#Load required libraries
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
liu <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/LIU_PROCESSED_TFM.rds")

# COMMAND ----------

head(liu@meta.data, 5)

# COMMAND ----------

# Create a new column 'Condition' with default value 'healthy'
liu$Condition <- "healthy"

# Update 'Condition' based on the 'Status' and 'Site' columns
liu$Condition[liu$Status == "Eczema" & liu$Site == "lesion"] <- "lesional"
liu$Condition[liu$Status == "Eczema" & liu$Site == "non_lesion"] <- "non lesional"

unique(liu$Condition)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Pseudobulk the counts based on the donor id

# COMMAND ----------

# pseudobulk the counts based on donor-condition-celltype
pseudo_liu <- AggregateExpression(liu, assays = "RNA", return.seurat = T, group.by = c("Condition", "donor_id", "h_celltype"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_liu))

# COMMAND ----------

pseudo_liu

# COMMAND ----------

pseudo_liu$celltype.cond <- paste(pseudo_liu$h_celltype, pseudo_liu$Condition, sep = "_")

# COMMAND ----------

Idents(pseudo_liu) <- "celltype.cond"

# COMMAND ----------

# MAGIC %md
# MAGIC ##  PCA

# COMMAND ----------

# Step 2: Normalize the data
pseudo_liu <- NormalizeData(pseudo_liu)

# Step 3: Find variable features
pseudo_liu <- FindVariableFeatures(pseudo_liu)

# Step 4: Scale the data
pseudo_liu <- ScaleData(pseudo_liu)

# Step 5: Run PCA
pseudo_liu <- RunPCA(pseudo_liu, features = VariableFeatures(object = pseudo_liu))

# COMMAND ----------

# Visualize PCA results with sample ID labels and color by celltype
plot1 <- DimPlot(pseudo_liu, reduction = "pca", group.by = "h_celltype", label = TRUE, repel = TRUE)

# Visualize PCA results with sample ID labels and color by sample id
plot2 <- DimPlot(pseudo_liu, reduction = "pca", group.by = "Condition", label = TRUE, repel = TRUE)

options(repr.plot.width=1800, repr.plot.height=1200)

plot1 + plot2

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tcells

# COMMAND ----------

bulk.tcell.de <- FindMarkers(object = pseudo_liu, 
                         ident.1 = "TC_lesional", 
                         ident.2 = "TC_healthy",
                         test.use = "DESeq2")
head(bulk.tcell.de, n = 15)


# COMMAND ----------

bulk.tcell.de.orig <- bulk.tcell.de
bulk.tcell.de.orig$gene <- rownames(bulk.tcell.de.orig)

# COMMAND ----------

tcell.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_LvsHC_tcell_allmarkers.xlsx")

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
                                merge_dat$p_val_adj.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val_adj.bulk > 0.05 & 
                                  merge_dat$p_val_adj.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val_adj.bulk < 0.05 & 
                                    merge_dat$p_val_adj.sc > 0.05)]
print(paste0('# Common: ', length(common)))
print(paste0('# Only in bulk: ', length(only_bulk)))
print(paste0('# Only in single-cell: ', length(only_sc)))
print(common)

# COMMAND ----------

tcell.de.sig <- tcell.de[tcell.de$p_val_adj.sc < 0.05 & tcell.de$avg_log2FC.sc > 1,]
bulk.tcell.de.sig <- bulk.tcell.de[bulk.tcell.de$p_val_adj.bulk < 0.05 & bulk.tcell.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

display(tcell.de.sig)

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
  main = "Common DEGs Liu Tcells" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

bulk.treg.de <- FindMarkers(object = pseudo_liu, 
                            ident.1 = "Treg_lesional", 
                            ident.2 = "Treg_healthy",
                            test.use = "DESeq2")
head(bulk.treg.de, n = 15)

# COMMAND ----------

bulk.treg.de.orig <- bulk.treg.de
bulk.treg.de.orig$gene <- rownames(bulk.treg.de.orig)

# COMMAND ----------

treg.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_treg_LvsHC_allmarkers.xlsx")
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
print(common)
print(only_bulk)

# COMMAND ----------

bulk.treg.de.sig <- bulk.treg.de[bulk.treg.de$p_val_adj.bulk < 0.05,]

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
  main = "Common DEGs Liu Tregs" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

bulk.macro.de <- FindMarkers(object = pseudo_liu, 
                             ident.1 = "Macro_lesional", 
                             ident.2 = "Macro_healthy",
                             test.use = "DESeq2")
head(bulk.macro.de, n = 15)

# COMMAND ----------

bulk.macro.de.orig <- bulk.macro.de
bulk.macro.de.orig$gene <- rownames(bulk.macro.de.orig)

# COMMAND ----------

macro.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_macro_LvsHC_allmarkers.xlsx")
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
print(common)
print(only_bulk)

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
  main = "Common DEGs Liu Macrophages" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

bulk.mono.de <- FindMarkers(object = pseudo_liu, 
                             ident.1 = "Mono_lesional", 
                             ident.2 = "Mono_healthy",
                             test.use = "DESeq2")
head(bulk.mono.de, n = 15)

# COMMAND ----------

bulk.mono.de.orig <- bulk.mono.de
bulk.mono.de.orig$gene <- rownames(bulk.mono.de.orig)

# COMMAND ----------

mono.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mono_LvsHC_allmarkers.xlsx")
rownames(mono.de) <- mono.de$gene
mono.de.orig <- mono.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.mono.de) <- paste0(names(bulk.mono.de), ".bulk")
bulk.mono.de$gene <- rownames(bulk.mono.de)

names(mono.de) <- paste0(names(mono.de), ".sc")
mono.de$gene <- rownames(mono.de)

merge_dat <- merge(mono.de, bulk.mono.de, by = "gene")
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
print(common)

# COMMAND ----------

mono.de.sig <- mono.de[mono.de$p_val_adj.sc < 0.05 & mono.de$avg_log2FC.sc > 1,]
bulk.mono.de.sig <- bulk.mono.de[bulk.mono.de$p_val_adj.bulk < 0.05 & bulk.mono.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(mono.de.sig),
  "Bulk" = rownames(bulk.mono.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Liu Monocytes" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Dendritic cells

# COMMAND ----------

bulk.dc.de <- FindMarkers(object = pseudo_liu, 
                         ident.1 = "DC_lesional", 
                         ident.2 = "DC_healthy",
                         test.use = "DESeq2")
head(bulk.dc.de, n = 15)

# COMMAND ----------

bulk.dc.de.orig <- bulk.dc.de
bulk.dc.de.orig$gene <- rownames(bulk.dc.de.orig)

# COMMAND ----------

dc.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_dc_LvsHC_allmarkers.xlsx")
rownames(dc.de) <- dc.de$gene
dc.de.orig <- dc.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.dc.de) <- paste0(names(bulk.dc.de), ".bulk")
bulk.dc.de$gene <- rownames(bulk.dc.de)

names(dc.de) <- paste0(names(dc.de), ".sc")
dc.de$gene <- rownames(dc.de)

merge_dat <- merge(dc.de, bulk.dc.de, by = "gene")
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
print(common)
print(only_bulk)

# COMMAND ----------

dc.de.sig <- dc.de[dc.de$p_val_adj.sc < 0.05 & dc.de$avg_log2FC.sc > 1,]
bulk.dc.de.sig <- bulk.dc.de[bulk.dc.de$p_val_adj.bulk < 0.05 & bulk.dc.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(dc.de.sig),
  "Bulk" = rownames(bulk.dc.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Liu Dendritic Cells", # Add title,
  cat.pos = c(0, 0), # Move the name of each circle higher
  cat.dist = c(0.02, 0.02) # Adjust the distance of the category names from the circles
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

bulk.ilc.de <- FindMarkers(object = pseudo_liu, 
                          ident.1 = "ILC_lesional", 
                          ident.2 = "ILC_healthy",
                          test.use = "DESeq2")
head(bulk.ilc.de, n = 15)

# COMMAND ----------

bulk.ilc.de.orig <- bulk.ilc.de
bulk.ilc.de.orig$gene <- rownames(bulk.ilc.de.orig)

# COMMAND ----------

ilc.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_ILC_LvsHC_allmarkers.xlsx")
rownames(ilc.de) <- ilc.de$gene
ilc.de.orig <- ilc.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.ilc.de) <- paste0(names(bulk.ilc.de), ".bulk")
bulk.ilc.de$gene <- rownames(bulk.ilc.de)

names(ilc.de) <- paste0(names(ilc.de), ".sc")
ilc.de$gene <- rownames(ilc.de)

merge_dat <- merge(ilc.de, bulk.ilc.de, by = "gene")
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
print(common)
print(only_bulk)

# COMMAND ----------

ilc.de.sig <- ilc.de[ilc.de$p_val_adj.sc < 0.05 & ilc.de$avg_log2FC.sc > 1,]
bulk.ilc.de.sig <- bulk.ilc.de[bulk.ilc.de$p_val_adj.bulk < 0.05 & bulk.ilc.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(ilc.de.sig),
  "Bulk" = rownames(bulk.ilc.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Liu ILCs", # Add title,
  cat.pos = c(0, 0), # Move the name of each circle higher
  cat.dist = c(0.02, 0.02) # Adjust the distance of the category names from the circles
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

bulk.mastc.de <- FindMarkers(object = pseudo_liu, 
                             ident.1 = "MastC_lesional", 
                             ident.2 = "MastC_healthy",
                             test.use = "DESeq2")
head(bulk.mastc.de, n = 15)

# COMMAND ----------

bulk.mastc.de.orig <- bulk.mastc.de
bulk.mastc.de.orig$gene <- rownames(bulk.mastc.de.orig)

# COMMAND ----------

mastc.de <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/liu_mast_LvsHC_allmarkers.xlsx")
rownames(mastc.de) <- mastc.de$gene
mastc.de.orig <- mastc.de

# COMMAND ----------

# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.mastc.de) <- paste0(names(bulk.mastc.de), ".bulk")
bulk.mastc.de$gene <- rownames(bulk.mastc.de)

names(mastc.de) <- paste0(names(mastc.de), ".sc")
mastc.de$gene <- rownames(mastc.de)

merge_dat <- merge(mastc.de, bulk.mastc.de, by = "gene")
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
print(common)
print(only_bulk)

# COMMAND ----------

mastc.de.sig <- mastc.de[mastc.de$p_val_adj.sc < 0.05 & mastc.de$avg_log2FC.sc > 1,]
bulk.mastc.de.sig <- bulk.mastc.de[bulk.mastc.de$p_val_adj.bulk < 0.05 & bulk.mastc.de$avg_log2FC.bulk > 1,]

# COMMAND ----------

# Create a list of your gene sets
gene_sets <- list(
  "SC" = rownames(mastc.de.sig),
  "Bulk" = rownames(bulk.mastc.de.sig)
)

# Plot the Venn diagram with colors and title
venn.plot <- venn.diagram(
  x = gene_sets,
  category.names = c("SC", "Bulk"),
  filename = NULL,  # Set to NULL to plot in RStudio
  output = TRUE,
  fill = c("red", "blue"), # Add colors
  main = "Common DEGs Liu Mastc" # Add title
)

# Display the plot
grid.newpage()
grid.draw(venn.plot)

# COMMAND ----------

# MAGIC %md
# MAGIC ##NK

# COMMAND ----------

table(liu@meta.data$Condition, liu@meta.data$h_celltype)

# COMMAND ----------

bulk.nk.de <- FindMarkers(object = pseudo_liu, 
                         ident.1 = "NK_lesional", 
                         ident.2 = "NK_healthy",
                         test.use = "DESeq2")
head(bulk.nk.de, n = 15)

# COMMAND ----------

bulk.nk.de.orig <- bulk.nk.de
bulk.nk.de.orig$gene <- rownames(bulk.nk.de.orig)

# COMMAND ----------

# MAGIC %md
# MAGIC Only 1 DEGs
# MAGIC
# MAGIC HSPH1   6.916528e-07  2.7193045  0.75   1.0 0.02345879
# MAGIC
# MAGIC (Same as Reynolds)
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC The rest of the cell types are not interesting or have very low number of cells

# COMMAND ----------

# MAGIC %md
# MAGIC ##Write and save

# COMMAND ----------

write.xlsx(bulk.treg.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_treg_LvsHC_bulk.xlsx")
write.xlsx(bulk.tcell.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_tcell_LvsHC_bulk.xlsx")
write.xlsx(bulk.macro.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_macro_LvsHC_bulk.xlsx")
write.xlsx(bulk.mono.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_mono_LvsHC_bulk.xlsx")
write.xlsx(bulk.dc.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_dc_LvsHC_bulk.xlsx")
write.xlsx(bulk.ilc.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_ilc_LvsHC_bulk.xlsx")
write.xlsx(bulk.nk.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_nk_LvsHC_bulk.xlsx")
write.xlsx(bulk.mastc.de.orig, file="/dbfs/mnt/sandbox/TFM_PAULA/Liu/LvsHC/pseudobulk/liu_mast_LvsHC_bulk.xlsx")
