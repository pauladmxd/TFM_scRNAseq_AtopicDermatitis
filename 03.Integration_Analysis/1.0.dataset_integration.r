# Databricks notebook source
# MAGIC %md
# MAGIC #Task 1: Dataset Integration
# MAGIC
# MAGIC Ensure cell type homogeneity of the current annotations across datasets
# MAGIC
# MAGIC - Integrate scRNA-seq datasets using **RPCA**, **CCA**,**Harmony** to 
# MAGIC check the comparability across datasets.
# MAGIC
# MAGIC - STACAS(Andreatta et al., 2024) is used in other script

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat_v2", .libPaths()))

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)


# COMMAND ----------

alkon <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")
reynolds <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

alkon$dataset <- "alkon"
reynolds$dataset <- "reynolds"

# COMMAND ----------

alkon$celltype_AR <- alkon$h_celltype_v4
reynolds$celltype_AR <- reynolds$h_celltype

# COMMAND ----------

unique(alkon$Condition)

# COMMAND ----------

alkon$Condition <- ifelse(alkon$Condition == "AD", "Lesional", alkon$Condition)

# COMMAND ----------

# # Create a new column 'Condition' with default value 'healthy'
reynolds$Condition <- "HC"

# Update 'Condition' based on the 'Status' and 'Site' columns
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "lesion"] <- "Lesional"
reynolds$Condition[reynolds$Status == "Eczema" & reynolds$Site == "non_lesion"] <- "Non_lesional"


# COMMAND ----------

unique(reynolds$Condition)

# COMMAND ----------

alkon$Condition_AR <- alkon$Condition
reynolds$Condition_AR <- reynolds$Condition

# COMMAND ----------

reynolds$Sample_id <- reynolds$donor_id

# COMMAND ----------

# MAGIC %md
# MAGIC ##Normalization and scale for PCA and UMAP

# COMMAND ----------

AR <- merge(alkon, y = reynolds, add.cell.ids = c("alkon", "reynolds"))

# COMMAND ----------

AR <- NormalizeData(AR)
AR <- FindVariableFeatures(AR)
AR <- ScaleData(AR)
AR <- RunPCA(AR)

# COMMAND ----------

AR <- FindNeighbors(AR, dims = 1:30, reduction = "pca")
AR <- FindClusters(AR, resolution = 2, cluster.name = "unintegrated_clusters")

# COMMAND ----------

# Save the merged Seurat object
 saveRDS(AR, file="/dbfs/mnt/sandbox/TFM_PAULA/MERGED_ARdatasets_TFM.rds")

# COMMAND ----------

AR <- RunUMAP(AR, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# COMMAND ----------

# visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(AR, reduction = "umap.unintegrated", group.by = c("dataset", "celltype_AR"))

# COMMAND ----------

DimPlot(AR, reduction = "pca", group.by = c("dataset", "celltype_AR"))

# COMMAND ----------

AR <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/MERGED_ARdatasets_TFM.rds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##SCTransform for the integration
# MAGIC It is performed in the datasets before the merging, as a normalization on each data sepparately is the correct way to do it.

# COMMAND ----------

# MAGIC %md
# MAGIC Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# MAGIC

# COMMAND ----------

alkon <- SCTransform(alkon, assay = "RNA", new.assay.name = "SCT", conserve.memory= TRUE) #Conserve memory to save memory

# COMMAND ----------

saveRDS(alkon, file="/dbfs/mnt/sandbox/TFM_PAULA/alkon_SCT_TFM.rds")

# COMMAND ----------

reynolds <- SCTransform(reynolds, assay = "RNA", new.assay.name = "SCT", conserve.memory= TRUE) #Conserve memory to save memory

# COMMAND ----------

saveRDS(reynolds, file="/dbfs/mnt/sandbox/TFM_PAULA/reynolds_SCT_TFM.rds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Read the sctransformed objects
# MAGIC

# COMMAND ----------

#Read the sctransform objects
reynolds <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/reynolds_SCT_TFM.rds")
alkon <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/alkon_SCT_TFM.rds")

# COMMAND ----------

AR <- merge(alkon, y = reynolds, add.cell.ids = c("alkon", "reynolds"))

# COMMAND ----------

AR

# COMMAND ----------

DefaultAssay(AR)

# COMMAND ----------

#I need the variable features and they are not stored in the merged object
r_features <- VariableFeatures(reynolds)
a_features <- VariableFeatures(alkon)
common_features <- intersect(a_features, r_features)

# COMMAND ----------

VariableFeatures(AR) <- common_features

# COMMAND ----------

AR

# COMMAND ----------

# MAGIC %md
# MAGIC ###PCA and UMAP before integration to compare later

# COMMAND ----------

AR <- RunPCA(AR)
AR <- FindNeighbors(AR, dims = 1:30, reduction = "pca")
AR <- FindClusters(AR, resolution = 2, cluster.name = "unintegrated_clusters")

# COMMAND ----------

# visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(AR, reduction = "pca", group.by = c("dataset", "celltype_AR"))

# COMMAND ----------

AR <- RunUMAP(AR, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# COMMAND ----------

# visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(AR, reduction = "umap.unintegrated", group.by = c("dataset", "celltype_AR"), label = TRUE)

# COMMAND ----------

AR

# COMMAND ----------

help(IntegrateLayers)

# COMMAND ----------

Reductions(AR)

# COMMAND ----------

DefaultAssay(AR) <- "SCT"

# COMMAND ----------

class(AR)

# COMMAND ----------

AR[["SCT"]]@data

# COMMAND ----------

VariableFeatures(AR)

# COMMAND ----------

# MAGIC %md
# MAGIC ###Integration with SEURAT using ANCHORS
# MAGIC - RPCA
# MAGIC - CCA

# COMMAND ----------

object_list <- list(alkon, reynolds)

# COMMAND ----------

# Step 1: Select integration features
features <- SelectIntegrationFeatures(object.list=object_list, nfeatures=2000)

# COMMAND ----------

alkon <- RunPCA(alkon, features= features)
reynolds <- RunPCA(reynolds, features= features)

# COMMAND ----------

object_list <- list(alkon, reynolds)

# COMMAND ----------

# Step 2: Prepare for SCT integration
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = features)

# COMMAND ----------

# MAGIC %md
# MAGIC ####CCA

# COMMAND ----------

# Step 3: Find integration anchors
anchors_cca <- FindIntegrationAnchors(
  object.list = object_list, 
  normalization.method = "SCT", 
  anchor.features = features, 
  reduction= "cca"
)

# COMMAND ----------

#Step 4: integration
integratedAR_cca <- IntegrateData(
  anchorset=anchors_cca,
  normalization.method="SCT"
)

# COMMAND ----------

integratedAR_cca <- RunPCA(integratedAR_cca, reduction.name = "pca.cca")

# COMMAND ----------

integratedAR_cca <- FindNeighbors(integratedAR_cca, dims = 1:30, reduction = "pca.cca")
integratedAR_cca <- FindClusters(integratedAR_cca, resolution = 2, cluster.name= "clusters.cca")

# COMMAND ----------

#Create a UMAP reduction fo rthe integrated data
integratedAR_cca <- RunUMAP(integratedAR_cca, dims = 1:30, reduction = "pca.cca", reduction.name = "umap.cca")

# COMMAND ----------

#After integration visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)

# COMMAND ----------

DimPlot(integratedAR_cca, reduction = "umap.cca", group.by = c("dataset", "celltype_AR"))

# COMMAND ----------

#Save results
saveRDS(integratedAR_cca, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_CCA_TFM.rds")

#Then read
integratedAR_cca <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_CCA_TFM.rds")

# COMMAND ----------

DimPlot(integratedAR_cca, reduction = "umap.cca", group.by = c("dataset", "celltype_AR"), label=TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC ####RPCA

# COMMAND ----------

# Step 3: Find integration anchors
anchors_rpca <- FindIntegrationAnchors(
  object.list = object_list, 
  normalization.method = "SCT", 
  anchor.features = features, 
  reduction= "rpca"
)

# COMMAND ----------

#Step 4: integration
integratedAR <- IntegrateData(
  anchorset=anchors_rpca,
  normalization.method="SCT"
)

# COMMAND ----------

Reductions(integratedAR)

# COMMAND ----------

integratedAR <- RunPCA(integratedAR, reduction.name = "pca.rpca")

# COMMAND ----------

integratedAR <- FindNeighbors(integratedAR, dims = 1:30, reduction = "pca.rpca")
integratedAR <- FindClusters(integratedAR, resolution = 2, cluster.name= "clusters.rpca")

# COMMAND ----------

#Create a UMAP reduction fo rthe integrated data
integratedAR <- RunUMAP(integratedAR, dims = 1:30, reduction = "pca.rpca", reduction.name = "umap.rpca")

# COMMAND ----------

#After integration visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(integratedAR, reduction = "umap.rpca", group.by = c("dataset", "celltype_AR"))

# COMMAND ----------

# #Save results
# saveRDS(integratedAR, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_RPCA_TFM.rds")
# #Read results
# integratedAR_RPCA <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_RPCA_TFM.rds")
# integratedAR_RPCA <- NULL

# COMMAND ----------

#After integration visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(integratedAR_RPCA, reduction = "umap.rpca", group.by = c("dataset", "celltype_AR"), label=TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Integration with Harmony in Seurat

# COMMAND ----------

#Read the objects scaled with SCTransform
reynolds <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/reynolds_SCT_TFM.rds")
alkon <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/alkon_SCT_TFM.rds")

# COMMAND ----------

alkon$dataset

# COMMAND ----------

AR <- merge(alkon, y = reynolds, add.cell.ids = c("alkon", "reynolds"))

# COMMAND ----------

AR[["SCT"]] <- split(AR[["SCT"]], f = AR$dataset) #Split the object layers counts and scale.counts ito the different layers that correspond to each dataset

# COMMAND ----------

#I need the variable features and they are not stored in the merged object
r_features <- VariableFeatures(reynolds)
a_features <- VariableFeatures(alkon)
common_features <- intersect(a_features, r_features)

# COMMAND ----------

VariableFeatures(AR) <- common_features

# COMMAND ----------

AR <- RunPCA(AR, assay="SCT")

# COMMAND ----------

AR

# COMMAND ----------

#I am having this error:
# Error in IntegrateLayers(object = AR, method = HarmonyIntegration, normalization.method = "SCT") : 
#   None of the features provided are found in this assay
# Error in `IntegrateLayers()`:
# Error in `IntegrateLayers()`:
# ! None of the features provided are found in this assay


# COMMAND ----------

DefaultAssay(AR) <- "SCT"

# COMMAND ----------

FindVariableFeatures(AR)

# COMMAND ----------

# Try to install harmony but problems with rcpp as with presto: Error in unloadNamespace(package) : namespace ‘Rcpp’ is imported by ‘SeuratObject’, ‘reticulate’, ‘uwot’, ‘later’, ‘RcppHNSW’, ‘spam’, ‘ggrepel’, ‘RcppAnnoy’, ‘reshape2’, ‘Rtsne’, ‘promises’, ‘httpuv’, ‘RSpectra’, ‘Seurat’, ‘plyr’ so cannot be unloaded
#  temp_lib <- tempfile()
# dir.create(temp_lib)
# install.packages(c("RhpcBLASctl", "Rcpp", "RcppArmadillo", "RcppProgress"), lib = temp_lib)
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools", lib = temp_lib)
# }
# devtools::install_github("immunogenomics/harmony", lib = temp_lib)
# library(harmony, lib.loc = temp_lib)

# COMMAND ----------

#Using SCT it does not work
# integrated_AR_harmony <- IntegrateLayers(object = AR, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = 'harmony',verbose = TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC As I cannot run the harmony integration using SCT data, I used RNA counts and proceed with the integration the same way, but adding the steps that SCTransform does not required as scale the data

# COMMAND ----------

# Try with assay RNA
AR[["RNA"]] <- split(AR[["RNA"]], f = AR$dataset) #Split the object layers counts and scale.counts ito the different layers that correspond to each dataset
DefaultAssay(AR) <- "RNA"
AR <- FindVariableFeatures(AR)
AR <- ScaleData(AR)
AR <- RunPCA(AR)

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/v13", .libPaths()))
library(harmony) #Required to load the package before runing

# COMMAND ----------

DefaultAssay(AR) <- "RNA"
IntegratedAR_harmony <- IntegrateLayers(object = AR, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony", verbose = TRUE)

# COMMAND ----------

# re-join layers after integration
IntegratedAR_harmony[["RNA"]] <- JoinLayers(IntegratedAR_harmony[["RNA"]])

# COMMAND ----------

IntegratedAR_harmony <- FindNeighbors(IntegratedAR_harmony, reduction = "integrated.harmony", dims = 1:30)
IntegratedAR_harmony <- FindClusters(IntegratedAR_harmony, resolution = 2, cluster.name="clusters.harmony")

# COMMAND ----------

IntegratedAR_harmony <- RunUMAP(IntegratedAR_harmony, dims = 1:30, reduction = "integrated.harmony", reduction.name="umap.harmony")


# COMMAND ----------

IntegratedAR_harmony <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_Harmony_TFM.rds")

# COMMAND ----------

#After integration visualize by batch and cell type annotation
options(repr.plot.width=1600, repr.plot.height=1200)

# COMMAND ----------

DimPlot(IntegratedAR_harmony, reduction = "umap.harmony", group.by = c("dataset", "celltype_AR"), label=TRUE)

# COMMAND ----------

#Save results
saveRDS(IntegratedAR_harmony, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_Harmony_TFM.rds")
