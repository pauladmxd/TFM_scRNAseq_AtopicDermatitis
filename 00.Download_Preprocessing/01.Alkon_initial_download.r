# Databricks notebook source
# MAGIC %md
# MAGIC # Load data from Alkon et al, 2023

# COMMAND ----------

## Append the library folder
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))

# COMMAND ----------

library(Seurat)
library(dplyr)
library(ggplot2)
library(dittoSeq)
#library(UCell)
#library(clusterProfiler)
#library(msigdbr)
library(stringr)
library(forcats)


# COMMAND ----------

AD1.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AD1/")
# Initialize the Seurat object with the raw (non-normalized data).
AD1 <- CreateSeuratObject(counts = AD1.counts)
AD1[["percent.mt"]] <- PercentageFeatureSet(AD1, pattern = "^MT-")
AD1_f <- subset(AD1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AD1_f <- NormalizeData(AD1_f, normalization.method = "LogNormalize", scale.factor = 10000)
AD1_f$Sample_id <- "AD1"
AD1_f$Condition <- "AD"
AD1_f


# COMMAND ----------

AD2.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AD2/")
# Initialize the Seurat object with the raw (non-normalized data).
AD2 <- CreateSeuratObject(counts = AD2.counts)
AD2[["percent.mt"]] <- PercentageFeatureSet(AD2, pattern = "^MT-")
AD2_f <- subset(AD2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AD2_f <- NormalizeData(AD2_f, normalization.method = "LogNormalize", scale.factor = 10000)
AD2_f$Sample_id <- "AD2"
AD2_f$Condition <- "AD"
AD2_f


# COMMAND ----------

AD3.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AD3/")
# Initialize the Seurat object with the raw (non-normalized data).
AD3 <- CreateSeuratObject(counts = AD3.counts)
AD3[["percent.mt"]] <- PercentageFeatureSet(AD3, pattern = "^MT-")
AD3_f <- subset(AD3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AD3_f <- NormalizeData(AD3_f, normalization.method = "LogNormalize", scale.factor = 10000)
AD3_f$Sample_id <- "AD3"
AD3_f$Condition <- "AD"
AD3_f


# COMMAND ----------

AD4.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AD4/")
# Initialize the Seurat object with the raw (non-normalized data).
AD4 <- CreateSeuratObject(counts = AD4.counts)
AD4[["percent.mt"]] <- PercentageFeatureSet(AD4, pattern = "^MT-")
AD4_f <- subset(AD4, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AD4_f <- NormalizeData(AD4_f, normalization.method = "LogNormalize", scale.factor = 10000)
AD4_f$Sample_id <- "AD4"
AD4_f$Condition <- "AD"
AD4_f


# COMMAND ----------

AD5.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AD5/")
# Initialize the Seurat object with the raw (non-normalized data).
AD5 <- CreateSeuratObject(counts = AD5.counts)
AD5[["percent.mt"]] <- PercentageFeatureSet(AD5, pattern = "^MT-")
AD5_f <- subset(AD5, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AD5_f <- NormalizeData(AD5_f, normalization.method = "LogNormalize", scale.factor = 10000)
AD5_f$Sample_id <- "AD5"
AD5_f$Condition <- "AD"
AD5_f


# COMMAND ----------

AP.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/AP/")
# Initialize the Seurat object with the raw (non-normalized data).
AP <- CreateSeuratObject(counts = AP.counts)
AP[["percent.mt"]] <- PercentageFeatureSet(AP, pattern = "^MT-")
AP_f <- subset(AP, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
AP_f <- NormalizeData(AP_f, normalization.method = "LogNormalize", scale.factor = 10000)
AP_f$Sample_id <- "AP"
AP_f$Condition <- "AP"
AP_f

# COMMAND ----------

PN1.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN1/")
# Initialize the Seurat object with the raw (non-normalized data).
PN1 <- CreateSeuratObject(counts = PN1.counts)
PN1[["percent.mt"]] <- PercentageFeatureSet(PN1, pattern = "^MT-")
PN1_f <- subset(PN1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN1_f <- NormalizeData(PN1_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN1_f$Sample_id <- "PN1"
PN1_f$Condition <- "PN"
PN1_f


# COMMAND ----------

PN2.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN2/")
# Initialize the Seurat object with the raw (non-normalized data).
PN2 <- CreateSeuratObject(counts = PN2.counts)
PN2[["percent.mt"]] <- PercentageFeatureSet(PN2, pattern = "^MT-")
PN2_f <- subset(PN2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN2_f <- NormalizeData(PN2_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN2_f$Sample_id <- "PN2"
PN2_f$Condition <- "PN"
PN2_f


# COMMAND ----------

PN3.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN3/")
# Initialize the Seurat object with the raw (non-normalized data).
PN3 <- CreateSeuratObject(counts = PN3.counts)
PN3[["percent.mt"]] <- PercentageFeatureSet(PN3, pattern = "^MT-")
PN3_f <- subset(PN3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN3_f <- NormalizeData(PN3_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN3_f$Sample_id <- "PN3"
PN3_f$Condition <- "PN"
PN3_f


# COMMAND ----------

PN4.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN4/")
# Initialize the Seurat object with the raw (non-normalized data).
PN4 <- CreateSeuratObject(counts = PN4.counts)
PN4[["percent.mt"]] <- PercentageFeatureSet(PN4, pattern = "^MT-")
PN4_f <- subset(PN4, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN4_f <- NormalizeData(PN4_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN4_f$Sample_id <- "PN4"
PN4_f$Condition <- "PN"
PN4_f


# COMMAND ----------

PN5.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN5/")
# Initialize the Seurat object with the raw (non-normalized data).
PN5 <- CreateSeuratObject(counts = PN5.counts)
PN5[["percent.mt"]] <- PercentageFeatureSet(PN5, pattern = "^MT-")
PN5_f <- subset(PN5, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN5_f <- NormalizeData(PN5_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN5_f$Sample_id <- "PN5"
PN5_f$Condition <- "PN"
PN5_f


# COMMAND ----------

PN6.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN6/")
# Initialize the Seurat object with the raw (non-normalized data).
PN6 <- CreateSeuratObject(counts = PN6.counts)
PN6[["percent.mt"]] <- PercentageFeatureSet(PN6, pattern = "^MT-")
PN6_f <- subset(PN6, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN6_f <- NormalizeData(PN6_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN6_f$Sample_id <- "PN6"
PN6_f$Condition <- "PN"
PN6_f


# COMMAND ----------

PN7.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/PN7/")
# Initialize the Seurat object with the raw (non-normalized data).
PN7 <- CreateSeuratObject(counts = PN7.counts)
PN7[["percent.mt"]] <- PercentageFeatureSet(PN7, pattern = "^MT-")
PN7_f <- subset(PN7, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
PN7_f <- NormalizeData(PN7_f, normalization.method = "LogNormalize", scale.factor = 10000)
PN7_f$Sample_id <- "PN7"
PN7_f$Condition <- "PN"
PN7_f


# COMMAND ----------

HC1.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/HC1/")
# Initialize the Seurat object with the raw (non-normalized data).
HC1 <- CreateSeuratObject(counts = HC1.counts)
HC1[["percent.mt"]] <- PercentageFeatureSet(HC1, pattern = "^MT-")
HC1_f <- subset(HC1, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
HC1_f <- NormalizeData(HC1_f, normalization.method = "LogNormalize", scale.factor = 10000)
HC1_f$Sample_id <- "HC1"
HC1_f$Condition <- "HC"
HC1_f


# COMMAND ----------

HC2.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/HC2/")
# Initialize the Seurat object with the raw (non-normalized data).
HC2 <- CreateSeuratObject(counts = HC2.counts)
HC2[["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^MT-")
HC2_f <- subset(HC2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
HC2_f <- NormalizeData(HC2_f, normalization.method = "LogNormalize", scale.factor = 10000)
HC2_f$Sample_id <- "HC2"
HC2_f$Condition <- "HC"
HC2_f


# COMMAND ----------

HC3.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/HC3/")
# Initialize the Seurat object with the raw (non-normalized data).
HC3 <- CreateSeuratObject(counts = HC3.counts)
HC3[["percent.mt"]] <- PercentageFeatureSet(HC3, pattern = "^MT-")
HC3_f <- subset(HC3, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
HC3_f <- NormalizeData(HC3_f, normalization.method = "LogNormalize", scale.factor = 10000)
HC3_f$Sample_id <- "HC3"
HC3_f$Condition <- "HC"
HC3_f


# COMMAND ----------

HC4.counts <- Read10X(data.dir = "/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023//data/HC4/")
# Initialize the Seurat object with the raw (non-normalized data).
HC4 <- CreateSeuratObject(counts = HC4.counts)
HC4[["percent.mt"]] <- PercentageFeatureSet(HC4, pattern = "^MT-")
HC4_f <- subset(HC4, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
HC4_f <- NormalizeData(HC4_f, normalization.method = "LogNormalize", scale.factor = 10000)
HC4_f$Sample_id <- "HC4"
HC4_f$Condition <- "HC"
HC4_f


# COMMAND ----------

# MAGIC %md
# MAGIC # Merge the data and perform integration (RCPA)

# COMMAND ----------

#Remove PN1, very few cells
all_cells <- merge(AD1_f, y = c(AD2_f, AD3_f, AD4_f, AD5_f, AP_f, PN1_f, PN2_f, PN3_f, PN4_f, PN5_f, PN6_f, PN7_f, HC1_f, HC2_f, HC3_f, HC4_f), add.cell.ids = c("AD1","AD2","AD3","AD4","AD5","AP", "PN1", "PN2", "PN3", "PN4", "PN5", "PN6", "PN7", "HC1", "HC2", "HC3", "HC4"))
all_cells

# COMMAND ----------

# split the dataset into a list
cells.list <- SplitObject(all_cells, split.by = "Sample_id")

# COMMAND ----------

# normalize and identify variable features for each dataset independently
cells.list <- lapply(X = cells.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# COMMAND ----------

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = cells.list)
cells.list <- lapply(X = cells.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# COMMAND ----------

anchors <- FindIntegrationAnchors(object.list = cells.list, anchor.features = features, reduction = "rpca")

# COMMAND ----------

obj.combined <- IntegrateData(anchorset = anchors)

# COMMAND ----------

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(obj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.combined <- ScaleData(obj.combined, verbose = FALSE)
obj.combined <- RunPCA(obj.combined, npcs = 50, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:50)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:50)
obj.combined <- FindClusters(obj.combined, resolution = 0.5)

# COMMAND ----------

DefaultAssay(obj.combined) <- "integrated"
obj.combined <- FindClusters(obj.combined, resolution = 0.7)
DefaultAssay(obj.combined) <- "RNA"


# COMMAND ----------

#Count number of cells above certain cutoff (RNA contamination)
assign_labels <- function(gene_counts_matrix) {
  S100A13_pos <- which(grepl("S100A13",rownames(gene_counts_matrix)))
  S100A16_pos <- which(grepl("S100A16",rownames(gene_counts_matrix)))
  S100A14_pos <- which(grepl("S100A14",rownames(gene_counts_matrix)))
  GSTP1_pos <- which(grepl("GSTP1",rownames(gene_counts_matrix)))
  GSN_pos <- which("GSN"==rownames(gene_counts_matrix))
  LGALS7_pos <- which("LGALS7"==rownames(gene_counts_matrix))
  DSP_pos <- which("DSP"==rownames(gene_counts_matrix))
  S100A4_pos <- which(grepl("S100A4",rownames(gene_counts_matrix)))
  S100A10_pos <- which(grepl("S100A10",rownames(gene_counts_matrix)))
  KRT2_pos <- which("KRT2"==rownames(gene_counts_matrix))
  LY6D_pos <- which(grepl("LY6D",rownames(gene_counts_matrix)))
  KRT6B_pos <- which(grepl("KRT6B",rownames(gene_counts_matrix)))
  KRT6C_pos <- which(grepl("KRT6C",rownames(gene_counts_matrix)))
  FABP5_pos <- which(grepl("FABP5",rownames(gene_counts_matrix)))
  KRT17_pos <- which(grepl("KRT17",rownames(gene_counts_matrix)))
  DMKN_pos <- which(grepl("DMKN",rownames(gene_counts_matrix)))
  AQP3_pos <- which(grepl("AQP3",rownames(gene_counts_matrix)))
  CXCL14_pos <- which(grepl("CXCL14",rownames(gene_counts_matrix)))
  PERP_pos <- which(grepl("PERP",rownames(gene_counts_matrix)))
  S100A2_pos <- which(grepl("S100A2",rownames(gene_counts_matrix)))
  SFN_pos <- which(grepl("SFN",rownames(gene_counts_matrix)))
  KRT1_pos <- which("KRT1"==rownames(gene_counts_matrix))
  KRT10_pos <- which(grepl("KRT10",rownames(gene_counts_matrix)))
  KRT5_pos <- which(grepl("KRT5",rownames(gene_counts_matrix)))
  KRT14_pos <- which(grepl("KRT14",rownames(gene_counts_matrix)))
  S100A6_pos <- which(grepl("S100A6",rownames(gene_counts_matrix)))
  S100A7_pos <- which("S100A7"==rownames(gene_counts_matrix))
  S100A8_pos <- which(grepl("S100A8",rownames(gene_counts_matrix)))
  S100A9_pos <- which(grepl("S100A9",rownames(gene_counts_matrix)))
  KRT6A_pos <- which(grepl("KRT6A",rownames(gene_counts_matrix)))
  KRT16_pos <- which(grepl("KRT16",rownames(gene_counts_matrix)))
  KRTDAP_pos <- which(grepl("KRTDAP",rownames(gene_counts_matrix)))
  LGALS7B_pos <- which(grepl("LGALS7B",rownames(gene_counts_matrix)))

  labels <- apply(gene_counts_matrix, 2, function(gene_expression) {
    if (gene_expression[S100A13_pos] > 2 &
    gene_expression[S100A16_pos] > 2 &    
    gene_expression[S100A14_pos] > 2.5 &
    gene_expression[GSTP1_pos] > 2.5 &
    gene_expression[GSN_pos] > 2.5 &
    gene_expression[LGALS7_pos] > 3 &
    gene_expression[DSP_pos] > 3 &
    gene_expression[S100A4_pos] > 3 &
    gene_expression[S100A10_pos] > 3 &
    gene_expression[KRT2_pos] > 3 &
    gene_expression[LY6D_pos] > 3 &
    gene_expression[KRT6B_pos] > 3 &
    gene_expression[KRT6C_pos] > 3 &
    gene_expression[FABP5_pos] > 3 &
    gene_expression[KRT17_pos] > 3 &
    gene_expression[DMKN_pos] > 3 &
    gene_expression[AQP3_pos] > 3 &
    gene_expression[CXCL14_pos] > 3 &
    gene_expression[PERP_pos] > 3 &
    gene_expression[S100A2_pos] > 3.5 &
    gene_expression[SFN_pos] > 3.5 &
    gene_expression[KRT1_pos] > 4 &
    gene_expression[KRT10_pos] > 4 &
    gene_expression[KRT5_pos] > 4 &
    gene_expression[KRT14_pos] > 4 &
    gene_expression[S100A6_pos] > 4 &
    gene_expression[S100A7_pos] > 4 &
    gene_expression[S100A8_pos] > 4 &
    gene_expression[S100A9_pos] > 4 &
    gene_expression[KRT6A_pos] > 4 &    
    gene_expression[KRT16_pos] > 4 &   
    gene_expression[KRTDAP_pos] > 4 &   
    gene_expression[LGALS7B_pos] > 4) {
      return("RNA contaminated")
    } else {
      return("Clean")
    }
  })
  return(labels)
}


# COMMAND ----------

labels <- assign_labels(obj.combined@assays$RNA@data)

# COMMAND ----------

table(labels)

# COMMAND ----------

DimPlot(obj.combined, group.by="Sample_id", raster=F) +
NoAxes()

# COMMAND ----------

DimPlot(obj.combined, group.by="Condition", raster=F) +
NoAxes()

# COMMAND ----------

DimPlot(obj.combined, raster=F, label=T) +
NoAxes()

# COMMAND ----------

# MAGIC %md
# MAGIC # Get the markers between the unsup clusters
# MAGIC

# COMMAND ----------

# Get the clusters
DefaultAssay(obj.combined) <- "RNA"
all.markers <- FindAllMarkers(obj.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)


# COMMAND ----------

#Save the assembled object
saveRDS(obj.combined, file="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_integrated_obj.rds")

# COMMAND ----------

# Save the markers
write.table(all.markers,file="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_DEG_unsup_clusters.0.5.tsv",sep="\t",quote=F,row.names = TRUE,
            col.names = TRUE)

# COMMAND ----------

obj.combined <- readRDS( file="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_integrated_obj.rds")

# COMMAND ----------

# Save the markers
all.markers <- read.table(file="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_DEG_unsup_clusters.0.5.tsv",sep="\t",header = TRUE)
table(all.markers$cluster)

# COMMAND ----------

options(repr.plot.width=500, repr.plot.height=300)


# COMMAND ----------

dittoBarPlot(obj.combined, "Condition", group.by="Condition",   scale = c( "count"))

# COMMAND ----------

table(obj.combined$Condition)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000)


# COMMAND ----------

DefaultAssay(obj.combined)

# COMMAND ----------

genes <- c("CD3D","CD4","CD8A","MKI67","IL7R","FOXP3","IFNG","IL13","CRTAM","IL9","KLRD1","KLRB1","CD79A","CD207","ITGAX","LAMP3","CLEC9A","CD163","TPSAB1","NRXN1","SCN7A","MLANA","ACTA2","LYVE1","PECAM1","COL1A1","KRT18","KRT19","KRT14","KRT10","PTPRC")
dittoDotPlot(obj.combined,vars=genes,group.by="seurat_clusters") +
    RotatedAxis()

# COMMAND ----------

# T cell: 4,18,25,30
#   CD4: 
#   CD8/NK: 18
#   Treg: 4,25,30
# Prolif: 6,25,27,32
# B: 24
# LC: 10
# DC: 21,28
# MacroPhages: 16
# Mast: 31
# Melanocytes: 19
# Peri: 23
# Smooth muscle: 12,14
# lymphatic endothelial: 15
# blood endothelial: 2,8,9,30
# fb: 3,5,11,12,20
# sweat gland: 26
# keratinocytes: 0,1,6,7,13,17,29,22

# COMMAND ----------

new.cluster.ids <- c("0 Kt","1 Kt","2 Blood end","3 Fb","4 CD4/Treg","5 Fb","6 Prolif","7 Kt","8 Blood end","9 Blood end","10 LC","11 Fb","12 Smooth","13 Kt","14 Smooth","15 Lymph end","16 Macro","17 Kt","18 CD8/NK","19 Mel","20 Fb","21 DC","22 Kt","23 Peric","24 B","25 Prolif","26 Sweat","27 Prolif","28 DC","29 Kt/Fb","30 Blood end","31 Mast","32 Prolif")
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
obj.combined$cell_population <- Idents(obj.combined)

# COMMAND ----------

DimPlot(obj.combined, raster=F, label=T) +
NoAxes()

# COMMAND ----------

new.cluster.ids <- c("Kt","Kt","Blood end","Fb","T cell","Fb","Kt Prolif","Kt","Blood end","Blood end","Myeloid","Fb","Smooth","Kt","Smooth","Lymph end","Myeloid","Kt","T cell","Mel","Fb","Myeloid","Kt","Peric","B","T cell Prolif","Sweat","Myeloid Prolif","Myeloid","Kt/Fb","Blood end","Mast","Kt Prolif")
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
obj.combined$cell_population_v2 <- Idents(obj.combined)
#Idents(obj.combined) <- obj.combined$cell_population

# COMMAND ----------

DimPlot(obj.combined, raster=F, label=T) +
NoAxes()

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=800)


# COMMAND ----------

FeaturePlot(obj.combined,features="FOXP3",raster=F)

# COMMAND ----------

FeaturePlot(obj.combined,features="IKZF2",raster=F)

# COMMAND ----------

#Save the assembled object
saveRDS(obj.combined, file="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_integrated_obj_annotated.rds")
