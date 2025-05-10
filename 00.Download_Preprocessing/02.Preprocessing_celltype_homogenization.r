# Databricks notebook source
# MAGIC %md
# MAGIC #DATASET PREPARATION 
# MAGIC
# MAGIC  Here I will change the annotations for the desired scRNA-seq datasets to have a homogeneus annotation and be able to compare all dataset information together. And also to have the relevant information for the study.
# MAGIC
# MAGIC

# COMMAND ----------

# MAGIC %md
# MAGIC Fisrt, it is important to analyze which are the current annotations in the datasets to see how is it possible to collapse them and make them homogeneus
# MAGIC
# MAGIC Then, I will remove samples of conditions that are not interesting for my study.

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))

library(Seurat)
library(dittoSeq)
library(dplyr)
library(openxlsx)

# COMMAND ----------

# MAGIC %md
# MAGIC ## AD/PS Reynolds et al, 2021
# MAGIC - Includes lesional and non-lesional samples

# COMMAND ----------

reynolds <- readRDS(file="/dbfs/mnt/sandbox/Reynolds/misc/4_obj_processed_allcells.rds")

# COMMAND ----------

table(reynolds$donor_id, reynolds$Status)

# COMMAND ----------

head(reynolds@meta.data, 5)

# COMMAND ----------

unique(reynolds$Status)

# COMMAND ----------

reynolds$h_celltype <- gsub("Undifferentiated_KC.*|Proliferating_KC|Differentiated_KC.*|Differentiated_KC", "KC", reynolds$final_clustering)
reynolds$h_celltype <- gsub("Melanocyte", "Melanocytes", reynolds$h_celltype)
reynolds$h_celltype <- gsub("LC_4|LC_3|LC_1|LC_2", "LC", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Th|Tc", "TC", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Treg", "Treg", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Mono|Inf_mono", "Mono", reynolds$h_celltype)
reynolds$h_celltype <- gsub("MigDC|moDC_1|moDC_2|moDC_3|DC1|DC2", "DC", reynolds$h_celltype)
reynolds$h_celltype <- gsub("ILC1_3|ILC1_NK|ILC2", "ILC", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Macro_1|Macro_2", "Macro", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Mast_cell", "MastC", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Pericyte_1_non_inflamm|Pericyte_2_inflamm", "Pericyte", reynolds$h_celltype)
reynolds$h_celltype <- gsub("Schwann1|Schwann2", "Schwann", reynolds$h_celltype)
reynolds$h_celltype <- gsub("F1|F2|F3", "Fibroblasts", reynolds$h_celltype)
reynolds$h_celltype <- gsub("VE1|VE2|VE3", "VE", reynolds$h_celltype)
reynolds$h_celltype <- gsub("LE1|LE2", "LE", reynolds$h_celltype)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200)
DimPlot(reynolds, group.by = "h_celltype", label = TRUE)

# COMMAND ----------

unique(reynolds$h_celltype)

# COMMAND ----------

# MAGIC %md
# MAGIC Here we have annotated KC as keratinocytes, Melanocytes,  LC, Fibroblasts, VE (vascular endothelium), ILC (innate) , LE (lymphoid endothelium) , TC (t cells) , Treg (regulatory T cells), Mono (monocytes) ,DC (dendritic cells) , NK,  Pericyte, Macro (macrophagues), Schwann, MastC,   Plasma   

# COMMAND ----------

unique(reynolds$Status) 
#In Status we can find to wich condition it belongs, Eczema is AD.
desiredR<- c("Eczema", "Healthy") #These are the relevant conditions I need

# COMMAND ----------

#Subset of the dataset with only AD and HC samples
subreynolds <- subset(reynolds,  subset = Status %in% desiredR)

# COMMAND ----------

#Save the dataset for next uses
#saveRDS(subreynolds, file="/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##AD/PN – Alkon et al, 2023
# MAGIC - Does not include non-lesional samples.
# MAGIC - 5 AD, 7 PN, 1 AP (atopic prurigo) and 4 healthy control​​

# COMMAND ----------

alkon <- readRDS(file ="/dbfs/mnt/sandbox/AD_PN/Alkon_et_al_2023/misc/1_integrated_obj_annotated.rds")

# COMMAND ----------

DefaultAssay(alkon) #if it is RNA we have the original counts (important to work with them as there was an integration performed)

# COMMAND ----------

head(alkon@meta.data, 5)

# COMMAND ----------

alkon$h_celltype <-alkon$cell_population
alkon$h_celltype <- gsub("21 DC", "Plasmacytoid DC", alkon$h_celltype) # These DC are clustered very separated from myeloid cells, according to literature it may be plasmocytoid DC
alkon$h_celltype <- gsub("^\\d+\\s*", "", alkon$h_celltype) # to remove the numbers followed by space

# COMMAND ----------

alkon$h_celltype <- gsub("Kt", "KC", alkon$h_celltype)
alkon$h_celltype <- gsub(".*Fb.*", "Fibroblasts", alkon$h_celltype)
alkon$h_celltype <- gsub("Blood end", "Endothelial blood", alkon$h_celltype)
alkon$h_celltype <- gsub("Smooth", "Smooth Muscle", alkon$h_celltype)
alkon$h_celltype <- gsub("Lymph end", "LE", alkon$h_celltype)
alkon$h_celltype <- gsub("Mel", "Melanocytes", alkon$h_celltype)
alkon$h_celltype <- gsub("Peric", "Pericyte", alkon$h_celltype)
alkon$h_celltype <- gsub("Sweat", "Sweat Gland", alkon$h_celltype)
alkon$h_celltype <- gsub("B", "B cells", alkon$h_celltype)
alkon$h_celltype <- gsub("Mast", "MastC", alkon$h_celltype)

# COMMAND ----------

Idents(alkon) <- alkon$h_celltype

# COMMAND ----------

DimPlot(alkon, group.by="cell_population_v2", label = TRUE)

# COMMAND ----------

DimPlot(alkon, group.by="h_celltype", label = TRUE)

# COMMAND ----------

unique(alkon$h_celltype)

# COMMAND ----------

# MAGIC %md
# MAGIC ### Macro subcluster

# COMMAND ----------

# MAGIC %md
# MAGIC I want to do a subcluster in Macro cells to separate them into Macrophagues and monocytes, to make it comparable to the other 2 datasets.

# COMMAND ----------

subset_macro <- subset(alkon, idents = "Macro") #First a subset with only macro cells

# COMMAND ----------

# MAGIC %md
# MAGIC Now it has to be performed the **normalization and dimentionality reduction** again in the subset to be able to subcluster this cell type.

# COMMAND ----------

subset_macro <- NormalizeData(subset_macro, normalization.method = "LogNormalize", scale.factor = 10000) #default values


# COMMAND ----------

Idents(subset_macro) <- subset_macro$h_celltype

# COMMAND ----------

subset_macro <- FindVariableFeatures(subset_macro, selection.method = "vst", nfeatures = 2000)

# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(subset_macro), 10)

# # plot variable features with labels
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2

# COMMAND ----------

#Scale data
all.genes <- rownames(subset_macro)
subset_macro <- ScaleData(subset_macro, features = all.genes) # si fuera un objeto muy grande tarda menos poner los high variable features

# COMMAND ----------

#Linear dimensionality reduction
subset_macro <- RunPCA(subset_macro, features = VariableFeatures(object = subset_macro))


# COMMAND ----------

ElbowPlot(subset_macro) #I just want 2 subclusters for macrophages and monocytes


# COMMAND ----------

# MAGIC %md
# MAGIC Find new clusters

# COMMAND ----------

subset_macro <- FindNeighbors(subset_macro, dims = 1:20)
subset_macro <- FindClusters(subset_macro, resolution = 0.03) #Use the seurat function to find new clusters

# COMMAND ----------

# MAGIC %md
# MAGIC Recompute UMAP with new clusters to see if they are correctly separated

# COMMAND ----------

#Recompute UMAP
subset_macro <- RunUMAP(subset_macro, dims = 1:20)

# COMMAND ----------

DimPlot(subset_macro, reduction = "umap", group.by = "seurat_clusters")

# COMMAND ----------

# MAGIC %md
# MAGIC Chenking in PanglaoDB I look for the common markers of Macro and Mono
# MAGIC
# MAGIC **Macro**: CD68, NAA, JAML, TYROBP --> _TOP in panglao by votes (not high UI markers for this cell type)_,
# MAGIC
# MAGIC **Mono**: RHOC, IFITM3, ZFP36L2,--> _markers with higher UI (ubiquitinous index, the specificity of these marker only in this cluster),_ 
# MAGIC
# MAGIC ** As it was so difficult to differentiate I also search markers in:
# MAGIC - https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2014.00514/full
# MAGIC - https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0279-4 
# MAGIC - https://www.oaepublish.com/articles/2574-1209.2019.04

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200)

DotPlot(subset_macro, features = c("CD68", "NAA", "JAML","TYROBP","RHOC", "IFITM3", "ZFP36L2"), group.by = "seurat_clusters")

# COMMAND ----------

# MAGIC %md
# MAGIC Now it is necessary to update the annotation with the new subclusters

# COMMAND ----------

#New annotation 
new.cluster.ids.macro <- c("Macro", "Mono")
names(new.cluster.ids.macro) <- levels(subset_macro)
subset_macro <- RenameIdents(subset_macro, new.cluster.ids.macro)

# COMMAND ----------

table(Idents(subset_macro))

# COMMAND ----------

subset_macro$h_celltype <- Idents(subset_macro) # Assign the new annotation to the celltypes in the subset
unique(subset_macro$h_celltype)

# COMMAND ----------

last_anno_macro <- as.character(alkon$h_celltype) # we have to save it as a character because it was a factor, and if not, it would not take the names correctly, it would save the factor name (1 and 2)
table(last_anno_macro)

# COMMAND ----------

subset_macro$h_celltype <- Idents(subset_macro)

# COMMAND ----------

last_anno_macro[which(colnames(alkon)%in%colnames(subset_macro))] <- as.character(subset_macro$h_celltype) #here the new annotation is assigned to the cells that were in the last annotation and now are in the subset
table(last_anno_macro)

# COMMAND ----------

#New annotations does not have the names, so I add them manually
alkon$h_celltype <- gsub("1", "Macro", alkon$h_celltype)
alkon$h_celltype <- gsub("2", "Mono", alkon$h_celltype)

# COMMAND ----------

#Now finally the cell types are updated in the original object
names(last_anno_macro) <- colnames(alkon)
alkon <- AddMetaData(alkon, last_anno_macro, col.name="h_celltype_v2")
table(alkon$h_celltype_v2) #to check

# COMMAND ----------

unique(alkon$h_celltype) #to check

# COMMAND ----------

options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(alkon, group.by="h_celltype", label = TRUE) + DimPlot(alkon, group.by="h_celltype_v2", label = TRUE)

# COMMAND ----------

table(alkon@meta.data$Condition, alkon@meta.data$h_celltype_v2)

# COMMAND ----------

# MAGIC %md
# MAGIC ### NK/CD8+ subcluster
# MAGIC Same procedure as before

# COMMAND ----------

subset_nk <- subset(alkon, idents = "CD8/NK") #First a subset with only NK/CD8+ cells

# COMMAND ----------

subset_nk <- NormalizeData(subset_nk, normalization.method = "LogNormalize", scale.factor = 10000) #default values


# COMMAND ----------

subset_nk <- FindVariableFeatures(subset_nk, selection.method = "vst", nfeatures = 2000)

# # Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(subset_nk), 10)
top10
# # plot variable features with labels
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2

# COMMAND ----------

#Scale data
all.genes <- rownames(subset_nk)
subset_nk <- ScaleData(subset_nk, features = all.genes)

# COMMAND ----------

subset_nk <- RunPCA(subset_nk, features = VariableFeatures(object = subset_nk))

# COMMAND ----------

ElbowPlot(subset_nk)

# COMMAND ----------

subset_nk <- FindNeighbors(subset_nk, dims = 1:20)
subset_nk <- FindClusters(subset_nk, resolution = 0.05) #Use the seurat function to find new clusters

# COMMAND ----------

# MAGIC %md
# MAGIC Recompute UMAP

# COMMAND ----------

#Recompute UMAP
subset_nk <- RunUMAP(subset_nk, dims = 1:20)

# COMMAND ----------

DimPlot(subset_nk, reduction = "umap", group.by = "seurat_clusters")

# COMMAND ----------

# MAGIC %md
# MAGIC In PanglaoDB I look for the common markers of NK and CD8 cells
# MAGIC
# MAGIC NK: TRDC, NKG7, KLRD1, KLRF1 
# MAGIC
# MAGIC T-cytotoxic (CD8): CD8A, TRAC, GZMB, PRF1

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200)

DotPlot(subset_nk, features = c("KLRD1", "KLRF1", "TRDC", "NKG7", "CD8A", "TRAC", "GZMB", "PRF1"), group.by = "seurat_clusters")

# COMMAND ----------

#New annotation 
new.cluster.ids.nk <- c("TC", "NK")
names(new.cluster.ids.nk) <- levels(subset_nk)
subset_nk <- RenameIdents(subset_nk, new.cluster.ids.nk)

# COMMAND ----------

table(Idents(subset_nk))

# COMMAND ----------

subset_nk$h_celltype <- Idents(subset_nk) # Assign the new annotation to the celltypes in the subset
unique(subset_nk$h_celltype)

# COMMAND ----------

subset_nk$h_celltype <- Idents(subset_nk)
# Update the main dataset with new annotations
alkon$h_celltype[names(subset_nk$h_celltype)] <- subset_nk$h_celltype

# COMMAND ----------

last_anno_nk <- as.character(alkon$h_celltype_v2)
table(last_anno_nk)

# COMMAND ----------

last_anno_nk[which(colnames(alkon)%in%colnames(subset_nk))] <- as.character(subset_nk$h_celltype)
table(last_anno_nk)

# COMMAND ----------

names(last_anno_nk) <- colnames(alkon)
alkon <- AddMetaData(alkon, last_anno_nk, col.name="h_celltype_v3")
table(alkon$h_celltype_v3)

# COMMAND ----------

options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(alkon, group.by="h_celltype_v2", label=T) + DimPlot(alkon, group.by="h_celltype_v3", label=T)

# COMMAND ----------

DimPlot(alkon, group.by="h_celltype_v2", label=T)

# COMMAND ----------

#New annotations does not have the names, so I add them manually
alkon$h_celltype <- gsub("1", "Tcit", alkon$h_celltype)
alkon$h_celltype <- gsub("2", "NK", alkon$h_celltype)

# COMMAND ----------

DimPlot(alkon, group.by="h_celltype", label = TRUE)

# COMMAND ----------

# MAGIC %md
# MAGIC ###CD4+/Treg subcluster

# COMMAND ----------

subset_cd4 <- subset(alkon, idents = "CD4/Treg") #First a subset with only CD4+/Treg cells

# COMMAND ----------

subset_cd4 <- NormalizeData(subset_cd4, normalization.method = "LogNormalize", scale.factor = 10000) #default values


# COMMAND ----------

subset_cd4 <- FindVariableFeatures(subset_cd4, selection.method = "vst", nfeatures = 2000)

# # Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(subset_cd4), 10)
top10
# # plot variable features with labels
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2

# COMMAND ----------

#Scale data
all.genes <- rownames(subset_cd4)
subset_cd4 <- ScaleData(subset_cd4, features = all.genes)

# COMMAND ----------

subset_cd4 <- RunPCA(subset_cd4, features = VariableFeatures(object = subset_cd4))

# COMMAND ----------

ElbowPlot(subset_cd4)

# COMMAND ----------

subset_cd4 <- FindNeighbors(subset_cd4, dims = 1:20)
subset_cd4 <- FindClusters(subset_cd4, resolution = 0.03) #Use the seurat function to find new clusters

# COMMAND ----------

# MAGIC %md
# MAGIC Recompute UMAP

# COMMAND ----------

#Recompute UMAP
subset_cd4 <- RunUMAP(subset_cd4, dims = 1:20)

# COMMAND ----------

DimPlot(subset_cd4, reduction = "umap", group.by = "seurat_clusters")

# COMMAND ----------

# MAGIC %md
# MAGIC In PanglaoDB I look for the common markers of Treg and Th cells
# MAGIC
# MAGIC Treg: FOXP3, IL2RA, CTLA4, IKZF2 
# MAGIC
# MAGIC Th: IL7R, CD3G, IL13, CD4, CD28

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200)

DotPlot(subset_cd4, features = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "IL7R", "CD3G", "IL13", "CD4", "CD28"), group.by = "seurat_clusters")

# COMMAND ----------

#New annotation 
new.cluster.ids.treg <- c("Treg", "TC")
names(new.cluster.ids.treg) <- levels(subset_cd4)
subset_cd4 <- RenameIdents(subset_cd4, new.cluster.ids.treg)

# COMMAND ----------

table(Idents(subset_cd4))

# COMMAND ----------

subset_cd4$h_celltype <- Idents(subset_cd4)
# Update the main dataset with new annotations
alkon$h_celltype[names(subset_cd4$h_celltype)] <- subset_cd4$h_celltype

# COMMAND ----------

subset_cd4$h_celltype <- Idents(subset_cd4) # Assign the new annotation to the celltypes in the subset
unique(subset_cd4$h_celltype)

# COMMAND ----------

#New annotations does not have the names, so I add them manually
alkon$h_celltype <- gsub("1", "Treg", alkon$h_celltype)
alkon$h_celltype <- gsub("2", "Th", alkon$h_celltype)

#Also group Th and Tcit in Tcells like in the other datasets
alkon$h_celltype <- gsub("Th|Tcit", "TC", alkon$h_celltype)

# COMMAND ----------

last_anno_treg <- as.character(alkon$h_celltype_v3)
table(last_anno_treg)

# COMMAND ----------

last_anno_treg[which(colnames(alkon)%in%colnames(subset_cd4))] <- as.character(subset_cd4$h_celltype)
table(last_anno_treg)

# COMMAND ----------

names(last_anno_treg) <- colnames(alkon)
alkon <- AddMetaData(alkon, last_anno_treg, col.name="h_celltype_v4")
table(alkon$h_celltype_v4)

# COMMAND ----------

options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(alkon, group.by="h_celltype_v3", label=T) + DimPlot(alkon, group.by="h_celltype_v4", label=T)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200)
# Define unique colors for each cell type
colors <- c("Macro" = "#1f77b4", "KC" = "#ff7f0e", "Prolif" = "#2ca02c", 
            "Treg" = "#aec7e8", "TC" = "#9467bd", "DC" = "#8c564b", 
            "Fibroblasts" = "#e377c2", "Endothelial blood" = "#7f7f7f", 
            "Smooth Muscle" = "#bcbd22", "LE" = "#17becf", "LC" = "#d62728", 
            "Melanocytes" = "#ffbb78", "Plasmacytoid DC" = "#98df8a", 
            "Pericyte" = "#ff9896", "Mono" = "#c5b0d5", "Sweat Gland" = "#c49c94", 
            "B cells" = "#f7b6d2", "MastC" = "#c7c7c7", "NK" = "#dbdb8d")

# Create DimPlot with the final annotation with specified colors to distinguish better
DimPlot(alkon, group.by = "h_celltype_v4", cols=colors, label =T)

# COMMAND ----------

# MAGIC %md
# MAGIC Here we have annotated Macro, Mono, LC, DC, KC, TC, Treg, NK, Fibroblasts,Endothelial blood, Smooth Muscle, LE, Melanocytes, Pericyte, KC/Fibroblasts, Sweat Gland, B cells, MastC 

# COMMAND ----------

table(alkon@meta.data$Condition, alkon@meta.data$h_celltype_v4)

# COMMAND ----------

# MAGIC %md
# MAGIC Now it is only saved the relevant conditions for my study, AD and HC.

# COMMAND ----------

unique(alkon$Condition)
desired <- c("AD", "HC")

# COMMAND ----------

subalkon <- subset(alkon, subset= Condition %in% desired)

# COMMAND ----------

#saveRDS(subalkon, file="/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")

# COMMAND ----------

# MAGIC %md
# MAGIC #Dataset Integration
# MAGIC I need to integrate the 2 datasets to be able to compare what is on each one.

# COMMAND ----------

alkon <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")
reynolds <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Check if the cell annotation is similar by looking the markers of each cell among the datasets

# COMMAND ----------

subalkon_cells <- c("KC", "TC", "Fibroblasts", "Melanocytes", "NK", "ILC", "DC", "Macro", "Mono", "Treg") #I did not write correctly fibroblasts so are not included
subalkon <- subset(alkon, h_celltype_v4 %in% subalkon_cells) #relevant cells for the study

# COMMAND ----------

unique(subalkon$h_celltype_v4)

# COMMAND ----------

Idents(subalkon) <- subalkon$h_celltype_v4
alkon_celltype_markers <- FindAllMarkers(subalkon)

# COMMAND ----------

library(openxlsx)

# Save alkon markers to Excel
write.xlsx(alkon_celltype_markers, "/dbfs/mnt/sandbox/TFM_PAULA/alkon_celltype_markers.xlsx")

# COMMAND ----------

#As I didnt save the fibroblast markers due to a writting error, I have to redo it
Idents(subalkon) <- subalkon$h_celltype_v4
alkon_fibroblast_markers <- FindMarkers(subalkon, ident.1 = "Fibroblasts", ident.2 = NULL)

# COMMAND ----------

alkon_fibroblast_markers$gene <- rownames(alkon_fibroblast_markers)

# COMMAND ----------

#Save fibroblasts markers
write.xlsx(alkon_fibroblast_markers, "/dbfs/mnt/sandbox/TFM_PAULA/alkon_fb_markers1.xlsx")

# COMMAND ----------

subreynolds <- subset(reynolds, h_celltype %in% subalkon_cells) #relevant cells for the study

# COMMAND ----------

Idents(subreynolds) <- subreynolds$h_celltype
reynolds_celltype_markers <- FindAllMarkers(subreynolds)

# COMMAND ----------

# Save reynolds markers to Excel
library(openxlsx)
write.xlsx(reynolds_celltype_markers, "/dbfs/mnt/sandbox/TFM_PAULA/reynolds_celltype_markers.xlsx")

# COMMAND ----------

#As I didnt save the fibroblast markers due to a writting error, I have to redo it
Idents(subreynolds) <- subreynolds$h_celltype
reynolds_fibroblast_markers <- FindMarkers(subreynolds, ident.1 = "Fibroblasts", ident.2 = NULL)

# COMMAND ----------

reynolds_fibroblast_markers$gene <- rownames(reynolds_fibroblast_markers)

# COMMAND ----------

# Save fibroblast markers
write.xlsx(reynolds_fibroblast_markers, "/dbfs/mnt/sandbox/TFM_PAULA/reynolds_fb_markers1.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Venn diagram of common markers 
# MAGIC Now read the markers and find the common ones across celltypes.

# COMMAND ----------

alkon_celltype_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/alkon_celltype_markers.xlsx")
reynolds_celltype_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/reynolds_celltype_markers.xlsx")

# COMMAND ----------

# Filter significant ones + absolute value minimum of 2
alkon_celltype_markers <- alkon_celltype_markers[alkon_celltype_markers$p_val_adj < 0.05 & abs(alkon_celltype_markers$avg_log2FC) >= 2, ]
reynolds_celltype_markers <- reynolds_celltype_markers[reynolds_celltype_markers$p_val_adj < 0.05 & abs(reynolds_celltype_markers$avg_log2FC) >= 2, ]

# COMMAND ----------

unique(reynolds_celltype_markers$cluster)

# COMMAND ----------

.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))

library(VennDiagram)
library(RColorBrewer)
library(grid)

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "TC"
alkon_tc_markers <- alkon_celltype_markers[alkon_celltype_markers$cluster == "TC", ]
reynolds_tc_markers <- reynolds_celltype_markers[alkon_celltype_markers$cluster == "TC", ]

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_tc_markers$gene,
    reynolds = reynolds_tc_markers$gene
  ),
  category.names = c( "Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log=FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-3, 3, 180),
  cat.dist = c(0.03, 0.03, 0.03),  # Adjusted distance for Reynolds Dataset
  cat.fontfamily = "sans",
  # Title
  main = "Tcell markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)
# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "Treg"
alkon_treg_markers <- alkon_celltype_markers[alkon_celltype_markers$cluster == "Treg", ]
reynolds_treg_markers <- reynolds_celltype_markers[reynolds_celltype_markers$cluster == "Treg", ]

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_treg_markers$gene,
    reynolds = reynolds_treg_markers$gene
  ),
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log=FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-3, 3, 180),
  cat.dist = c(0.03, 0.03, 0.03),  # Adjusted distance for Reynolds Dataset
  cat.fontfamily = "sans",
  # Title
  main = "Treg markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)
# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "NK"
alkon_nk_markers <- alkon_celltype_markers[alkon_celltype_markers$cluster == "NK", ]
reynolds_nk_markers <- reynolds_celltype_markers[reynolds_celltype_markers$cluster == "NK", ]

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_nk_markers$gene,
    reynolds = reynolds_nk_markers$gene
  ),
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log=FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-3, 3, 180),
  cat.dist = c(0.03, 0.03, 0.03),  # Adjusted distance for Reynolds Dataset
  cat.fontfamily = "sans",
  # Title
  main = "NK markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)
# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "KC"
alkon_kc_markers <- alkon_celltype_markers[alkon_celltype_markers$cluster == "KC", ]
reynolds_kc_markers <- reynolds_celltype_markers[reynolds_celltype_markers$cluster == "KC", ]

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_kc_markers$gene,
    reynolds = reynolds_kc_markers$gene
  ),
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log = FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-3, 3),
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans",
  
  # Title
  main = "Keratinocytes markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)

# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "Melanocytes"
alkon_melanocytes_markers <- alkon_celltype_markers[alkon_celltype_markers$cluster == "Melanocytes", ]
reynolds_melanocytes_markers <- reynolds_celltype_markers[reynolds_celltype_markers$cluster == "Melanocytes", ]

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_melanocytes_markers$gene,
    reynolds = reynolds_melanocytes_markers$gene
  ),
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log = FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-3, 3),
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans",
  
  # Title
  main = "Melanocytes markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)

# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

alkon_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/alkon_fb_markers1.xlsx")
reynolds_fb_markers <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/reynolds_fb_markers1.xlsx")

# COMMAND ----------

alkon_fb_markers <- alkon_fb_markers[alkon_fb_markers$p_val_adj < 0.05 & abs(alkon_fb_markers$avg_log2FC) >= 2, ]
reynolds_fb_markers <- reynolds_fb_markers[reynolds_fb_markers$p_val_adj < 0.05 & abs(reynolds_fb_markers$avg_log2FC) >= 2, ]

# COMMAND ----------

# Prepare a palette of 3 colors with R colorbrewer
myCol <- brewer.pal(3, "Pastel2")

# Create a Venn diagram for celltype "fibroblasts"

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = list(
    alkon = alkon_fb_markers$gene,
    reynolds = reynolds_fb_markers$gene
  ),
  category.names = c("Alkon", "Reynolds"),
  filename = NULL,  # Set filename to NULL to avoid saving to file
  output = TRUE,
  log = FALSE,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  
  # Numbers
  cex = 1,  # Increased size
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,  # Increased size
  cat.fontface = "bold",
  cat.default.pos = "text",
  cat.pos = c(-3, 3),
  cat.dist = c(0.03, 0.03),
  cat.fontfamily = "sans",
  
  # Title
  main = "Fibroblasts markers",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans"
)

# Display the Venn diagram
grid.draw(venn_plot)

# COMMAND ----------

# MAGIC %md
# MAGIC Now, lets check if the common markers are correct comparing with literature

# COMMAND ----------

#Save the common genes in variables
treg_common_genes <-  intersect(alkon_treg_markers$gene, reynolds_treg_markers$gene)
tcell_common_genes <- intersect(alkon_tc_markers$gene, reynolds_tc_markers$gene)
keratinocytes_common_genes <-intersect(alkon_kc_markers$gene, reynolds_kc_markers$gene)
melanocytes_common_genes <- intersect(alkon_melanocytes_markers$gene, reynolds_melanocytes_markers$gene)
fb_common_genes <- intersect(alkon_fb_markers$gene, reynolds_fb_markers$gene)

# COMMAND ----------

#According to literature (panglaodb + juanlu markers)
treg_markers <- c("FOXP3")  # Treg
tcell_markers <- c("CXCR4", "ICOS", "CCR7", "CD3G" ,"CD4", "PTPRC", "IL7R", "CD3E","TRBC2", "CD3D")  # T cell CD4
tcell_cd8_markers <- c("CD8A", "TBX21")  # T cell CD8 
ilc_markers <- c("KLRB1", "IL1R1")  # ILC
keratinocyte_markers <- c("LCE3C", "COL17A1", "KRT10", "KRT15")  # Keratinocyte
nk_markers <- c("NKG7", "KLRD1", "KLRF1")  # NK cell
melanocyte_markers <- c("MITF", "TYR", "DCT")  # Melanocyte
fb_markers <- c("VIM", "SERPINH1", "PDGFRB", "FAP")  # Fibroblast


# COMMAND ----------

# Compare the results with the literature and print with cell type names
treg_common <- treg_markers[treg_markers %in% treg_common_genes]
cat("Treg genes in literature:", treg_common, "\n")

tcell_common <- tcell_markers[tcell_markers %in% tcell_common_genes]
cat("T cell CD4 genes in literature:", tcell_common, "\n")

tcell_cd8_common <- tcell_cd8_markers[tcell_cd8_markers %in% tcell_common_genes]
cat("T cell genes in literature:", tcell_cd8_common, "\n")

nk_in_literature <- nk_markers[nk_markers %in% nk_common_genes]
cat("NK genes in literature:", nk_in_literature, "\n")

keratinocyte_common <- keratinocyte_markers[keratinocyte_markers %in% keratinocytes_common_genes]
cat("Keratinocyte genes in literature:", keratinocyte_common, "\n")

melanocyte_common <- melanocyte_markers[melanocyte_markers %in% melanocytes_common_genes]
cat("Melanocyte genes in literature:", melanocyte_common, "\n")

fb_common <- fb_markers[fb_markers %in% fb_common_genes]
cat("Fibroblast genes in literature:", fb_common, "\n")

# COMMAND ----------

cell_markers <- c(
  "FOXP3",  # Treg
  "CXCR5", "ICOS", "CCR7", "CD4",  # T cell CD4
  "CD8A", "TBX21",  # T cell CD8
  "KLRB1", "IL1R1",  # ILC
  "LCE3C", "COL17A1",  # Keratinocyte
  "CNN1", "TAGLN",  # Smooth muscle cell
  "PDGFRB", "ACTA2", "CSPG4",  # Pericyte
  "NKG7", "KLRD1", "KLRF1",  # NK cell
  "MITF", "TYR",  # Melanocyte
  "KIT",  # Mast cell
  "CD14",  # Macrophage
  "FAP",  # Fibroblast
  "CD34",  # Endothelial cell
  "CD86" , #DC
  "MS4A1", "CD27", "SPN", "CD19", "CD38"  # B cell
)

# COMMAND ----------

options(repr.plot.width=2800, repr.plot.height=1500)

DotPlot(alkon, features = cell_markers, group.by = "h_celltype_v4")

# COMMAND ----------

options(repr.plot.width=2500, repr.plot.height=1200)

DotPlot(reynolds, features = cell_markers, group.by = "h_celltype")
