# Databricks notebook source
# MAGIC %md
# MAGIC #Task 1: Dataset Integration
# MAGIC
# MAGIC
# MAGIC Ensure cell type homogeneity of the current annotations across datasets.
# MAGIC
# MAGIC - Integrate scRNA-seq datasets using **STACAS(Andreatta et al., 2024)** 

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))

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
  #remotes::install_bioc(pkg_name, lib=temp_library, upgrade=FALSE)
  remotes::install_github(pkg_name, lib = temp_library, upgrade=FALSE)
  for (x in list.files(temp_library)) {
    file.copy(
      file.path(temp_library, x),
      my_library,
      recursive=TRUE
    )
  }
}

# COMMAND ----------

if (!requireNamespace("STACAS")) install_from_github("carmonalab/STACAS")

# COMMAND ----------

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library("STACAS")
options(future.globals.maxSize = 1e9)

# COMMAND ----------

AR <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/MERGED_ARdatasets_TFM.rds")

# COMMAND ----------

Features(AR)

# COMMAND ----------

nfeatures <- 1000
ndim <- 20

# COMMAND ----------

#One step integration
  AR_integrated <- AR %>% SplitObject(split.by = "dataset") %>%
      Run.STACAS(dims = 1:ndim, anchor.features = nfeatures) %>%
      RunUMAP(dims = 1:ndim) 

DimPlot(AR_integrated, group.by = "dataset")

# COMMAND ----------

options(repr.plot.width=1600, repr.plot.height=1200)
DimPlot(AR_integrated,  group.by = c("dataset", "celltype"), label=TRUE)

# COMMAND ----------

#Save results
saveRDS(AR_integrated, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_STACAS_TFM.rds")

# COMMAND ----------

#SEMISUPERVISED
AR.semisup <- NormalizeData(AR) |>
    SplitObject(split.by = "dataset")|>
    Run.STACAS(cell.labels = "celltype")

AR.semisup <- RunUMAP(AR.semisup, dims = 1:30)


# COMMAND ----------

DimPlot(AR.semisup, group.by = c("dataset", "celltype"), label=TRUE)

# COMMAND ----------

saveRDS(AR.semisup, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_ssSTACAS_TFM.rds")

# COMMAND ----------

AR.semisup <- AR.semisup %>% RunUMAP(dims=1:ndim)

p1_ss <- DimPlot(AR.semisup, group.by = "dataset") +
  theme(aspect.ratio = 1) +
  ggtitle("Dataset after semi-supervised integration")
p2_ss <- DimPlot(AR.semisup, group.by = "celltype", label=T, label.size = 5) + 
  NoLegend() + theme(aspect.ratio = 1) + ggtitle("Cell labels after semi-supervised integration")

p1_ss | p2_ss

# COMMAND ----------

# MAGIC %md
# MAGIC ##Now it is performed again but **blocking the highest variable genes** (recomended by the paper)

# COMMAND ----------

#First step is to calculate highly variable genes
if (!requireNamespace("SignatuR")) install_from_github("carmonalab/SignatuR")
library(SignatuR)
#Retrieve full list of signatures for human
hs.sign <- GetSignature(SignatuR$Hs)

# COMMAND ----------

obj.list <- SplitObject(AR, split.by = "dataset")

# COMMAND ----------

my.genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
                        GetSignature(SignatuR$Hs$Compartments))

AR_integrated_blockList <- Run.STACAS(obj.list, genesBlockList = my.genes.blocklist,
                                          dims = 1:ndim, anchor.features = nfeatures)

# COMMAND ----------

AR_integrated_blockList
AR_integrated_blockList <- RunUMAP(AR_integrated_blockList, dims = 1:ndim)

# COMMAND ----------

DimPlot(AR_integrated_blockList , group.by = c("dataset", "celltype"), label=TRUE)

# COMMAND ----------

#SEMISUPERVISED
AR.semisup_block <- obj.list %>%
  Run.STACAS(dims = 1:ndim, anchor.features = nfeatures, cell.labels = "celltype", genesBlockList = my.genes.blocklist)

# COMMAND ----------

AR.semisup_block <- RunUMAP(AR.semisup_block, dims = 1:ndim)

# COMMAND ----------

DimPlot(AR.semisup_block , group.by = c("dataset", "celltype"), label=TRUE)

# COMMAND ----------

saveRDS(AR.semisup_block, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_ssSTACAS_block_TFM.rds")
saveRDS(AR_integrated_blockList, file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_STACAS_block_TFM.rds")
