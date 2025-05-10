# Databricks notebook source
# MAGIC %md
# MAGIC #Task 1: Dataset Integration
# MAGIC ##1.2: Evaluation and validation
# MAGIC
# MAGIC Ensure cell type homogeneity of the current annotations across datasets.
# MAGIC
# MAGIC Here the **evaluation of the integration methods** will be performed.

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

# if (!requireNamespace("STACAS")) install_from_github("carmonalab/STACAS")
if (!requireNamespace("scIntegrationMetrics")) install_from_github("carmonalab/scIntegrationMetrics")

# COMMAND ----------

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
# library("STACAS")
library("scIntegrationMetrics")
library(openxlsx)
options(future.globals.maxSize = 1e9)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Batch mixing evaluation
# MAGIC Batch mixing measures whether similar cells from different batches are well mixed after integration. Frequently used metrics of batch mixing are entropy, kBET, and integration LISI (iLISI)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Biological variance preservation evaluation
# MAGIC Preservation of biological variance can be quantified by how close to each other cells of the same type are, and how separated from each other cells of different types are in the joint integrated embeddings. Commonly-used metrics include average silhouette width (ASW), average Rand index (ARI), and cluster LISI (cLISI)

# COMMAND ----------

# MAGIC %md
# MAGIC ### RPCA

# COMMAND ----------

integratedAR_RPCA <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_RPCA_TFM.rds")

# COMMAND ----------

integratedAR_RPCA

# COMMAND ----------

identical(integratedAR_RPCA[["RNA"]]@counts@x, integratedAR_RPCA[["RNA"]]@data@x)
#It is already scaled and normalized
VariableFeatures(integratedAR_RPCA)
#FindVariableFeatures already done

# COMMAND ----------

?getIntegrationMetrics

# COMMAND ----------

RPCA_metrics <- getIntegrationMetrics(integratedAR_RPCA, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca.rpca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("norm_cLISI", "iLISI", "CiLISI")) 

unlist(RPCA_metrics)

# COMMAND ----------

RPCA_CiLISI <- getIntegrationMetrics(integratedAR_RPCA, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca.rpca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("CiLISI_means")) 

unlist(RPCA_CiLISI)

# COMMAND ----------

# MAGIC %md
# MAGIC CiLISI: Macro 0.07  KC 0.1  Prolif NaN  Treg 0.07  TC 0.02  DC 0.02  Fibroblasts 0.09  Endothelial blood NaN  Smooth Muscle NaN  LE 0.07  LC 0.01  Melanocytes 0.07  Plasmacytoid DC NaN  Pericyte 0  Mono 0  Sweat Gland NaN  B cells NaN  MastC 0.01  NK 0  nan NaN  VE NaN  ILC NaN  Schwann NaN  Plasma NaN  

# COMMAND ----------

# MAGIC %md
# MAGIC ##CCA

# COMMAND ----------

integratedAR_CCA <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_CCA_TFM.rds")

# COMMAND ----------

integratedAR_CCA

# COMMAND ----------

CCA_metrics <- getIntegrationMetrics(integratedAR_CCA, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca.cca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("norm_cLISI", "iLISI","CiLISI"))  

unlist(CCA_metrics)

# COMMAND ----------

CCA_CiLISI <- getIntegrationMetrics(integratedAR_CCA, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca.cca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("CiLISI_means"))  

unlist(CCA_CiLISI)

# COMMAND ----------

# MAGIC %md
# MAGIC CiLISI: Macro 0.08  KC 0.15  Prolif NaN  Treg 0.12  TC 0.03  DC 0.02  Fibroblasts 0.11  Endothelial blood NaN  Smooth Muscle NaN  LE 0.05  LC 0.01  Melanocytes 0.1  Plasmacytoid DC NaN  Pericyte 0  Mono 0  Sweat Gland NaN  B cells NaN  MastC 0  NK 0  nan NaN  VE NaN  ILC NaN  Schwann NaN  Plasma NaN  

# COMMAND ----------

# MAGIC %md
# MAGIC ##HARMONY

# COMMAND ----------

integratedAR_HARM <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_Harmony_TFM.rds")

# COMMAND ----------

integratedAR_HARM

# COMMAND ----------

# MAGIC %md
# MAGIC Error : cannot allocate vector of size 572.0 Gb
# MAGIC Error: cannot allocate vector of size 572.0 Gb That error happen when I try to calculate celltype_ASW

# COMMAND ----------

harm_metrics <- getIntegrationMetrics(integratedAR_HARM, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("norm_cLISI", "iLISI")) 

unlist(harm_metrics)

# COMMAND ----------

harm_metrics <- getIntegrationMetrics(integratedAR_HARM, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("CiLISI_means")) 

unlist(harm_metrics)

# COMMAND ----------

# MAGIC %md
# MAGIC CiLISI: Macro 0.06  KC 0.01  Prolif NaN  Treg 0.04  TC 0.02  DC 0.01  Fibroblasts 0.05  Endothelial blood NaN  Smooth Muscle NaN  LE 0.08  LC 0.01  Melanocytes 0.17  Plasmacytoid DC NaN  Pericyte 0  Mono 0  Sweat Gland NaN  B cells NaN  MastC 0.01  NK 0  nan NaN  VE NaN  ILC NaN  Schwann NaN  Plasma NaN 

# COMMAND ----------

# MAGIC %md
# MAGIC ##STACAS

# COMMAND ----------

integratedAR_stacas <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_STACAS_block_TFM.rds")

# COMMAND ----------

stacas_metrics <- getIntegrationMetrics(integratedAR_stacas, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("norm_cLISI", "iLISI", "CiLISI_means")) 

unlist(stacas_metrics)

# COMMAND ----------

# MAGIC %md
# MAGIC CiLISI: Macro 0.1  KC 0.1  Prolif NaN  Treg 0.06  TC 0.02  DC 0.01  Fibroblasts 0.13  Endothelial blood NaN  Smooth Muscle NaN  LE 0.15  LC 0.02  Melanocytes 0.2  Plasmacytoid DC NaN  Pericyte 0  Mono 0  Sweat Gland NaN  B cells NaN  MastC 0.01  NK 0  nan NaN  VE NaN  ILC NaN  Schwann NaN  Plasma NaN  
# MAGIC
# MAGIC iLISI  :    1.0826853 
# MAGIC        
# MAGIC norm_cLISI : 0.9936791

# COMMAND ----------

# MAGIC %md
# MAGIC ##ssSTACAS

# COMMAND ----------

integratedAR_ssstacas <- readRDS(file="/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integrated_AR_ssSTACAS_block_TFM.rds")

# COMMAND ----------

ssStacas_metrics <- getIntegrationMetrics(integratedAR_ssstacas, 
                                      meta.label = "celltype",
                                      meta.batch = "dataset",
                                      method.reduction = "pca",
                                      iLISI_perplexity = 20,
                                      cLISI_perplexity = 20,
                                      metrics = c("norm_cLISI", "iLISI", "CiLISI_means")) 

unlist(ssStacas_metrics)

# COMMAND ----------

# MAGIC %md
# MAGIC CiLISI: Macro 0.1  KC 0.12  Prolif NaN  Treg 0.05  TC 0.02  DC 0.01  Fibroblasts 0.12  Endothelial blood NaN  Smooth Muscle NaN  LE 0.15  LC 0.02  Melanocytes 0.21  Plasmacytoid DC NaN  Pericyte 0  Mono 0  Sweat Gland NaN  B cells NaN  MastC 0.01  NK 0  nan NaN  VE NaN  ILC NaN  Schwann NaN  Plasma NaN  
# MAGIC
# MAGIC iLISI : 1.0740949  
# MAGIC
# MAGIC norm_cLISI :  0.9944726 

# COMMAND ----------

# MAGIC %md
# MAGIC ##Save results

# COMMAND ----------

library(openxlsx)

# COMMAND ----------

# metrics_list <- list(#ssStacas_metrics = ssStacas_metrics, 
#                     #  stacas_metrics = stacas_metrics, 
#                     #  harm_metrics = harm_metrics,
#                     #  CCA_metrics = CCA_metrics,
#                      RPCA_metrics= RPCA_metrics)

write.xlsx(RPCA_metrics, file = "/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/RPCA_metrics_output.xlsx")

# COMMAND ----------

write.xlsx(CCA_metrics, file = "/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/CCA_metrics_output.xlsx")

# COMMAND ----------

write.xlsx(harm_metrics, file = "/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/harmony_metrics_output.xlsx")

# COMMAND ----------

write.xlsx(ssStacas_metrics, file = "/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/ssStacas_metrics_output.xlsx")

# COMMAND ----------

write.xlsx(stacas_metrics, file = "/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/Stacas_metrics_output.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##All the metrics

# COMMAND ----------

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(ggplot2)
library(openxlsx)

# COMMAND ----------

# Set options to avoid scientific notation
options(scipen = 999)

# Read the XLSX file
all_metrics <- list(
 ssSTACAS = read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integration_metrics.xlsx", sheet = 1),
 harmony = read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integration_metrics.xlsx", sheet = 2),
 RPCA = read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integration_metrics.xlsx", sheet = 3),
 CCA = read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/integrated_objects/integration_metrics.xlsx", sheet = 4)
 )
# Convert the column to numeric
all_metrics$ssSTACAS$CiLISI <- as.numeric(all_metrics$ssSTACAS$CiLISI) #It is done to have the decimal numbers correctly saved
all_metrics$harmony$CiLISI <- as.numeric(all_metrics$harmony$CiLISI)
all_metrics$RPCA$CiLISI <- as.numeric(all_metrics$RPCA$CiLISI)
all_metrics$CCA$CiLISI <- as.numeric(all_metrics$CCA$CiLISI)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Harmony

# COMMAND ----------

all_metrics$harmony$iLISI[1]

# COMMAND ----------

options(repr.plot.width=700, repr.plot.height=700)

# COMMAND ----------

#Keratinocyes CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = iLISI[1], y = CiLISI[3]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = iLISI[1], y = CiLISI[3], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = iLISI[1], y = CiLISI[3]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = iLISI[1], y = CiLISI[3], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = iLISI[1], y = CiLISI[3]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = iLISI[1], y = CiLISI[3], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = iLISI[1], y = CiLISI[3]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = iLISI[1], y = CiLISI[3], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "iLISI vs CiLISI Keratinocytes", x = "iLISI", y = "CiLISI")

# COMMAND ----------

#Keratinocyes CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[3]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[3], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[3]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[3], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[3]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[3], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[3]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[3], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs CiLISI Keratinocytes", x = "norm_cLISI", y = "CiLISI")

# COMMAND ----------

#Tcell CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[6]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[6], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[6]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[6], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[6]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[6], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[6]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[6], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs CiLISI Tcell", x = "norm_cLISI", y = "CiLISI")

# COMMAND ----------

#Fibroblasts CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[8]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[8], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[8]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[8], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[8]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[8], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[8]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[8], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs CiLISI Fibroblasts", x = "norm_cLISI", y = "CiLISI")

# COMMAND ----------

options(scipen = 999)
all_metrics$RPCA$CiLISI

# COMMAND ----------

#Melanocytes CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[13]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = CiLISI[13], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[13]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = CiLISI[13], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[13]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = CiLISI[13], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[13]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = CiLISI[13], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs CiLISI Melanocytes", x = "norm_cLISI", y = "CiLISI")

# COMMAND ----------

average_CiLISI <- list(
  harmony = mean(all_metrics$harmony$CiLISI, na.rm = TRUE),
  ssSTACAS = mean(all_metrics$ssSTACAS$CiLISI, na.rm = TRUE),
  CCA = mean(all_metrics$CCA$CiLISI, na.rm = TRUE),
  RPCA = mean(all_metrics$RPCA$CiLISI, na.rm = TRUE)
)

# COMMAND ----------

$harmony
[1] 0.03538462

$ssSTACAS
[1] 0.06230769

$CCA
[1] 0.05153846

$RPCA
[1] 0.04076923

# COMMAND ----------

#Average CiLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = average_CiLISI$harmony), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = average_CiLISI$harmony, label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = average_CiLISI$ssSTACAS), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = average_CiLISI$ssSTACAS, label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = average_CiLISI$CCA), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = average_CiLISI$CCA, label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = average_CiLISI$RPCA), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = average_CiLISI$RPCA, label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs average CiLISI ", x = "norm_cLISI", y = "avg_CiLISI")

# COMMAND ----------

#iLISI and normcLISI
ggplot() +
  geom_point(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = iLISI[1]), color = "blue") +
  geom_text(data = all_metrics$harmony, aes(x = norm_cLISI[1], y = iLISI[1], label = "harmony"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = iLISI[1]), color = "green") +
  geom_text(data = all_metrics$ssSTACAS, aes(x = norm_cLISI[1], y = iLISI[1], label = "ssSTACAS"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = iLISI[1]), color = "red") +
  geom_text(data = all_metrics$CCA, aes(x = norm_cLISI[1], y = iLISI[1], label = "CCA"), vjust = -1, size = 2.5) +
  geom_point(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = iLISI[1]), color = "cyan") +
  geom_text(data = all_metrics$RPCA, aes(x = norm_cLISI[1], y = iLISI[1], label = "RPCA"), vjust = -1, size = 2.5) +
  theme_minimal() +
  labs(title = "cLISI vs iLISI", x = "norm_cLISI", y = "iLISI")

# COMMAND ----------

# MAGIC %md
# MAGIC #Integration Evaluation: Randomized CiLISI
# MAGIC
# MAGIC Sayols suggests checking this parameter by performing a test where you randomize the labels of the cell types and see how this evaluates (the CiLISI should be very low, but this might be more for checking the cLISI).
# MAGIC
# MAGIC **1. Downsample the number of cells per identity class**
# MAGIC
# MAGIC **2. Randomize the batches for all cells.**
# MAGIC
# MAGIC **3.  Then, for about 100 randomizations, perform the integration and calculate the CiLISI.**
# MAGIC
# MAGIC **4. Use a loop in Databricks to save the CiLISIs, testing with about 5 to see how long it takes.**
# MAGIC
# MAGIC Sayols suggests creating subsets and not using all cells for this.
# MAGIC
# MAGIC Perform downsampling for this because it takes a long time and you need to test all methods.
# MAGIC
# MAGIC subset(x = pbmc, downsample = 100)
# MAGIC Check how long the integration takes.

# COMMAND ----------


