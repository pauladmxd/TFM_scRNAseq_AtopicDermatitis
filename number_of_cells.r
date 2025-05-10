# Databricks notebook source
.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))

library(Seurat)
library(dittoSeq)
library(dplyr)

# COMMAND ----------

alkon <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/ALKON_PROCESSED_TFM.rds")
liu <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/LIU_PROCESSED_TFM.rds")
reynolds <- readRDS("/dbfs/mnt/sandbox/TFM_PAULA/REYNOLDS_PROCESSED_TFM.rds")

# COMMAND ----------

print(paste("Nº cells alkon:", ncol(alkon)))
print(paste("Nº cells reynolds:", ncol(reynolds)))
print(paste("Nª cells liu:", ncol(liu)))
