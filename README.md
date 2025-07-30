# Integrated Analysis of scRNA-seq Datasets Reveals Underexplored Pathways in Atopic Dermatitis
This repository contains the notebooks for the joint analysis of two single-cell RNA-sequencing datasets of atopic dermatitis (AD) skin samples. The goal was to identify conserved disease signatures across patient cohorts and uncover cell-type specific pathways and gene expression patterns particularly in keratinocytes. At the same time, different validated methods where tested to compare the differences in biological conclusions.

## Repository Structure
**Main folders:**

##### 00.Download_Preprocessing
      Notebooks related to the downloading of the datasets and their preprocessing for cell type homogeneization

##### 01.LvsHC_findmarkers

      Notebooks related to the analysis (Differential Expression Analysis (DEA) and Gene Set Enrichment Analysis (GSEA)) of the datasets independently using the FindAllMarkers approaches default (wilkoxon) and MAST tests.
###### ---- DEFAULT
               DEA using FinAllMarkers with default test.
###### ---- MAST  
               DEA using FinAllMarkers with MAST test.
##### 02.LvsHC_pseudobulk
      Notebooks related to the analysis (DEA and GSEA) of the datasets independently using the pseudo-bulk aggregated approach with two methods.
###### ---- DESeq2
               DEA using DESeq2 package.
###### ---- Limma 
               DEA using limma-voom.
##### 03.Integration_Analysis
      Notebooks related to the exploration of the annotation consistency between datasets doign different integration approaches and evaluating them.

##### 04.Pseudobulk_Merged_Cov_Adjusted

      Notebooks related to the analysis (DEA and GSEA) of the merged datasets, adjusting or not the covariable and using the pseudo-bulk aggregated approach with two methods.
###### ---- DESeq2
               DEA using DESeq2 package.
###### ---- Limma 
               DEA using limma-voom.
##### 05.Celltypist_Reannotation
###### ---- 01.Reannotation
               Notebooks related to the reannotation of Alkon dataset using celltypist package and Reynolds cell type annotation as a reference.
###### ---- 02.Reannotated_LvsHC_Pseudobulk
               Notebooks related to the analysis (DEA and GSEA) of the reannotated datasets (merged or independently) using the pseudo-bulk aggregated approach in:
               DESeq2 and Limma, as sepparated folders.
##### 06.Report_Figures
      Notebooks related to the generation of plots and figures for comparisons to be included in the paper report.

## Data Avaliability
**_Alkon et al, 2023:_** 10X Genomics data sets are publicly available via Gene Expression Omnibus GSE222840. Data from HC samples are available under Gene Expression Omnibus GSE173205 from a previously published data.

**_Reynolds et al, 2021:_** The raw sequencing data, expression count data with cell classifications are deposited at ArrayExpress: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8142.

## Key findings

- Pseudo-bulk differential expression improved robustness and biological relevance of results.
- Joint batch correction and detailed cell-type reannotation enabled detection of conserved disease signatures across AD datasets.
- Meaningful insights were derived even from modest, heterogeneous public data.
- Limited overlap between pseudobulk approaches due to methodologic differences when pathway enrichment analysis are performed.
- Limitations in pathway databases and clinical metadata were acknowledged, but did not prevent identification of core biological signals.
- GSVA is a robust alternative method to GSEA for functional enrichment, with limited overlap in the results and less number of significant ones.
- The study underscores the power of data reuse for generating new hypotheses without additional data collection.

## Contribution
Master thesis project performed in Almirall S.A. and Pompeu Fabra University.
Author: Paula Delgado Manzano
Supervisor: Juan Luis Trincado

## Contact
pauladelgadomanzano@gmailcom // paula.delgadomanzano@almirall.com // paula.delgado01@estudiant.upf.edu
