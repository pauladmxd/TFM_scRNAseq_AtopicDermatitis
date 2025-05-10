# Batch-Adjusted Joint Analysis of scRNA-seq Datasets Reveals Underexplored Pathways in Atopic Dermatitis
This repository contains the notebooks for the joint analysis of multiple single-cell RNA-sequencing datasets related to atopic dermatitis (AD). The goal was to identify conserved disease signatures across patient cohorts and uncover cell-type specific pathways and gene expression patterns.

## Repository Structure
**Main folders:**

##### 00.Download_Preprocessing
      Notebooks related to the downloading of the datasets and their preprocessing for cell type homogeneization

##### 01.LvsHC_findmarkers
      Notebooks related to the analysis (Differential Expression Analysis (DEA) and Gene Set Enrichment Analysis (GSEA)) of the datasets independently using the FindAllMarkers approach.

##### 02.LvsHC_pseudobulk
      Notebooks related to the analysis (DEA and GSEA) of the datasets independently using the pseudo-bulk aggregated approach.

##### 03.Integration_Analysis
      Notebooks related to the exploration of the annotation consistency between datasets doign different integration approaches and evaluating them.

##### 04.Pseudobulk_Merged_Cov_Adjusted
      Notebooks related to the analysis (DEA and GSEA) of the merged datasets, adjusting or not the covariable and using the pseudo-bulk aggregated approach.
##### 05.Celltypist_Reannotation
###### ---- 01.Reannotation
               Notebooks related to the reannotation of Alkon dataset using celltypist package and Reynolds cell type annotation as a reference.
###### ---- 02.Reannotated_LvsHC_Pseudobulk
               Notebooks related to the analysis (DEA and GSEA) of the reannotated datasets (merged or independently) using the pseudo-bulk aggregated approach.
##### 06.Report_Figures
      Notebooks related to the generation of plots and figures for comparisons to be included in the paper report.

## Data Avaliability
**_Alkon et al, 2023:_** 10X Genomics data sets are publicly available via Gene Expression Omnibus GSE222840. Data from HC samples are available under Gene Expression Omnibus GSE173205 from a previously published data.

**_Reynolds et al, 2021:_** The raw sequencing data, expression count data with cell classifications are deposited at ArrayExpress: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8142.

##Key findings

- Pseudo-bulk differential expression improved robustness and biological relevance of results.
- Joint batch correction and detailed cell-type reannotation enabled detection of conserved disease signatures across AD datasets.
- Meaningful insights were derived even from modest, heterogeneous public data.
- Limitations in pathway databases and clinical metadata were acknowledged, but did not prevent identification of core biological signals.
- The study underscores the power of data reuse for generating new hypotheses without additional data collection.

##Contact
paula.delgadomanzano@almirall.com // paula.delgado01@estudiant.upf.edu