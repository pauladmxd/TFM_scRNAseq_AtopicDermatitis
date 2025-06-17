# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC AR merged covariable dataset adjusted DEGs
# MAGIC In order to interpret the DEA genes found and understand the biological context:
# MAGIC
# MAGIC Here I will perform the pathway enrichment with different databases (KEGG, Reactome...) in the common markers of the datasets.
# MAGIC
# MAGIC Also GO enrichment.

# COMMAND ----------

# MAGIC %sh
# MAGIC apt-get -y install libglpk-dev #Correct igraph - Need t load before, compatible for 14.3 LTS

# COMMAND ----------

# Load libraries
## Append the library folder
.libPaths(c("/dbfs/home/boriol@almirall.com/my_r_packages/bulkRNASeq_PBMCs_R4.3", .libPaths()))

# Load libraries
library(clusterProfiler)
library(ReactomePA) 
library(msigdbr)
library(DOSE)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)

.libPaths(c("/dbfs/home/jtrincado@almirall.com/my_r_packages/Seurat", .libPaths()))
library(openxlsx)
library(enrichplot)
library(ggplot2)

# COMMAND ----------

# MAGIC %md
# MAGIC ##Functions, G_list and HALLMARKS databases needed

# COMMAND ----------

prepare_gene_list <- function(res, G_list) {

  # Rank by pvalue
  res_rankedlist <- res %>%
    mutate(adj.P.Val = ifelse(adj.P.Val == 0, 1e-300, adj.P.Val)) %>%
    mutate(rank = -log10(adj.P.Val) * sign(logFC))
  
  # Assigning the 'gene' column values to the 'hgnc_symbol'
  res_rankedlist$hgnc_symbol <- res_rankedlist$gene
    display(res_rankedlist)

  gene_list <- left_join(res_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
  gene_list_entrez_id <- gene_list %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()
  # gene_list_entrez_id <- gene_list %>% dplyr::select(entrezgene_id, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()
  
  geneList_entrez <- gene_list_entrez_id$rank
  # geneList_entrez <- gene_list_entrez_id$rank2
  names(geneList_entrez) <- gene_list_entrez_id$entrezgene_id
  
  geneList_entrez <- geneList_entrez[is.finite(geneList_entrez)]
  geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)
  geneList_entrez
  return(geneList_entrez)
}

# COMMAND ----------

library(viridis)
#Fuction to generate a dotplot of the top pathway results
# Arguments: results, number of pathways (N with NES <0 and N with NES>0) and tittle of the plot
dotplot_pathway_generator <- function(results, N, title) {
  aux_pos <- arrange(results[which(results$NES > 0),], -NES)
  aux_neg <- arrange(results[which(results$NES < 0),], -NES)
  top_results <- unique(rbind(head(aux_pos, n = N), tail(aux_neg, n = N)))
  
  ggplot(top_results, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_c(option="plasma") +
    labs(x = "NES", y = "Pathway", size = "SetSize", color = "P-value adjust") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 13),
          plot.title = element_text(hjust = 1, family = "Arial")) +
    scale_size_continuous(range = c(2, 10))  # Adjust the range to make the dots bigger
}

# COMMAND ----------

get_enrichments <- function(gene_list, pval_cutoff, term2gene_hallmark, term2gene_hallmark2 = NULL) {
  set.seed(123)
  # KEGG enrichment
  kegg_result <- gseKEGG(gene_list, pvalueCutoff = pval_cutoff, organism = "hsa", verbose = FALSE, eps=0)
  res_kegg <- kegg_result@result
  
  # Reactome enrichment
  reactome_result <- gsePathway(gene_list, pAdjustMethod = "BH", pvalueCutoff = pval_cutoff, organism = "human", verbose = FALSE,  eps=0)
  res_reactome <- reactome_result@result
  
  # Hallmark enrichment
  hallmark_result <- GSEA(gene_list, TERM2GENE = term2gene_hallmark, pvalueCutoff = pval_cutoff,  eps=0)
  res_hallmark <- hallmark_result@result
  
  if (!is.null(term2gene_hallmark2)) {
    hallmark_result2 <- GSEA(gene_list, TERM2GENE = term2gene_hallmark2, pvalueCutoff = pval_cutoff,  eps=0)
    res_hallmark2 <- hallmark_result2@result
    return(list(KEGG = res_kegg, Reactome = res_reactome, Hallmark1 = res_hallmark, Hallmark2 = res_hallmark2))
  }
  
  return(list(KEGG = res_kegg, Reactome = res_reactome, Hallmark = res_hallmark))
}

# COMMAND ----------

#Fuction to generate a barplot of the top pathway results
# Arguments: results, number of pathways (N with NES <0 and N with NES>0) and tittle of the plot
barplot_gsea <- function(res, N, title) {
  aux_pos <- arrange(res[which(res$NES > 0),], -NES)
  aux_neg <- arrange(res[which(res$NES < 0),], -NES)
  GSEA_f <- unique(rbind(head(aux_pos, n = N), tail(aux_neg, n = N)))

  options(repr.plot.width = 1400, repr.plot.height = 1000, echo = FALSE)
  
  plot <- ggplot(GSEA_f, aes(NES, fct_reorder(Description, NES), fill = p.adjust), showCategory = 2*N) + 
    geom_col(orientation = 'y') + 
    scale_fill_continuous(low = 'red', high = 'blue', guide = guide_colorbar(reverse = TRUE)) + 
    theme_minimal() + ylab(NULL) + xlab("NES") +
    ggtitle(title)
  
  return(plot)
}

# COMMAND ----------

G_list <- readRDS("/dbfs/mnt/sandbox/RNASeq/PBMCs_IL4/pathways/G_list20240710.rds") #List with all translations to other id names
G_list <-  G_list %>% dplyr::filter(transcript_biotype == "protein_coding")

# COMMAND ----------

#HALLMARKS

library(msigdbr)
msigdbr_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
msigdbr_collections()

C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C2_t2g) #This collection includes gene sets curated from various sources such as online pathway databases and the biomedical literature. 

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g) # These gene sets summarize and represent specific well-defined biological states or processes. They are designed to reduce noise and redundancy, providing a clearer biological context.


CP_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>% 
  dplyr::select(gs_name, entrez_gene)
head(CP_t2g) #CP (Canonical Pathways): Gene sets from pathway databases representing canonical biological processes

C7_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C7_t2g) #Immunologic Signatures: Gene sets representing expression signatures of immune cell states, cell types, and perturbations

# COMMAND ----------

# MAGIC %md
# MAGIC ##DifKC Keratinocytes

# COMMAND ----------

res_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/Limma/alkon_limma_results_DifKC.xlsx")

# COMMAND ----------

geneList_entrez_kc <- prepare_gene_list(res_kc, G_list)

# COMMAND ----------

res_kc_0.05 <- get_enrichments(geneList_entrez_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

  options(repr.plot.width = 1400, repr.plot.height = 1300, echo = FALSE)

# COMMAND ----------

barplot_gsea(res_kc_0.05$Reactome, 15, "Top 30 - Diff KC - Merged Cov. Adj. - Reactome")

# COMMAND ----------

  options(repr.plot.width = 1500, repr.plot.height = 800, echo = FALSE)

# COMMAND ----------

dotplot_pathway_generator(res_kc_0.05$Reactome, 10, "Top 20 - Diff KC - Merged Cov. Adj. - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

# MAGIC %sh
# MAGIC  mkdir /dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/Limma

# COMMAND ----------

#Save results
write.xlsx(res_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/Limma/res_kc_0.05_AR.xlsx")
# write.xlsx(kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_results/GSEA/res_kc_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC check

# COMMAND ----------

res_kc_0.05$Reactome[grepl("R-HSA-6785807", res_kc_0.05$Reactome$ID), ]

# COMMAND ----------

# MAGIC %md
# MAGIC ##UndifKC Keratinocytes

# COMMAND ----------

res_undifkc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/Limma/alkon_limma_results_UndifKC.xlsx")

# COMMAND ----------

geneList_entrez_undifkc <- prepare_gene_list(res_undifkc, G_list)

# COMMAND ----------

res_undifkc_0.05 <- get_enrichments(geneList_entrez_undifkc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

  options(repr.plot.width = 1400, repr.plot.height = 1300, echo = FALSE)

# COMMAND ----------

barplot_gsea(res_undifkc_0.05$Reactome, 15, "Top 30 - Undif KC - Merged Cov. Adj. - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_undifkc_0.05$Reactome, 10, "Top 20 - Undif KC - Merged Cov. Adj. - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_undifkc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/Limma/res_undifkc_0.05_AR.xlsx")
# write.xlsx(kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_results/GSEA/res_kc_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC check

# COMMAND ----------

res_kc_0.05$Reactome[grepl("R-HSA-6785807", res_kc_0.05$Reactome$ID), ]

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tc

# COMMAND ----------

res_tc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/Limma/alkon_limma_results_Tc.xlsx")

# COMMAND ----------

res_tc

# COMMAND ----------

geneList_entrez_Tc <- prepare_gene_list(res_tc, G_list)

# COMMAND ----------

# MAGIC %md
# MAGIC ###ENRICHMENT

# COMMAND ----------

res_tc_0.05 <- get_enrichments(geneList_entrez_Tc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_tc_0.05$Reactome, 20, "Top 40 Tc AR - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_tc_0.05$Reactome, 15, "Top 30 - Tc - Merged Cov. Adj. - Reactome")

# COMMAND ----------

#Save results
write.xlsx(res_tc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/Limma/res_tc_0.05_AR_.xlsx")
# write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_results/GSEA/res_tcell_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Th

# COMMAND ----------

res_th <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/Limma/alkon_limma_results_Th.xlsx")

# COMMAND ----------

res_th

# COMMAND ----------

geneList_entrez_Th <- prepare_gene_list(res_th, G_list)

# COMMAND ----------

# MAGIC %md
# MAGIC ###ENRICHMENT

# COMMAND ----------

res_th_0.05 <- get_enrichments(geneList_entrez_Th, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_th_0.05$Reactome, 20, "Top 40 Th AR - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_th_0.05$Reactome, 15, "Top 30 - Th - Merged Cov. Adj. - Reactome")

# COMMAND ----------

#Save results
write.xlsx(res_th_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/Limma/res_th_0.05_AR_.xlsx")
# write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_results/GSEA/res_tcell_AR_GO.xlsx")
