# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC Reynolds markers
# MAGIC To compare the results when using only one dataset and using the merged dataset

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
    mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)) %>%
    mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%
    mutate(rank2 = (1 + avg_log2FC) * -log10(p_val_adj))
  
  # Assigning the 'gene' column values to the 'hgnc_symbol'
  res_rankedlist$hgnc_symbol <- res_rankedlist$gene
    display(res_rankedlist)

  gene_list <- left_join(res_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
  gene_list_entrez_id <- gene_list %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()
  
  geneList_entrez <- gene_list_entrez_id$rank
  names(geneList_entrez) <- gene_list_entrez_id$entrezgene_id
  
  geneList_entrez <- geneList_entrez[is.finite(geneList_entrez)]
  geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)
  geneList_entrez
  return(geneList_entrez)
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

#Fuction to generate a dotplot of the top pathway results
# Arguments: results, number of pathways (N with NES <0 and N with NES>0) and tittle of the plot
dotplot_pathway_generator <- function(results, N, title) {
  aux_pos <- arrange(results[which(results$NES > 0),], -NES)
  aux_neg <- arrange(results[which(results$NES < 0),], -NES)
  top_results <- unique(rbind(head(aux_pos, n = N), tail(aux_neg, n = N)))
  
  ggplot(top_results, aes(x = NES, y = reorder(Description, NES), size = setSize, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "red", high = "darkblue") +
    labs(x = "NES", y = "Pathway", size = "SetSize", color = "P-value adjust") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
}

# COMMAND ----------

#Fuction to generate a barplot of the top pathway results
# Arguments: results, number of pathways (N with NES <0 and N with NES>0) and tittle of the plot
barplot_gsea <- function(res, N, title) {
  aux_pos <- arrange(res[which(res$NES > 0),], -NES)
  aux_neg <- arrange(res[which(res$NES < 0),], -NES)
  GSEA_f <- unique(rbind(head(aux_pos, n = N), tail(aux_neg, n = N)))

  options(repr.plot.width = 1000, repr.plot.height = 1000, echo = FALSE)
  
  plot <- ggplot(GSEA_f, aes(NES, fct_reorder(Description, NES), fill = p.adjust)) + 
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
# MAGIC ##Keratinocytes (UNDIFF)

# COMMAND ----------

res_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_reynolds_bulk_undif_kc_de.xlsx")

# COMMAND ----------

geneList_entrez_kc <- prepare_gene_list(res_kc, G_list)

# COMMAND ----------

res_kc_0.05 <- get_enrichments(geneList_entrez_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

res_kc_0.05$Hallmark2

# COMMAND ----------

barplot_gsea(res_kc_0.05$Hallmark2, 20,"Keratinocytes Undiff - Reynolds- Hallmarks C2")

# COMMAND ----------

# res_kc_1 <- get_enrichments(geneList_entrez_kc, 1, H_t2g, C2_t2g)

# COMMAND ----------

set.seed(123)

kc_GO <- enrichGO(gene = names(geneList_entrez_kc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_kc_GO <- kc_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_kc_GO <- res_kc_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_kc_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(kc_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("Reynolds KC Undiff - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/pseudobulk/res_undif_kc_0.05_reynolds.xlsx")
write.xlsx(kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/pseudobulk/res_undif_kc_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (DIFF)

# COMMAND ----------

res_kc_dif <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/only_reynolds_bulk_dif_kc_de.xlsx")

# COMMAND ----------

geneList_entrez_kc_dif <- prepare_gene_list(res_kc_dif, G_list)

# COMMAND ----------

res_dif_kc_0.05 <- get_enrichments(geneList_entrez_kc_dif, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

res_dif_kc_0.05$Hallmark2

# COMMAND ----------

barplot_gsea(res_dif_kc_0.05$Hallmark2, 20,"Keratinocytes DIFF - Reynolds- Hallmarks C2")

# COMMAND ----------

# res_kc_1 <- get_enrichments(geneList_entrez_kc, 1, H_t2g, C2_t2g)

# COMMAND ----------

set.seed(123)

kc_GO_DIF <- enrichGO(gene = names(geneList_entrez_kc_dif),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_kc_GO_DIF <- kc_GO_DIF@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_kc_GO_DIF <- kc_GO_DIF %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_kc_GO_DIF)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(kc_GO_DIF, showCategory=25, label_format=100, font.size=9) + ggtitle("Reynolds KC DIFF - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_dif_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/pseudobulk/res_dif_kc_0.05_reynolds.xlsx")
write.xlsx(kc_GO_DIF, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/pseudobulk/res_dif_kc_reynolds_GO.xlsx")
