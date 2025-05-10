# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC alkon markers
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

#Fuction to generate a dotplot of the pathway results
# Arguments: results and tittle of the plot
dotplot_pathway_generator <- function(results, title) {
  ggplot(results, aes(x = NES, y = reorder(Description, NES), size = SetSize, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    scale_color_gradient(low = "red", high = "darkblue") +
    labs(x = "NES", y = "Pathway", size = "SetSize", color = "P-value adjust") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
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
# MAGIC ##Tcell

# COMMAND ----------

res_tcell <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_tcell_LvsHC_bulk_v2.xlsx")

# COMMAND ----------

geneList_entrez_Tcell <- prepare_gene_list(res_tcell, G_list)

# COMMAND ----------

# MAGIC %md
# MAGIC ###ENRICHMENT

# COMMAND ----------

res_tcell_0.05 <- get_enrichments(geneList_entrez_Tcell, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

set.seed(123)
tcell_GO <- enrichGO(gene = names(geneList_entrez_Tcell),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_tcell_GO <- tcell_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_tcell_GO <- res_tcell_GO %>%
  filter(p.adjust < 0.05)
display(filtered_res_tcell_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(tcell_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Alkon- Tcell - L vs HC - GO")

# COMMAND ----------

#Save results
# write.xlsx(res_tcell_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_tcell_0.05_alkon.xlsx")
write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_tcell_alkon_GO_v2.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

res_macro <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_macro_LvsHC_bulk_v2.xlsx")

# COMMAND ----------

geneList_entrez_macro <- prepare_gene_list(res_macro, G_list)

# COMMAND ----------

res_macro_0.05 <- get_enrichments(geneList_entrez_macro, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

set.seed(123)

macro_GO <- enrichGO(gene = names(geneList_entrez_macro),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_macro_GO <- macro_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_macro_GO <- res_macro_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_macro_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(macro_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Alkon- Macro - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
# write.xlsx(res_macro_0.05$Hallmark, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_macro_0.05_alkon_Hallmark.xlsx")
write.xlsx(macro_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_macro_alkon_GO_v2.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

res_treg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_treg_LvsHC_bulk_v2.xlsx")

# COMMAND ----------

geneList_entrez_treg <- prepare_gene_list(res_treg, G_list)

# COMMAND ----------

res_treg_0.05 <- get_enrichments(geneList_entrez_treg, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

treg_GO <- enrichGO(gene = names(geneList_entrez_treg),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_treg_GO <- treg_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_treg_GO <- res_treg_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_treg_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(treg_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Alkon- Treg - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

write.xlsx(treg_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_treg_alkon_GO_v2.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

res_fb <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_fb_LvsHC_bulk_v2.xlsx")

# COMMAND ----------

geneList_entrez_fb <- prepare_gene_list(res_fb, G_list)

# COMMAND ----------

res_fb_0.05 <- get_enrichments(geneList_entrez_fb, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

fb_GO <- enrichGO(gene = names(geneList_entrez_fb),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_fb_GO <- fb_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_fb_GO <- res_fb_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_fb_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(fb_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("Alkon- FB - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
# write.xlsx(res_fb_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/res_fb_0.05_alkon.xlsx")
write.xlsx(fb_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_fb_alkon_GO_v2.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

res_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/pseudobulk/alkon_kc_LvsHC_bulk_v2.xlsx")

# COMMAND ----------

geneList_entrez_kc <- prepare_gene_list(res_kc, G_list)

# COMMAND ----------

res_kc_0.05 <- get_enrichments(geneList_entrez_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_kc_0.05$Reactome, 10, "Top KC Alkon - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

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
barplot(kc_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("Alkon KC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_kc_0.05_alkon_v2.xlsx")
write.xlsx(kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Alkon/GSEA_individually/pseudobulk/res_kc_alkon_GO_v2.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC check

# COMMAND ----------

res_kc_0.05$Reactome[grepl("R-HSA-6785807", res_kc_0.05$Reactome$ID), ]
