# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC AR merged covariable dataset adjusted DEGs
# MAGIC
# MAGIC ## (NEW ANNOTATION CELLTYPIST)
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
    mutate(padj = ifelse(padj == 0, 1e-300, padj)) %>%
    mutate(rank = -log10(padj) * sign(log2FoldChange)) %>%
    mutate(rank2 = stat)
  
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

# DBTITLE 1,dotplot
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
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 1, family = "Arial")) +
    scale_size_continuous(range = c(2, 10))  # Adjust the range to make the dots bigger
}

# COMMAND ----------

library(viridis)
# Function to generate a dotplot of the top upregulated pathway results
# Arguments: results, number of pathways (N with NES >0) and title of the plot
dotplot_pathway_generator_up <- function(results, N, title) {
  aux_pos <- arrange(results[which(results$NES > 0),], -NES)
  top_results <- head(aux_pos, n = N)
  
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
# MAGIC ##Keratinocytes (Differentiated)

# COMMAND ----------

res_dif_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/bulk.dif.kc.de.cov.xlsx")

# COMMAND ----------

geneList_entrez_dif_kc <- prepare_gene_list(res_dif_kc, G_list)

# COMMAND ----------

res_dif_kc_0.05 <- get_enrichments(geneList_entrez_dif_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

options(repr.plot.width=1300, repr.plot.height=1200,  echo= F)

# COMMAND ----------

dotplot_pathway_generator_up(res_dif_kc_0.05$Reactome, 30, "Top 30 UP - Differentiated KC - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_dif_kc_0.05$Reactome, 15, "Top 30 - Differentiated KC - Reactome")

# COMMAND ----------

barplot_gsea(res_dif_kc_0.05$Reactome, 20, "Top - Differentiated KC AR - Reactome")

# COMMAND ----------

barplot_gsea(res_dif_kc_0.05$Hallmark1, 20, "Top - Differentiated KC AR - Hallmark H2")

# COMMAND ----------

barplot_gsea(res_dif_kc_0.05$Hallmark2, 20, "Top - Differentiated KC AR - Hallmark C2tg")

# COMMAND ----------

barplot_gsea(res_dif_kc_0.05$KEGG, 20, "Top - Differentiated KC AR - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

dif_kc_GO <- enrichGO(gene = names(geneList_entrez_dif_kc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_dif_kc_GO <- dif_kc_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_dif_kc_GO <- res_dif_kc_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_dif_kc_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(dif_kc_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("AR DIFF KC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_dif_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_dif_kc_0.05_AR.xlsx")
write.xlsx(dif_kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_dif_kc_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes (Undifferentiated)

# COMMAND ----------

res_undif_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/bulk.undif.kc.de.cov.xlsx")

# COMMAND ----------

geneList_entrez_undif_kc <- prepare_gene_list(res_undif_kc, G_list)

# COMMAND ----------

res_undif_kc_0.05 <- get_enrichments(geneList_entrez_undif_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_undif_kc_0.05$Reactome, 20, "Top - Undifferentiated KC AR - Reactome")

# COMMAND ----------

barplot_gsea(res_undif_kc_0.05$Hallmark1, 20, "Top - Undifferentiated KC AR - Hallmark H2")

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1000,  echo= F)

# COMMAND ----------

dotplot_pathway_generator_up(res_undif_kc_0.05$Reactome, 15, "Top 15 UP - Undifferentiated KC - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_undif_kc_0.05$Reactome, 15, "Top 30 - Undiff KC - Reactome (Merged-celltypist)")

# COMMAND ----------

barplot_gsea(res_undif_kc_0.05$Hallmark2, 20, "Top - Undifferentiated KC AR - Hallmark C2tg")

# COMMAND ----------

barplot_gsea(res_undif_kc_0.05$KEGG, 20, "Top - Undifferentiated KC AR - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

undif_kc_GO <- enrichGO(gene = names(geneList_entrez_undif_kc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_undif_kc_GO <- undif_kc_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_undif_kc_GO <- res_undif_kc_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_undif_kc_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(undif_kc_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("AR Undiff KC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_undif_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_undif_kc_0.05_AR.xlsx")
write.xlsx(undif_kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_undif_kc_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Tcell

# COMMAND ----------

res_tcell <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/bulk.tcell.de.cov.xlsx")

# COMMAND ----------

res_tcell

# COMMAND ----------

geneList_entrez_Tcell <- prepare_gene_list(res_tcell, G_list)

# COMMAND ----------

# MAGIC %md
# MAGIC ###ENRICHMENT

# COMMAND ----------

res_tcell_0.05 <- get_enrichments(geneList_entrez_Tcell, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_tcell_0.05$Hallmark1, 20, "Top 40 Tcell AR - Hallmark H")

# COMMAND ----------

options(repr.plot.width=1100, repr.plot.height=1000)

# COMMAND ----------

dotplot_pathway_generator_up(res_tcell_0.05$Reactome, 15, "Top 15 UP - Tcell (celltypist) - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_tcell_0.05$Reactome, 15, "Top 30 - Tcell - Reactome (Merged-celltypist)")

# COMMAND ----------

barplot_gsea(res_tcell_0.05$Reactome, 20, "Top 40 Tcell AR - Reactome")

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)

# COMMAND ----------

barplot_gsea(res_tcell_0.05$KEGG, 20, "Top 40 Tcell AR - KEGG")

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
barplot(tcell_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("AR- Tcell - L vs HC - GO")

# COMMAND ----------

# MAGIC %sh
# MAGIC mkdir /dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA

# COMMAND ----------

#Save results
write.xlsx(res_tcell_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_tcell_0.05_AR_.xlsx")
write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_tcell_AR_GO.xlsx")

# COMMAND ----------



# COMMAND ----------

# MAGIC %md
# MAGIC ##Th

# COMMAND ----------

res_th <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/bulk.th.de.cov.xlsx")

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

barplot_gsea(res_th_0.05$Hallmark1, 20, "Top 40 Th AR - Hallmark H")

# COMMAND ----------

options(repr.plot.width=1100, repr.plot.height=1000,  echo= F)

# COMMAND ----------

dotplot_pathway_generator_up(res_th_0.05$Reactome, 15, "Top 15 UP - Th (celltypist) - Reactome")

# COMMAND ----------

dotplot_pathway_generator(res_th_0.05$Reactome, 15,"Top 30 - Th - Reactome (Merged-celltypist)")

# COMMAND ----------

barplot_gsea(res_th_0.05$Reactome, 20, "Top 40 Th AR - Reactome")

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)

# COMMAND ----------

barplot_gsea(res_th_0.05$KEGG, 20, "Top 40 Th AR - KEGG")

# COMMAND ----------

set.seed(123)
th_GO <- enrichGO(gene = names(geneList_entrez_Th),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_th_GO <- th_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_th_GO <- res_th_GO %>%
  filter(p.adjust < 0.05)
display(filtered_res_th_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(tcell_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("AR- Th - L vs HC - GO")

# COMMAND ----------

#Save results
write.xlsx(res_th_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_th_0.05_AR_.xlsx")
write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_th_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

res_treg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/DEGs/bulk.treg.de.cov.xlsx")

# COMMAND ----------

geneList_entrez_treg <- prepare_gene_list(res_treg, G_list)

# COMMAND ----------

res_treg_0.05 <- get_enrichments(geneList_entrez_treg, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
barplot_gsea(res_treg_0.05$Reactome, 20, "Top Treg AR - Reactome")

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)

# COMMAND ----------



# COMMAND ----------

dotplot_pathway_generator(res_treg_0.05$Reactome, 15, "Top 30 - Treg - Reactome (Merged-celltypist)")

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot_gsea(res_treg_0.05$KEGG, 20, "Top Treg AR - KEGG")

# COMMAND ----------

options(repr.plot.width=1800, repr.plot.height=1000,  echo= F)
barplot_gsea(res_treg_0.05$Hallmark2, 20, "Top Treg AR - Hallmark C2")

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
barplot(treg_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("AR- Treg - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

write.xlsx(res_treg_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_treg_0.05_AR.xlsx")
write.xlsx(treg_GO, "/dbfs/mnt/sandbox/TFM_PAULA/merged_AR_celltypist_results/GSEA/res_treg_AR_GO.xlsx")
