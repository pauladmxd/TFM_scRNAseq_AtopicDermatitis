# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC Reynolds markers
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
  res <- res %>%
    mutate(avg_log2FC = ifelse(cluster %in% c("healthy", "HC"), -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign
  
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
# MAGIC ##Tcell

# COMMAND ----------

res_tcell <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_LvsHC_tcell_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_Tcell <- prepare_gene_list(res_tcell, G_list)

# COMMAND ----------

# MAGIC %md
# MAGIC ###ENRICHMENT

# COMMAND ----------

res_tcell_0.05 <- get_enrichments(geneList_entrez_Tcell, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

res_tcell_1 <- get_enrichments(geneList_entrez_Tcell, 1, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_tcell_0.05$Hallmark1, 20, "Top 40 Tcell Reynolds - Hallmark H")

# COMMAND ----------

dotplot_pathway_generator(res_tcell_0.05$KEGG, 20, "Top Tcell Reynolds - KEGG")

# COMMAND ----------

dotplot_pathway_generator(res_tcell_0.05$Hallmark1,20, "Top Tcell Reynolds - Hallmark H")

# COMMAND ----------

barplot_gsea(res_tcell_0.05$Reactome, 20, "Top Tcell Reynolds - Reactome")

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
barplot(tcell_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Reynolds- Tcell - L vs HC - GO")

# COMMAND ----------

#Save results
write.xlsx(res_tcell_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_tcell_0.05_reynolds.xlsx")
write.xlsx(tcell_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_tcell_reynolds_GO.xlsx")

# COMMAND ----------

res_tcell_1$Reactome[grepl("R-HSA-6785807", res_tcell_1$Reactome$ID), ]

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

res_macro <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_macro <- prepare_gene_list(res_macro, G_list)

# COMMAND ----------

res_macro_0.05 <- get_enrichments(geneList_entrez_macro, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# res_macro_1 <- get_enrichments(geneList_entrez_macro, 1, H_t2g)

# COMMAND ----------

barplot_gsea(res_macro_0.05$Hallmark1, 20, "Top 40 Macro Reynolds - Hallmark")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
barplot_gsea(res_macro_0.05$Hallmark2, 20, "Top 40 Macro Reynolds - Hallmark C2")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
barplot_gsea(res_macro_0.05$Reactome, 20, "Top 40 Macro Reynolds - Reactome")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
barplot_gsea(res_macro_0.05$KEGG, 20, "Top 40 Macro Reynolds - KEGG")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
dotplot_pathway_generator(res_macro_0.05$Hallmark2, 20,  "Top Macro Reynolds - Hallmarks C2")

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
barplot(macro_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Reynolds- Macro - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_macro_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_macro_0.05_reynolds.xlsx")
write.xlsx(macro_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_macro_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

res_treg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_treg <- prepare_gene_list(res_treg, G_list)

# COMMAND ----------

res_treg_0.05 <- get_enrichments(geneList_entrez_treg, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# res_treg_1 <- get_enrichments(geneList_entrez_treg, 1, H_t2g)

# COMMAND ----------

options(repr.plot.width=1800, repr.plot.height=1000,  echo= F)
barplot_gsea(res_treg_0.05$Hallmark2, 30, "Top Treg Reynolds - Hallmark C2")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1000,  echo= F)
barplot_gsea(res_treg_0.05$Reactome, 30, "Top Treg Reynolds - Reactome")

# COMMAND ----------

  options(repr.plot.width = 1300, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_treg_0.05$Reactome, 20, "Top Treg Reyolds - Reactome")

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
barplot(treg_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("Reynolds- Treg - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

write.xlsx(res_treg_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_treg_0.05_reynolds.xlsx")
write.xlsx(treg_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_treg_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

res_fb <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

res_fb1 <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Alkon/alkon_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

display(res_fb)

# COMMAND ----------

  # res_fb <- res_fb %>%
  #   mutate(avg_log2FC = ifelse(cluster == "healthy", -abs(avg_log2FC), avg_log2FC)) #It is done because the analysis was done only with positive results and to be able to differentiate healthy and lesional I assign that sign
  
  # # Rank by pvalue
  # res_rankedlist <- res_fb %>%
  #   mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)) %>%
  #   mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%
  #   mutate(rank2 = (1 + avg_log2FC) * -log10(p_val_adj))
  # display(res_rankedlist)

# COMMAND ----------

geneList_entrez_fb <- prepare_gene_list(res_fb, G_list)

# COMMAND ----------

geneList_entrez_fb1 <- prepare_gene_list(res_fb1, G_list)

# COMMAND ----------

res_fb_0.05 <- get_enrichments(geneList_entrez_fb, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

res_fb_0.051 <- get_enrichments(geneList_entrez_fb1, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_fb_0.05$Hallmark1, 20, "Top 40 Fibroblasts Reynolds - Hallmark Ht2g")

# COMMAND ----------

barplot_gsea(res_fb_0.051$Hallmark1, 20, "Top 40 Fibroblasts Alkon - Hallmark Ht2g")

# COMMAND ----------

# MAGIC %md
# MAGIC **it is opposite in alkon**

# COMMAND ----------

  options(repr.plot.width = 1400, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_fb_0.05$Reactome, 30, "Top Fibroblasts Reynolds - Reactome")

# COMMAND ----------

  options(repr.plot.width = 1400, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_fb_0.05$Hallmark1, 30, "Top Fibroblasts Reynolds - Hallmark H")

# COMMAND ----------

barplot_gsea(res_fb_0.05$KEGG, 20, "Top 40 Fibroblasts Reynolds - KEGG")

# COMMAND ----------

barplot_gsea(res_fb_0.05$Reactome, 20, "Top 40 Fibroblasts Reynolds - Reactome")

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
barplot(fb_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("Reynolds- FB - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_fb_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_0.05_reynolds.xlsx")
write.xlsx(fb_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_fb_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

res_kc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_kc <- prepare_gene_list(res_kc, G_list)

# COMMAND ----------

res_kc_0.05 <- get_enrichments(geneList_entrez_kc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

res_kc_1 <- get_enrichments(geneList_entrez_kc, 1, H_t2g, C2_t2g)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot_gsea(res_kc_0.05$Reactome, 15, "Top 30 KC Reynolds - Reactome")

# COMMAND ----------

barplot_gsea(res_kc_0.05$Hallmark1, 20, "Top 40 Keratinocytes Reynolds - Hallmark Ht2g")

# COMMAND ----------

dotplot_pathway_generator(res_kc_0.05$Hallmark2, 20, "Top 40 KC Reynolds - Hallmark Ct2g")

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
barplot(kc_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("Reynolds KC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_kc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_kc_0.05_reynolds.xlsx")
write.xlsx(kc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_kc_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC Check

# COMMAND ----------

res_kc_1$Reactome[grepl("R-HSA-6785807", res_kc_1$Reactome$ID), ]

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

res_ilc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_ILC_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_ilc <- prepare_gene_list(res_ilc, G_list)

# COMMAND ----------

res_ilc_0.05 <- get_enrichments(geneList_entrez_ilc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# MAGIC %md
# MAGIC No results

# COMMAND ----------

set.seed(123)

ilc_GO <- enrichGO(gene = names(geneList_entrez_ilc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_ilc_GO <- ilc_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_ilc_GO <- res_ilc_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_ilc_GO)

# COMMAND ----------

#Save results
write.xlsx(ilc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_ilc_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

res_mastc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mast_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_mast <- prepare_gene_list(res_mastc, G_list)

# COMMAND ----------

res_mast_0.05 <- get_enrichments(geneList_entrez_mast, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

# MAGIC %md
# MAGIC No results

# COMMAND ----------

set.seed(123)

mast_GO <- enrichGO(gene = names(geneList_entrez_mast),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_mast_GO <- mast_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_mast_GO <- res_mast_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_mast_GO)

# COMMAND ----------

#Save results
write.xlsx(mast_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_mast_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

res_mono <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_mono <- prepare_gene_list(res_mono, G_list)

# COMMAND ----------

res_mono_0.05 <- get_enrichments(geneList_entrez_mono, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

barplot_gsea(res_mono_0.05$Hallmark1, 20, "Top 40 Mono Reynolds - Hallmark")

# COMMAND ----------

  options(repr.plot.width = 1400, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_mono_0.05$Hallmark2, 20, "Top 40 Mono Reynolds - Hallmark c2")

# COMMAND ----------

set.seed(123)

mono_GO <- enrichGO(gene = names(geneList_entrez_mono),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_mono_GO <- mono_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_mono_GO <- res_mono_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_mono_GO)

# COMMAND ----------

#Save results
write.xlsx(res_mono_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_mono_0.05_reynolds.xlsx")
write.xlsx(mono_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_mono_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##DC

# COMMAND ----------

res_dc <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/reynolds_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

geneList_entrez_dc <- prepare_gene_list(res_dc, G_list)

# COMMAND ----------

res_dc_0.05 <- get_enrichments(geneList_entrez_dc, 0.05, H_t2g, C2_t2g)

# COMMAND ----------

options(repr.plot.width = 1400, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_dc_0.05$Hallmark2, 20, "Top 40 DC Reynolds - Hallmark c2")

# COMMAND ----------

options(repr.plot.width = 1200, repr.plot.height = 1000, echo = FALSE)
dotplot_pathway_generator(res_dc_0.05$Reactome, 30, "Top 60 DC Reynolds - Reactome")

# COMMAND ----------

set.seed(123)

dc_GO <- enrichGO(gene = names(geneList_entrez_dc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_dc_GO <- dc_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_dc_GO <- res_dc_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_dc_GO)

# COMMAND ----------

#Save results
write.xlsx(res_dc_0.05, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_dc_0.05_reynolds.xlsx")
write.xlsx(dc_GO, "/dbfs/mnt/sandbox/TFM_PAULA/Reynolds/LvsHC/GSEA_individually/res_dc_reynolds_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##NK
