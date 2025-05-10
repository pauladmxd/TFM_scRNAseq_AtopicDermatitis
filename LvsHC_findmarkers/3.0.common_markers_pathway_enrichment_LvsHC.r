# Databricks notebook source
# MAGIC %md
# MAGIC # Pathway enrichment in LvsHC common markers
# MAGIC In order to interpret the DEA genes found and understand the biological context:
# MAGIC
# MAGIC Here I will perform the pathway enrichment with different databases (KEGG, Reactome...) in the common markers of the datasets.
# MAGIC
# MAGIC Also GO enrichment.

# COMMAND ----------

# MAGIC %sh
# MAGIC apt-get -y install libglpk-dev #Correct igraph - Need t load before

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

G_list <- readRDS("/dbfs/mnt/sandbox/RNASeq/PBMCs_IL4/pathways/G_list20240710.rds") #List with all translations to other id names

# COMMAND ----------

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

ARL_tcell <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/ARL_Tcell_LvsHC_allmarkers.xlsx")
RL_tcell <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_Tcell_LvsHC_allmarkers.xlsx")

# COMMAND ----------

ARL_tcell <- ARL_tcell %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))
RL_tcell <- RL_tcell %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(ARL_tcell)

# COMMAND ----------

display(RL_tcell)

# COMMAND ----------

#Rank by pvalue
ARL_tcell_rankedlist <- ARL_tcell %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+avg_log2FC) * -log10(avg_pvalue))

# COMMAND ----------

display(ARL_tcell_rankedlist)

# COMMAND ----------

#R-L ranking
RL_tcell_rankedlist <- RL_tcell %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+ avg_log2FC)  * -log10(avg_pvalue))

# COMMAND ----------

display(RL_tcell_rankedlist)

# COMMAND ----------

colnames(ARL_tcell_rankedlist)

# COMMAND ----------

G_list <- readRDS("/dbfs/mnt/sandbox/RNASeq/PBMCs_IL4/pathways/G_list20240710.rds") #List with all translations to other id names

# COMMAND ----------

G_list <-  G_list %>% dplyr::filter(transcript_biotype == "protein_coding")

# COMMAND ----------

# Assigning the 'gene' column values to the 'hgnc_symbol'
ARL_tcell_rankedlist$hgnc_symbol <- ARL_tcell_rankedlist$gene
RL_tcell_rankedlist$hgnc_symbol <- RL_tcell_rankedlist$gene

# COMMAND ----------

#Merge by symbol and select relevant columns for gene list
#ARL
gene_list_ARL_Tcell <- left_join(ARL_tcell_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_ARL_Tcell <- gene_list_ARL_Tcell %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

#RL
gene_list_RL_Tcell <- left_join(RL_tcell_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_Tcell <- gene_list_RL_Tcell %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_ARL_Tcell  <- gene_list_entrez_id_ARL_Tcell$rank #add the rank to the geneList
geneList_entrez_RL_Tcell  <- gene_list_entrez_id_RL_Tcell$rank 

#entrez id as the names in geneList
names(geneList_entrez_ARL_Tcell) <- gene_list_entrez_id_ARL_Tcell$entrezgene_id 
names(geneList_entrez_RL_Tcell) <- gene_list_entrez_id_RL_Tcell$entrezgene_id 

#The same procedure but with the second type of ranking rank2
geneList_entrez_ARL_Tcell_1  <- gene_list_entrez_id_ARL_Tcell$rank2 #add the rank to the geneList
geneList_entrez_RL_Tcell_1  <- gene_list_entrez_id_RL_Tcell$rank2

names(geneList_entrez_ARL_Tcell_1) <- gene_list_entrez_id_ARL_Tcell$entrezgene_id 
names(geneList_entrez_RL_Tcell_1) <- gene_list_entrez_id_RL_Tcell$entrezgene_id

# COMMAND ----------

# MAGIC %md
# MAGIC Clean and sort the gene lists to use them in gse

# COMMAND ----------

# Clean and sort the ARL geneList
geneList_entrez_ARL_Tcell <- geneList_entrez_ARL_Tcell[is.finite(geneList_entrez_ARL_Tcell)]
geneList_entrez_ARL_Tcell <- sort(geneList_entrez_ARL_Tcell, decreasing = TRUE)

# Clean and sort the ARL geneList with rank2
geneList_entrez_ARL_Tcell_1 <- geneList_entrez_ARL_Tcell_1[is.finite(geneList_entrez_ARL_Tcell_1)]
geneList_entrez_ARL_Tcell_1 <- sort(geneList_entrez_ARL_Tcell_1, decreasing = TRUE)

geneList_entrez_ARL_Tcell
geneList_entrez_ARL_Tcell_1

# COMMAND ----------

# Clean and sort the R-L geneList
geneList_entrez_RL_Tcell <- geneList_entrez_RL_Tcell[is.finite(geneList_entrez_RL_Tcell)]
geneList_entrez_RL_Tcell <- sort(geneList_entrez_RL_Tcell, decreasing = TRUE)

# Clean and sort the R-L geneList with rank2
geneList_entrez_RL_Tcell_1 <- geneList_entrez_RL_Tcell_1[is.finite(geneList_entrez_RL_Tcell_1)]
geneList_entrez_RL_Tcell_1 <- sort(geneList_entrez_RL_Tcell_1, decreasing = TRUE)

geneList_entrez_RL_Tcell
geneList_entrez_RL_Tcell_1

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

set.seed(123)
tcell_ARL_KEGG <- gseKEGG(geneList_entrez_ARL_Tcell,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          verbose      = FALSE,
                          minGSSize = 10,
                          maxGSSize = 500)
                           #,scoreType = "pos") # because I only have upregulated genes

# Extract the result from the KEGG enrichment analysis
res_tcell_ARL_KEGG <- tcell_ARL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_tcell_ARL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_tcell_ARL_KEGG)
}

# COMMAND ----------

set.seed(123)

tcell_ARL_KEGG_1 <- gseKEGG(geneList_entrez_ARL_Tcell_1,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          verbose      = FALSE,
                          minGSSize = 10,
                          maxGSSize = 500)                       
                          #,scoreType = "pos") # because I only have upregulated genes

# Extract the result from the KEGG enrichment analysis
res_tcell_ARL_KEGG_1 <- tcell_ARL_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_tcell_ARL_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_tcell_ARL_KEGG_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC no significant results with rank2

# COMMAND ----------

set.seed(123)

#Pathway enrichment with KEGG in R-L common markers
tcell_RL_KEGG <- gseKEGG(geneList_entrez_RL_Tcell,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)
                      #,scoreType = "pos", #because I only have upregulated genes


res_tcell_RL_KEGG <- tcell_RL_KEGG@result
# Check if there are any enriched pathways
if (nrow(res_tcell_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_tcell_RL_KEGG)
}

# COMMAND ----------

set.seed(123)
#Pathway enrichment with KEGG in R-L common markers
tcell_RL_KEGG_1 <- gseKEGG(geneList_entrez_RL_Tcell_1,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)
                      #scoreType = "pos") #because I only have upregulated genes


res_tcell_RL_KEGG_1 <- tcell_RL_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_tcell_RL_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_tcell_RL_KEGG_1)
}

# COMMAND ----------

colnames(res_tcell_ARL_KEGG)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
dotplot(tcell_RL_KEGG, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Tcell - Lvs HC - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#ARL Reactome Enrichment
tcell_ARL_Reactome <- gsePathway(geneList_entrez_ARL_Tcell,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_tcell_ARL_Reactome <- tcell_ARL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_tcell_ARL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_ARL_Reactome)
}

# COMMAND ----------


set.seed(123)
#R-L Reactome enrichment
tcell_RL_Reactome <- gsePathway(geneList_entrez_RL_Tcell,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_tcell_RL_Reactome <- tcell_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_tcell_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(tcell_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Tcell - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

# MAGIC %md
# MAGIC R-L hallmarks enrichment

# COMMAND ----------

set.seed(123)
tcell_RL_HALLMARKS <- GSEA(geneList_entrez_RL_Tcell, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)
res_tcell_RL_HALLMARKS <- tcell_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_tcell_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)

#HALLMARKS C2
tcell_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_Tcell, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                            minGSSize=10,
                          maxGSSize=500)
res_tcell_RL_HALLMARKS_C2 <- tcell_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_tcell_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_RL_HALLMARKS_C2)
}


# COMMAND ----------

# MAGIC %md
# MAGIC RL using rank2

# COMMAND ----------

set.seed(123)
#HALLMARKS C2 but with type ranking rank2
tcell_RL_HALLMARKS_C2_1 <- GSEA(geneList_entrez_RL_Tcell_1, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)
res_tcell_RL_HALLMARKS_C2_1 <- tcell_RL_HALLMARKS_C2_1@result

# Check if there are any enriched pathways
if (nrow(res_tcell_RL_HALLMARKS_C2_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_RL_HALLMARKS_C2_1)
}


# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(tcell_RL_HALLMARKS, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Tcell - L vs HC - hallmarks")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(tcell_RL_HALLMARKS_C2, showCategory= 25, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Tcell - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC A-R-L HALLMARKS enrichment

# COMMAND ----------

set.seed(123)
tcell_ARL_HALLMARKS <- GSEA(geneList_entrez_ARL_Tcell, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize = 10,
                          maxGSSize=500)

res_tcell_ARL_HALLMARKS <- tcell_ARL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_tcell_ARL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_ARL_HALLMARKS)
}


# COMMAND ----------

set.seed(123)
tcell_ARL_HALLMARKS_C2 <- GSEA(geneList_entrez_ARL_Tcell, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)

res_tcell_ARL_HALLMARKS_C2 <- tcell_ARL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_tcell_ARL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_tcell_ARL_HALLMARKS_C2)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No pathways in ARL

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

tcell_ARL_GO <- enrichGO(gene = names(geneList_entrez_ARL_Tcell),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_tcell_ARL_GO <- tcell_ARL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_tcell_ARL_GO <- res_tcell_ARL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_tcell_ARL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(tcell_ARL_GO, showCategory=25, label_format=50, font.size=9) + ggtitle("ARL- Tcell - L vs HC - GO")

# COMMAND ----------

tcell_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_Tcell),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_tcell_RL_GO <- tcell_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_tcell_RL_GO <- res_tcell_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(res_tcell_RL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(tcell_RL_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("RL - Tcell - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results RL
write.xlsx(res_tcell_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_PE/pathway_enrich_tcell_RL_KEGG.xlsx")
write.xlsx(res_tcell_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_PE/pathway_enrich_tcell_RL_REACTOME.xlsx")
write.xlsx(res_tcell_RL_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_PE/pathway_enrich_tcell_RL_HALLMARKS.xlsx")
write.xlsx(res_tcell_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_PE/pathway_enrich_tcell_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_tcell_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/RL_PE/functional_enrich_tcell_RL_GO.xlsx")

#Save results ARL
write.xlsx(res_tcell_ARL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/ARL_PE/pathway_enrich_tcell_ARL_KEGG.xlsx")
write.xlsx(res_tcell_ARL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/ARL_PE/pathway_enrich_tcell_ARL_REACTOME.xlsx")
write.xlsx(res_tcell_ARL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Tcell/ARL_PE/functional_enrich_tcell_ARL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##NK

# COMMAND ----------

RL_nk <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/NK/RL_nk_LvsHC_allmarkers.xlsx")

# COMMAND ----------

# Here I do not have healthy markers so  skip this part
# RL_nk <- RL_nk %>%
#   mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(RL_nk)

# COMMAND ----------

# MAGIC %md
# MAGIC Just to take into account, NK only have DEA upregulated in lesional samples. It is something interesting. Also, in Alkon, it was seen NK were only sequenced in AD patient samples.

# COMMAND ----------

#Rank
# Replace 0 p-values with a small number to be able to calculate the rank properly
RL_nk <- RL_nk %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue))

RL_nk_rankedlist <- RL_nk %>%
  mutate(rank = -log10(avg_pvalue))%>%
  mutate(rank2 = -log10(avg_pvalue) * avg_log2FC)

# COMMAND ----------

# MAGIC %md
# MAGIC Most of the genes have same adj p value. Also there are few genes

# COMMAND ----------

RL_nk_rankedlist$hgnc_symbol <- RL_nk_rankedlist$gene
gene_list_RL_NK <- left_join(RL_nk_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_NK <- gene_list_RL_NK %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_RL_NK  <- gene_list_entrez_id_RL_NK$rank 

#entrez id as the names in geneList
names(geneList_entrez_RL_NK) <- gene_list_entrez_id_RL_NK$entrezgene_id

# COMMAND ----------

# Clean and sort the geneList
geneList_entrez_RL_NK <- geneList_entrez_RL_NK[is.finite(geneList_entrez_RL_NK)]

geneList_entrez_RL_NK <- sort(geneList_entrez_RL_NK, decreasing = TRUE)

geneList_entrez_RL_NK

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

#Pathway enrichment with KEGG in NK cell common markers
nk_RL_KEGG <- gseKEGG(geneList_entrez_RL_NK,
                      organism     = 'hsa',
                      pvalueCutoff = 1,
                      verbose      = FALSE,
                      scoreType = "pos")

res_nk_RL_KEGG <- nk_RL_KEGG@result

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough statistical power to perform the pathway enrichment.

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

#R-L Reactome enrichment
nk_RL_Reactome <- gsePathway(geneList_entrez_RL_NK,
                pvalueCutoff = 1,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                scoreType = "pos")
res_nk_RL_Reactome <- nk_RL_Reactome@result

#display(res_nk_RL_Reactome)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough statistical power to perform the pathway enrichment.

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

nk_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_NK, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 1, 
                          scoreType = "pos")
res_nk_RL_HALLMARKS_C2 <- nk_RL_HALLMARKS_C2@result

#display(res_nk_RL_HALLMARKS_C2)

# COMMAND ----------

# MAGIC %md
# MAGIC Not enough statistical power to perform the pathway enrichment. I do not even try with H2g

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

nk_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_NK),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_nk_RL_GO <- nk_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_nk_RL_GO <- res_nk_RL_GO %>%
  filter(p.adjust < 0.05)

# No significant results


# COMMAND ----------

# MAGIC %md
# MAGIC Again not egough statistical power, no enriched GO term with low pvalue

# COMMAND ----------

# MAGIC %md
# MAGIC ##ILC

# COMMAND ----------

RL_ilc <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/RL_ilc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
RL_ilc <- RL_ilc %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(RL_ilc)

# COMMAND ----------

#Rank
# Replace 0 p-values with a small number to be able to calculate the rank properly
RL_ilc <- RL_ilc %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue))

RL_ilc_rankedlist <- RL_ilc %>%
  mutate(rank = -log10(avg_pvalue))%>%
  mutate(rank2 = -log10(avg_pvalue) * avg_log2FC)

# COMMAND ----------

RL_ilc_rankedlist$hgnc_symbol <- RL_ilc_rankedlist$gene
gene_list_RL_ILC <- left_join(RL_ilc_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_ILC <- gene_list_RL_ILC %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_RL_ILC <- gene_list_entrez_id_RL_ILC$rank 

#entrez id as the names in geneList
names(geneList_entrez_RL_ILC) <- gene_list_entrez_id_RL_ILC$entrezgene_id

# COMMAND ----------

# Clean and sort the geneList
geneList_entrez_RL_ILC <- geneList_entrez_RL_ILC[is.finite(geneList_entrez_RL_ILC)]

geneList_entrez_RL_ILC <- sort(geneList_entrez_RL_ILC, decreasing = TRUE)

geneList_entrez_RL_ILC

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

# Pathway enrichment with KEGG in ILC cell common markers
ilc_RL_KEGG <- gseKEGG(geneList_entrez_RL_ILC,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_ilc_RL_KEGG <- ilc_RL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_ilc_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_ilc_RL_KEGG)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(ilc_RL_KEGG, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - ILC - L vs HC - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#RL Reactome Enrichment for ILC
ilc_RL_Reactome <- gsePathway(geneList_entrez_RL_ILC,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500,
                eps=0)
res_ilc_RL_Reactome <- ilc_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_ilc_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_ilc_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(ilc_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - ILC - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

set.seed(123)
ilc_RL_HALLMARKS <- GSEA(geneList_entrez_RL_ILC, 
                        TERM2GENE = H_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500)
res_ilc_RL_HALLMARKS <- ilc_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_ilc_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_ilc_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
ilc_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_ILC, 
                            TERM2GENE = C2_t2g, 
                            pvalueCutoff = 0.05, 
                            minGSSize = 10,
                            maxGSSize = 500) 
res_ilc_RL_HALLMARKS_C2 <- ilc_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_ilc_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_ilc_RL_HALLMARKS_C2)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(ilc_RL_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - ILC - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

ilc_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_ILC),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_ilc_RL_GO <- ilc_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_ilc_RL_GO <- res_ilc_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_ilc_RL_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(ilc_RL_GO, showCategory=50, label_format=100, font.size=9) + ggtitle("RL- ILC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_ilc_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/pathway_enrich_ilc_RL_KEGG.xlsx")
write.xlsx(res_ilc_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/pathway_enrich_ilc_RL_REACTOME.xlsx")
write.xlsx(res_ilc_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/pathway_enrich_ilc_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_ilc_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/ILC/functional_enrich_ilc_RL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Macrophages

# COMMAND ----------

ARL_macro <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/ARL_macro_LvsHC_allmarkers.xlsx")
RL_macro <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_macro_LvsHC_allmarkers.xlsx")

# COMMAND ----------

ARL_macro <- ARL_macro %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))
RL_macro <- RL_macro %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(ARL_macro)

# COMMAND ----------

#Rank by pvalue ARL
ARL_macro_rankedlist <- ARL_macro %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+avg_log2FC) * -log10(avg_pvalue))

# COMMAND ----------

#R-L ranking
RL_macro_rankedlist <- RL_macro %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+ avg_log2FC)  * -log10(avg_pvalue))

# COMMAND ----------

display(RL_macro_rankedlist)

# COMMAND ----------

# Assigning the 'gene' column values to the 'hgnc_symbol'
ARL_macro_rankedlist$hgnc_symbol <- ARL_macro_rankedlist$gene
RL_macro_rankedlist$hgnc_symbol <- RL_macro_rankedlist$gene

# COMMAND ----------

#Merge by symbol and select relevant columns for gene list
#ARL
gene_list_ARL_macro <- left_join(ARL_macro_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_ARL_macro <- gene_list_ARL_macro %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()
#RL
gene_list_RL_macro <- left_join(RL_macro_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_macro <- gene_list_RL_macro %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_ARL_macro  <- gene_list_entrez_id_ARL_macro$rank #add the rank to the geneList
geneList_entrez_RL_macro  <- gene_list_entrez_id_RL_macro$rank 

#entrez id as the names in geneList
names(geneList_entrez_ARL_macro) <- gene_list_entrez_id_ARL_macro$entrezgene_id 
names(geneList_entrez_RL_macro) <- gene_list_entrez_id_RL_macro$entrezgene_id 

#The same procedure but with the second type of ranking rank2
geneList_entrez_ARL_macro_1  <- gene_list_entrez_id_ARL_macro$rank2 #add the rank to the geneList
geneList_entrez_RL_macro_1  <- gene_list_entrez_id_RL_macro$rank2

names(geneList_entrez_ARL_macro_1) <- gene_list_entrez_id_ARL_macro$entrezgene_id 
names(geneList_entrez_RL_macro_1) <- gene_list_entrez_id_RL_macro$entrezgene_id

# COMMAND ----------

# MAGIC
# MAGIC %md
# MAGIC Clean and sort the gene lists to use them in gse

# COMMAND ----------

# Clean and sort the ARL geneList
geneList_entrez_ARL_macro <- geneList_entrez_ARL_macro[is.finite(geneList_entrez_ARL_macro)]
geneList_entrez_ARL_macro <- sort(geneList_entrez_ARL_macro, decreasing = TRUE)

# Clean and sort the ARL geneList with rank2
geneList_entrez_ARL_macro_1 <- geneList_entrez_ARL_macro_1[is.finite(geneList_entrez_ARL_macro_1)]
geneList_entrez_ARL_macro_1 <- sort(geneList_entrez_ARL_macro_1, decreasing = TRUE)

geneList_entrez_ARL_macro
geneList_entrez_ARL_macro_1

# COMMAND ----------

# Clean and sort the R-L geneList
geneList_entrez_RL_macro <- geneList_entrez_RL_macro[is.finite(geneList_entrez_RL_macro)]
geneList_entrez_RL_macro <- sort(geneList_entrez_RL_macro, decreasing = TRUE)

# Clean and sort the R-L geneList with rank2
geneList_entrez_RL_macro_1 <- geneList_entrez_RL_macro_1[is.finite(geneList_entrez_RL_macro_1)]
geneList_entrez_RL_macro_1 <- sort(geneList_entrez_RL_macro_1, decreasing = TRUE)

geneList_entrez_RL_macro
geneList_entrez_RL_macro_1

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

set.seed(123)
macro_ARL_KEGG <- gseKEGG(geneList_entrez_ARL_macro,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          verbose      = FALSE,
                          minGSSize = 10,
                          maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_macro_ARL_KEGG <- macro_ARL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_macro_ARL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_macro_ARL_KEGG)
}

# COMMAND ----------

set.seed(123)

#Pathway enrichment with KEGG in R-L common markers
macro_RL_KEGG <- gseKEGG(geneList_entrez_RL_macro,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)


res_macro_RL_KEGG <- macro_RL_KEGG@result
# Check if there are any enriched pathways
if (nrow(res_macro_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_macro_RL_KEGG)
}

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
dotplot(macro_RL_KEGG, x = "NES", showCategory= 30, size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - macro - L vs HC- KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#ARL Reactome Enrichment
macro_ARL_Reactome <- gsePathway(geneList_entrez_ARL_macro,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_macro_ARL_Reactome <- macro_ARL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_macro_ARL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_macro_ARL_Reactome)
}

# COMMAND ----------

set.seed(123)
#R-L Reactome enrichment
macro_RL_Reactome <- gsePathway(geneList_entrez_RL_macro,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500,
                eps=0)
res_macro_RL_Reactome <- macro_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_macro_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_macro_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(macro_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Macro - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

# MAGIC %md
# MAGIC A-R-L HALLMARKS enrichment

# COMMAND ----------

set.seed(123)
macro_ARL_HALLMARKS <- GSEA(geneList_entrez_ARL_macro, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize = 10,
                          maxGSSize=500)

res_macro_ARL_HALLMARKS <- macro_ARL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_macro_ARL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_macro_ARL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
macro_ARL_HALLMARKS_C2 <- GSEA(geneList_entrez_ARL_macro, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.2, 
                          minGSSize=10,
                          maxGSSize=500)

res_macro_ARL_HALLMARKS_C2 <- macro_ARL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_macro_ARL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the  enrichment analysis result 
  display(res_macro_ARL_HALLMARKS_C2)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results in ARL Treg

# COMMAND ----------

# MAGIC %md
# MAGIC R-L hallmarks enrichment

# COMMAND ----------

set.seed(123)
macro_RL_HALLMARKS <- GSEA(geneList_entrez_RL_macro, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)
res_macro_RL_HALLMARKS <- macro_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_macro_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_macro_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)

#HALLMARKS C2
macro_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_macro, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500,
                          eps=0)#To control low pvalue warning, it does not use approximations
res_macro_RL_HALLMARKS_C2 <- macro_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_macro_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_macro_RL_HALLMARKS_C2)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(macro_RL_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Macro - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

macro_ARL_GO <- enrichGO(gene = names(geneList_entrez_ARL_macro),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_macro_ARL_GO <- macro_ARL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_macro_ARL_GO <- res_macro_ARL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_macro_ARL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(macro_ARL_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("ARL- Macro - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC Low counts

# COMMAND ----------

macro_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_macro),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_macro_RL_GO <- macro_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_macro_RL_GO <- res_macro_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(res_macro_RL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(macro_RL_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("RL - Macro - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results RL
write.xlsx(res_macro_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_PE/pathway_enrich_macro_RL_KEGG.xlsx")
write.xlsx(res_macro_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_PE/pathway_enrich_macro_RL_REACTOME.xlsx")
write.xlsx(res_macro_RL_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_PE/pathway_enrich_macro_RL_HALLMARKS.xlsx")
write.xlsx(res_macro_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_PE/pathway_enrich_macro_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_macro_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/RL_PE/functional_enrich_macro_RL_GO.xlsx")

#Save results ARL
write.xlsx(res_macro_ARL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Macrophage/ARL_PE/functional_enrich_macro_ARL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Monocytes

# COMMAND ----------

RL_mono <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/RL_mono_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
RL_mono <- RL_mono %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(RL_mono)

# COMMAND ----------

#Rank
# Replace 0 p-values with a small number to be able to calculate the rank properly
RL_mono <- RL_mono %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue))

RL_mono_rankedlist <- RL_mono %>%
  mutate(rank = -log10(avg_pvalue))%>%
  mutate(rank2 = -log10(avg_pvalue) * avg_log2FC)

# COMMAND ----------

RL_mono_rankedlist$hgnc_symbol <- RL_mono_rankedlist$gene
gene_list_RL_MONO <- left_join(RL_mono_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_MONO <- gene_list_RL_MONO %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_RL_MONO  <- gene_list_entrez_id_RL_MONO$rank 

#entrez id as the names in geneList
names(geneList_entrez_RL_MONO) <- gene_list_entrez_id_RL_MONO$entrezgene_id

# COMMAND ----------

# Clean and sort the geneList
geneList_entrez_RL_MONO <- geneList_entrez_RL_MONO[is.finite(geneList_entrez_RL_MONO)]

geneList_entrez_RL_MONO <- sort(geneList_entrez_RL_MONO, decreasing = TRUE)

geneList_entrez_RL_MONO

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

# Pathway enrichment with KEGG in MONO cell common markers
mono_RL_KEGG <- gseKEGG(geneList_entrez_RL_MONO,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        verbose      = FALSE,
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps=0) #no pvalue aproximations to control the warning

# Extract the result from the KEGG enrichment analysis
res_mono_RL_KEGG <- mono_RL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_mono_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_mono_RL_KEGG)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(mono_RL_KEGG, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - MONO - L vs HC - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#RL Reactome Enrichment for MONO
mono_RL_Reactome <- gsePathway(geneList_entrez_RL_MONO,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500,
                eps = 0)#BECAUSE OF A WARNING OF EXTREMELY LOW P VALUES
res_mono_RL_Reactome <- mono_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_mono_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mono_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(mono_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - MONO - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

set.seed(123)
mono_RL_HALLMARKS <- GSEA(geneList_entrez_RL_MONO, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)
res_mono_RL_HALLMARKS <- mono_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_mono_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mono_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
mono_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_MONO, 
                        TERM2GENE = C2_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500,
                        eps=0)#again warning of low pvalues
res_mono_RL_HALLMARKS_C2 <- mono_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_mono_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mono_RL_HALLMARKS_C2)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(mono_RL_HALLMARKS, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - MONO - L vs HC - HALLMARKS")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(mono_RL_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - MONO - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

mono_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_MONO),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_mono_RL_GO <- mono_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_mono_RL_GO <- res_mono_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_mono_RL_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(mono_RL_GO, showCategory=30, label_format=100, font.size=9) + ggtitle("RL- MONO - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_mono_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/pathway_enrich_mono_RL_KEGG.xlsx")
write.xlsx(res_mono_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/pathway_enrich_mono_RL_REACTOME.xlsx")
write.xlsx(res_mono_RL_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/pathway_enrich_mono_RL_HALLMARKS.xlsx")
write.xlsx(res_mono_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/pathway_enrich_mono_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_mono_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Monocyte/functional_enrich_mono_RL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##DC

# COMMAND ----------

RL_dc <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/RL_dc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
RL_dc <- RL_dc %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(RL_dc)

# COMMAND ----------

#Rank
# Replace 0 p-values with a small number to be able to calculate the rank properly
RL_dc <- RL_dc %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue))

RL_dc_rankedlist <- RL_dc %>%
  mutate(rank = -log10(avg_pvalue))%>%
  mutate(rank2 = -log10(avg_pvalue) * avg_log2FC)

# COMMAND ----------

RL_dc_rankedlist$hgnc_symbol <- RL_dc_rankedlist$gene
gene_list_RL_DC <- left_join(RL_dc_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_DC <- gene_list_RL_DC %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_RL_DC  <- gene_list_entrez_id_RL_DC$rank 

#entrez id as the names in geneList
names(geneList_entrez_RL_DC) <- gene_list_entrez_id_RL_DC$entrezgene_id

# COMMAND ----------

# Clean and sort the geneList
geneList_entrez_RL_DC <- geneList_entrez_RL_DC[is.finite(geneList_entrez_RL_DC)]

geneList_entrez_RL_DC <- sort(geneList_entrez_RL_DC, decreasing = TRUE)

geneList_entrez_RL_DC

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

# Pathway enrichment with KEGG in DC cell common markers
dc_RL_KEGG <- gseKEGG(geneList_entrez_RL_DC,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_dc_RL_KEGG <- dc_RL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_dc_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_dc_RL_KEGG)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(dc_RL_KEGG, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - DC - L vs HC - KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#RL Reactome Enrichment for DC
dc_RL_Reactome <- gsePathway(geneList_entrez_RL_DC,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500,
                eps = 0)#BECAUSE OF A WARNING OF EXTREMELY LOW P VALUES
res_dc_RL_Reactome <- dc_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_dc_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_dc_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(dc_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - DC - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

set.seed(123)
dc_RL_HALLMARKS <- GSEA(geneList_entrez_RL_DC, 
                        TERM2GENE = H_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500)
res_dc_RL_HALLMARKS <- dc_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_dc_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_dc_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
dc_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_DC, 
                        TERM2GENE = C2_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500,
                        eps=0)#again  warning of low pvalues
res_dc_RL_HALLMARKS_C2 <- dc_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_dc_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_dc_RL_HALLMARKS_C2)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(dc_RL_HALLMARKS, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - DC - L vs HC - HALLMARKS")

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(dc_RL_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - DC - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

dc_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_DC),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_dc_RL_GO <- dc_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_dc_RL_GO <- res_dc_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_dc_RL_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(dc_RL_GO, showCategory=50, label_format=100, font.size=9) + ggtitle("RL- DC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_dc_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/pathway_enrich_dc_RL_KEGG.xlsx")
write.xlsx(res_dc_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/pathway_enrich_dc_RL_REACTOME.xlsx")
write.xlsx(res_dc_RL_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/pathway_enrich_dc_RL_HALLMARKS.xlsx")
write.xlsx(res_dc_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/pathway_enrich_dc_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_dc_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/DC/functional_enrich_dc_RL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Treg

# COMMAND ----------

ARL_treg <- read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/ARL_treg_LvsHC_allmarkers.xlsx")
RL_treg <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_treg_LvsHC_allmarkers.xlsx")

# COMMAND ----------

ARL_treg <- ARL_treg %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))
RL_treg <- RL_treg %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(ARL_treg)

# COMMAND ----------

display(RL_treg)

# COMMAND ----------

#Rank by pvalue ARL
ARL_treg_rankedlist <- ARL_treg %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+avg_log2FC) * -log10(avg_pvalue))

# COMMAND ----------

#R-L ranking
RL_treg_rankedlist <- RL_treg %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+ avg_log2FC)  * -log10(avg_pvalue))

# COMMAND ----------

display(RL_treg_rankedlist)

# COMMAND ----------

# Assigning the 'gene' column values to the 'hgnc_symbol'
ARL_treg_rankedlist$hgnc_symbol <- ARL_treg_rankedlist$gene
RL_treg_rankedlist$hgnc_symbol <- RL_treg_rankedlist$gene

# COMMAND ----------

#Merge by symbol and select relevant columns for gene list
#ARL
gene_list_ARL_Treg <- left_join(ARL_treg_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_ARL_Treg <- gene_list_ARL_Treg %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()
#RL
gene_list_RL_Treg <- left_join(RL_treg_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_Treg <- gene_list_RL_Treg %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_ARL_Treg  <- gene_list_entrez_id_ARL_Treg$rank #add the rank to the geneList
geneList_entrez_RL_Treg  <- gene_list_entrez_id_RL_Treg$rank 

#entrez id as the names in geneList
names(geneList_entrez_ARL_Treg) <- gene_list_entrez_id_ARL_Treg$entrezgene_id 
names(geneList_entrez_RL_Treg) <- gene_list_entrez_id_RL_Treg$entrezgene_id 

#The same procedure but with the second type of ranking rank2
geneList_entrez_ARL_Treg_1  <- gene_list_entrez_id_ARL_Treg$rank2 #add the rank to the geneList
geneList_entrez_RL_Treg_1  <- gene_list_entrez_id_RL_Treg$rank2

names(geneList_entrez_ARL_Treg_1) <- gene_list_entrez_id_ARL_Treg$entrezgene_id 
names(geneList_entrez_RL_Treg_1) <- gene_list_entrez_id_RL_Treg$entrezgene_id

# COMMAND ----------

# MAGIC
# MAGIC %md
# MAGIC Clean and sort the gene lists to use them in gse

# COMMAND ----------

# Clean and sort the ARL geneList
geneList_entrez_ARL_Treg <- geneList_entrez_ARL_Treg[is.finite(geneList_entrez_ARL_Treg)]
geneList_entrez_ARL_Treg <- sort(geneList_entrez_ARL_Treg, decreasing = TRUE)

# Clean and sort the ARL geneList with rank2
geneList_entrez_ARL_Treg_1 <- geneList_entrez_ARL_Treg_1[is.finite(geneList_entrez_ARL_Treg_1)]
geneList_entrez_ARL_Treg_1 <- sort(geneList_entrez_ARL_Treg_1, decreasing = TRUE)

geneList_entrez_ARL_Treg
geneList_entrez_ARL_Treg_1

# COMMAND ----------

# Clean and sort the R-L geneList
geneList_entrez_RL_Treg <- geneList_entrez_RL_Treg[is.finite(geneList_entrez_RL_Treg)]
geneList_entrez_RL_Treg <- sort(geneList_entrez_RL_Treg, decreasing = TRUE)

# Clean and sort the R-L geneList with rank2
geneList_entrez_RL_Treg_1 <- geneList_entrez_RL_Treg_1[is.finite(geneList_entrez_RL_Treg_1)]
geneList_entrez_RL_Treg_1 <- sort(geneList_entrez_RL_Treg_1, decreasing = TRUE)

geneList_entrez_RL_Treg
geneList_entrez_RL_Treg_1

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

set.seed(123)
treg_ARL_KEGG <- gseKEGG(geneList_entrez_ARL_Treg,
                         organism     = 'hsa',
                         pvalueCutoff = 0.05,
                         verbose      = FALSE,
                         minGSSize = 10,
                         maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_treg_ARL_KEGG <- treg_ARL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_treg_ARL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_ARL_KEGG)
}

# COMMAND ----------

set.seed(123)

treg_ARL_KEGG_1 <- gseKEGG(geneList_entrez_ARL_Treg_1,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05,
                          verbose      = FALSE,
                          minGSSize = 10,
                          maxGSSize = 500)                       
                          #,scoreType = "pos") # because I only have upregulated genes

# Extract the result from the KEGG enrichment analysis
res_treg_ARL_KEGG_1 <- treg_ARL_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_treg_ARL_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_ARL_KEGG_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results in ARL common markers

# COMMAND ----------

set.seed(123)

#Pathway enrichment with KEGG in R-L common markers
treg_RL_KEGG <- gseKEGG(geneList_entrez_RL_Treg,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)


res_treg_RL_KEGG <- treg_RL_KEGG@result
# Check if there are any enriched pathways
if (nrow(res_treg_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_RL_KEGG)
}

# COMMAND ----------

set.seed(123)
#Pathway enrichment with KEGG in R-L common markers using RANK2
treg_RL_KEGG_1 <- gseKEGG(geneList_entrez_RL_Treg_1,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)

res_treg_RL_KEGG_1 <- treg_RL_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_treg_RL_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_RL_KEGG_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results in RL treg common markers pathway enrichment using rank2

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
dotplot(treg_RL_KEGG, x = "NES", showCategory= 30, size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Treg - L vs HC- KEGG")

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#ARL Reactome Enrichment
treg_ARL_Reactome <- gsePathway(geneList_entrez_ARL_Treg,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_treg_ARL_Reactome <- treg_ARL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_treg_ARL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_treg_ARL_Reactome)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant pathway enrichment in ARL Treg

# COMMAND ----------

set.seed(123)
#R-L Reactome enrichment
treg_RL_Reactome <- gsePathway(geneList_entrez_RL_Treg,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_treg_RL_Reactome <- treg_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_treg_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_treg_RL_Reactome)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(treg_RL_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Treg - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

# MAGIC %md
# MAGIC A-R-L HALLMARKS enrichment

# COMMAND ----------

set.seed(123)
treg_ARL_HALLMARKS <- GSEA(geneList_entrez_ARL_Treg, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize = 10,
                          maxGSSize=500)

res_treg_ARL_HALLMARKS <- treg_ARL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_treg_ARL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_treg_ARL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
treg_ARL_HALLMARKS_C2 <- GSEA(geneList_entrez_ARL_Treg, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.2, 
                          minGSSize=10,
                          maxGSSize=500)

res_treg_ARL_HALLMARKS_C2 <- treg_ARL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_treg_ARL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the  enrichment analysis result 
  display(res_treg_ARL_HALLMARKS_C2)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results in ARL Treg

# COMMAND ----------

# MAGIC %md
# MAGIC R-L hallmarks enrichment

# COMMAND ----------

set.seed(123)
treg_RL_HALLMARKS <- GSEA(geneList_entrez_RL_Treg, 
                          TERM2GENE = H_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500)
res_treg_RL_HALLMARKS <- treg_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_treg_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)

#HALLMARKS C2
treg_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_Treg, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                          minGSSize=10,
                          maxGSSize=500,
                          eps=0)#To control low pvalue warning, it does not use approximations
res_treg_RL_HALLMARKS_C2 <- treg_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_treg_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_treg_RL_HALLMARKS_C2)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(treg_RL_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("RL - Treg - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

treg_ARL_GO <- enrichGO(gene = names(geneList_entrez_ARL_Treg),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_treg_ARL_GO <- treg_ARL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_treg_ARL_GO <- res_treg_ARL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_treg_ARL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(treg_ARL_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("ARL- Treg - L vs HC - GO")

# COMMAND ----------

set.seed(123)

treg_ARL_GO_1 <- enrichGO(gene = names(geneList_entrez_ARL_Treg_1),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_treg_ARL_GO_1 <- treg_ARL_GO_1@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_treg_ARL_GO_1 <- res_treg_ARL_GO_1 %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_treg_ARL_GO_1)

# COMMAND ----------

# MAGIC %md
# MAGIC Same results if using rank2

# COMMAND ----------

treg_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_Treg),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)

res_treg_RL_GO <- treg_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_treg_RL_GO <- res_treg_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(res_treg_RL_GO)

# COMMAND ----------

options(repr.plot.width=1000, repr.plot.height=1000,  echo= F)
barplot(treg_RL_GO, showCategory=30, label_format=50, font.size=9) + ggtitle("RL - Treg - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results RL
write.xlsx(res_treg_RL_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_PE/pathway_enrich_treg_RL_KEGG.xlsx")
write.xlsx(res_treg_RL_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_PE/pathway_enrich_treg_RL_REACTOME.xlsx")
write.xlsx(res_treg_RL_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_PE/pathway_enrich_treg_RL_HALLMARKS.xlsx")
write.xlsx(res_treg_RL_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_PE/pathway_enrich_treg_RL_HALLMARKS_C2.xlsx")
write.xlsx(res_treg_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/RL_PE/functional_enrich_treg_RL_GO.xlsx")

#Save results ARL
write.xlsx(res_treg_ARL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Treg/ARL_PE/functional_enrich_treg_ARL_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Fibroblasts

# COMMAND ----------

AR_fb <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/AR_fb_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
AR_fb <- AR_fb %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(AR_fb)

# COMMAND ----------

#Rank by pvalue
AR_fb_rankedlist <- AR_fb %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+avg_log2FC) * -log10(avg_pvalue))

# COMMAND ----------

display(AR_fb_rankedlist)

# COMMAND ----------

#Add new column with the gene names as hgnc symbol
AR_fb_rankedlist$hgnc_symbol <- AR_fb_rankedlist$gene 

# COMMAND ----------

gene_list_AR_fb <- left_join(AR_fb_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_AR_fb <- gene_list_AR_fb %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_AR_fb  <- gene_list_entrez_id_AR_fb$rank #add the rank to the geneList
#entrez id as the names in geneList
names(geneList_entrez_AR_fb) <- gene_list_entrez_id_AR_fb$entrezgene_id 

#The same procedure but with the second type of ranking rank2
geneList_entrez_AR_fb_1  <- gene_list_entrez_id_AR_fb$rank2 #add the rank to the geneList
names(geneList_entrez_AR_fb_1) <- gene_list_entrez_id_AR_fb$entrezgene_id

# COMMAND ----------

# Clean and sort the AR geneList
geneList_entrez_AR_fb <- geneList_entrez_AR_fb[is.finite(geneList_entrez_AR_fb)]
geneList_entrez_AR_fb <- sort(geneList_entrez_AR_fb, decreasing = TRUE)

# Clean and sort the AR geneList with rank2
geneList_entrez_AR_fb_1 <- geneList_entrez_AR_fb_1[is.finite(geneList_entrez_AR_fb_1)]
geneList_entrez_AR_fb_1 <- sort(geneList_entrez_AR_fb_1, decreasing = TRUE)

geneList_entrez_AR_fb
geneList_entrez_AR_fb_1

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

set.seed(123)
fb_AR_KEGG <- gseKEGG(geneList_entrez_AR_fb,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_fb_AR_KEGG <- fb_AR_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_fb_AR_KEGG)
}

# COMMAND ----------

set.seed(123)
fb_AR_KEGG_1 <- gseKEGG(geneList_entrez_AR_fb_1,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        verbose      = FALSE,
                        minGSSize = 10,
                        maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_fb_AR_KEGG_1 <- fb_AR_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_fb_AR_KEGG_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results with KEGG pathways

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#AR Reactome Enrichment
fb_AR_Reactome <- gsePathway(geneList_entrez_AR_fb,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_fb_AR_Reactome <- fb_AR_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_fb_AR_Reactome)
}

# COMMAND ----------

set.seed(123)
#AR Reactome Enrichment with rank 2
fb_AR_Reactome_1 <- gsePathway(geneList_entrez_AR_fb_1,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_fb_AR_Reactome_1 <- fb_AR_Reactome_1@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_Reactome_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_fb_AR_Reactome_1)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(fb_AR_Reactome, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("AR - FB - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

set.seed(123)
fb_AR_HALLMARKS <- GSEA(geneList_entrez_AR_fb, 
                        TERM2GENE = H_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500)
res_fb_AR_HALLMARKS <- fb_AR_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_fb_AR_HALLMARKS)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant pathway enrichment

# COMMAND ----------

set.seed(123)

# HALLMARKS C2
fb_AR_HALLMARKS_C2 <- GSEA(geneList_entrez_AR_fb,
                          TERM2GENE = C2_t2g,
                          pvalueCutoff = 0.05,
                          minGSSize = 5,
                          maxGSSize = 500) #No results with minimum 10
res_fb_AR_HALLMARKS_C2_5minimum <- fb_AR_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_fb_AR_HALLMARKS_C2_5minimum) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_fb_AR_HALLMARKS_C2_5minimum)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(fb_AR_HALLMARKS_C2, showCategory= 30, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("AR - KC - L vs HC - hallmarks c2 (5 minimum)")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

fb_AR_GO <- enrichGO(gene = names(geneList_entrez_AR_fb),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_fb_AR_GO <- fb_AR_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_fb_AR_GO <- res_fb_AR_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_fb_AR_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(fb_AR_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("AR- FB - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_fb_AR_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/pathway_enrich_fb_AR_KEGG.xlsx")
write.xlsx(res_fb_AR_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/pathway_enrich_fb_AR_REACTOME.xlsx")
write.xlsx(res_fb_AR_HALLMARKS_C2_5minimum, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/pathway_enrich_fb_AR_HALLMARKS_C2_min5.xlsx")
write.xlsx(res_fb_AR_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Fibroblast/functional_enrich_fb_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##Keratinocytes

# COMMAND ----------

AR_kc <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/AR_kc_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
AR_kc <- AR_kc %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(AR_kc)

# COMMAND ----------

#Rank by pvalue
AR_kc_rankedlist <- AR_kc %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue)) %>%
  mutate(rank = -log10(avg_pvalue)*sign(avg_log2FC)) %>%
  mutate(rank2 = (1+avg_log2FC) * -log10(avg_pvalue))

# COMMAND ----------

display(AR_kc_rankedlist)

# COMMAND ----------

#Add new column with the gene names as hgnc symbol
AR_kc_rankedlist$hgnc_symbol <- AR_kc_rankedlist$gene 

# COMMAND ----------

gene_list_AR_kc <- left_join(AR_kc_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_AR_kc <- gene_list_AR_kc %>% dplyr::select(entrezgene_id, rank, rank2) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_AR_kc  <- gene_list_entrez_id_AR_kc$rank #add the rank to the geneList
#entrez id as the names in geneList
names(geneList_entrez_AR_kc) <- gene_list_entrez_id_AR_kc$entrezgene_id 

#The same procedure but with the second type of ranking rank2
geneList_entrez_AR_kc_1  <- gene_list_entrez_id_AR_kc$rank2 #add the rank to the geneList
names(geneList_entrez_AR_kc_1) <- gene_list_entrez_id_AR_kc$entrezgene_id 


# COMMAND ----------

# Clean and sort the AR geneList
geneList_entrez_AR_kc <- geneList_entrez_AR_kc[is.finite(geneList_entrez_AR_kc)]
geneList_entrez_AR_kc <- sort(geneList_entrez_AR_kc, decreasing = TRUE)

# Clean and sort the AR geneList with rank2
geneList_entrez_AR_kc_1 <- geneList_entrez_AR_kc_1[is.finite(geneList_entrez_AR_kc_1)]
geneList_entrez_AR_kc_1 <- sort(geneList_entrez_AR_kc_1, decreasing = TRUE)

geneList_entrez_AR_kc
geneList_entrez_AR_kc_1

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

set.seed(123)
kc_AR_KEGG <- gseKEGG(geneList_entrez_AR_kc,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose      = FALSE,
                      minGSSize = 10,
                      maxGSSize = 500)
                     #,scoreType = "pos") # because I only have upregulated genes

# Extract the result from the KEGG enrichment analysis
res_kc_AR_KEGG <- kc_AR_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_kc_AR_KEGG)
}

# COMMAND ----------

set.seed(123)
kc_AR_KEGG_1 <- gseKEGG(geneList_entrez_AR_kc_1,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        verbose      = FALSE,
                        minGSSize = 10,
                        maxGSSize = 500)
                       #,scoreType = "pos") # because I only have upregulated genes

# Extract the result from the KEGG enrichment analysis
res_kc_AR_KEGG_1 <- kc_AR_KEGG_1@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_KEGG_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_kc_AR_KEGG_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant results with rank2

# COMMAND ----------

# MAGIC %md
# MAGIC ###REACTOME

# COMMAND ----------

set.seed(123)
#AR Reactome Enrichment
kc_AR_Reactome <- gsePathway(geneList_entrez_AR_kc,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_kc_AR_Reactome <- kc_AR_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_kc_AR_Reactome)
}

# COMMAND ----------

set.seed(123)
#AR Reactome Enrichment
kc_AR_Reactome_1 <- gsePathway(geneList_entrez_AR_kc_1,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)
res_kc_AR_Reactome_1 <- kc_AR_Reactome_1@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_Reactome_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_kc_AR_Reactome_1)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No significant enrichment with rank2

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1500,  echo= F)
dotplot(kc_AR_Reactome, showCategory= 50, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("AR - KC - L vs HC - Reactome")

# COMMAND ----------

# MAGIC %md
# MAGIC ###HALLMARKS

# COMMAND ----------

set.seed(123)
kc_AR_HALLMARKS <- GSEA(geneList_entrez_AR_kc, 
                        TERM2GENE = H_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=10,
                        maxGSSize=500)
res_kc_AR_HALLMARKS <- kc_AR_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_kc_AR_HALLMARKS)
}

# COMMAND ----------

set.seed(123)

#HALLMARKS C2
kc_AR_HALLMARKS_C2 <- GSEA(geneList_entrez_AR_kc, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                            minGSSize=10,
                          maxGSSize=500)
res_kc_AR_HALLMARKS_C2 <- kc_AR_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_kc_AR_HALLMARKS_C2)
}

# COMMAND ----------

set.seed(123)

#HALLMARKS C2 but with type ranking rank2
kc_AR_HALLMARKS_C2_1 <- GSEA(geneList_entrez_AR_kc_1, 
                          TERM2GENE = C2_t2g, 
                          pvalueCutoff = 0.05, 
                            minGSSize=10,
                          maxGSSize=500)
res_kc_AR_HALLMARKS_C2_1 <- kc_AR_HALLMARKS_C2_1@result

# Check if there are any enriched pathways
if (nrow(res_kc_AR_HALLMARKS_C2_1) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_kc_AR_HALLMARKS_C2_1)
}

# COMMAND ----------

options(repr.plot.width=1500, repr.plot.height=1500,  echo= F)
dotplot(kc_AR_HALLMARKS_C2, showCategory= 25, x = "NES", size = "GeneRatio", font.size = 9, label_format = 50) + ggtitle("AR - KC - L vs HC - HALLMARKS C2")

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

kc_AR_GO <- enrichGO(gene = names(geneList_entrez_AR_kc),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_kc_AR_GO <- kc_AR_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_kc_AR_GO <- res_kc_AR_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_kc_AR_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(kc_AR_GO, showCategory=25, label_format=100, font.size=9) + ggtitle("AR- KC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

#Save results
write.xlsx(res_kc_AR_KEGG, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/pathway_enrich_kc_AR_KEGG.xlsx")
write.xlsx(res_kc_AR_Reactome, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/pathway_enrich_kc_AR_REACTOME.xlsx")
write.xlsx(res_kc_AR_HALLMARKS, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/pathway_enrich_kc_AR_HALLMARKS.xlsx")
write.xlsx(res_kc_AR_HALLMARKS_C2, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/pathway_enrich_kc_AR_HALLMARKS_C2.xlsx")
write.xlsx(res_kc_AR_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/Keratinocytes/functional_enrich_kc_AR_GO.xlsx")

# COMMAND ----------

# MAGIC %md
# MAGIC ##MastC

# COMMAND ----------

RL_mast <-read.xlsx("/dbfs/mnt/sandbox/TFM_PAULA/common_markers/MastC/RL_mast_LvsHC_allmarkers.xlsx")

# COMMAND ----------

#Make negative healthy avglog2FC to be able to see which are from healthy and which from lesional in the  pathway plots
RL_mast <- RL_mast %>%
  mutate(avg_log2FC = ifelse(Condition == "healthy", -abs(avg_log2FC), avg_log2FC))

# COMMAND ----------

display(RL_mast)

# COMMAND ----------

# MAGIC %md
# MAGIC Few DEGs

# COMMAND ----------

#Rank
# Replace 0 p-values with a small number to be able to calculate the rank properly
RL_mast <- RL_mast %>%
  mutate(avg_pvalue = ifelse(avg_pvalue == 0, 1e-300, avg_pvalue))

RL_mast_rankedlist <- RL_mast %>%
  mutate(rank = -log10(avg_pvalue))%>%
  mutate(rank2 = -log10(avg_pvalue) * avg_log2FC)

# COMMAND ----------

RL_mast_rankedlist$hgnc_symbol <- RL_mast_rankedlist$gene
gene_list_RL_MAST <- left_join(RL_mast_rankedlist, G_list, by = "hgnc_symbol") %>% distinct_all()
gene_list_entrez_id_RL_MAST <- gene_list_RL_MAST %>% dplyr::select(entrezgene_id, rank) %>% distinct(entrezgene_id, .keep_all = TRUE) %>% drop_na()

# COMMAND ----------

#Prepare the gene List
geneList_entrez_RL_MAST  <- gene_list_entrez_id_RL_MAST$rank 

#entrez id as the names in geneList
names(geneList_entrez_RL_MAST) <- gene_list_entrez_id_RL_MAST$entrezgene_id

# COMMAND ----------

# Clean and sort the geneList
geneList_entrez_RL_MAST <- geneList_entrez_RL_MAST[is.finite(geneList_entrez_RL_MAST)]

geneList_entrez_RL_MAST <- sort(geneList_entrez_RL_MAST, decreasing = TRUE)

geneList_entrez_RL_MAST

# COMMAND ----------

# MAGIC %md
# MAGIC ###KEGG ENRICHMENT

# COMMAND ----------

# Pathway enrichment with KEGG in MAST cell common markers
mast_RL_KEGG <- gseKEGG(geneList_entrez_RL_MAST,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05,
                        verbose      = FALSE,
                        minGSSize = 10,
                        maxGSSize = 500)

# Extract the result from the KEGG enrichment analysis
res_mast_RL_KEGG <- mast_RL_KEGG@result

# Check if there are any enriched pathways
if (nrow(res_mast_RL_KEGG) == 0) {
  print("No pathways enriched")
} else {
  # Display the KEGG enrichment analysis result 
  display(res_mast_RL_KEGG)
}

# COMMAND ----------

set.seed(123)
#RL Reactome Enrichment for MAST
mast_RL_Reactome <- gsePathway(geneList_entrez_RL_MAST,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE,
                organism = "human",
                minGSSize=10,
                maxGSSize=500)

res_mast_RL_Reactome <- mast_RL_Reactome@result

# Check if there are any enriched pathways
if (nrow(res_mast_RL_Reactome) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mast_RL_Reactome)
}

# COMMAND ----------

set.seed(123)
mast_RL_HALLMARKS <- GSEA(geneList_entrez_RL_MAST, 
                        TERM2GENE = H_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=1,
                        maxGSSize=500)
res_mast_RL_HALLMARKS <- mast_RL_HALLMARKS@result

# Check if there are any enriched pathways
if (nrow(res_mast_RL_HALLMARKS) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mast_RL_HALLMARKS)
}

# COMMAND ----------

set.seed(123)
mast_RL_HALLMARKS_C2 <- GSEA(geneList_entrez_RL_MAST, 
                        TERM2GENE = C2_t2g, 
                        pvalueCutoff = 0.05, 
                        minGSSize=1,
                        maxGSSize=500,
                        eps=0)#again warning of low pvalues
res_mast_RL_HALLMARKS_C2 <- mast_RL_HALLMARKS_C2@result

# Check if there are any enriched pathways
if (nrow(res_mast_RL_HALLMARKS_C2) == 0) {
  print("No pathways enriched")
} else {
  # Display the enrichment analysis result 
  display(res_mast_RL_HALLMARKS_C2)
}

# COMMAND ----------

# MAGIC %md
# MAGIC No pathways enriched with that few genes.

# COMMAND ----------

# MAGIC %md
# MAGIC ###GO

# COMMAND ----------

set.seed(123)

mast_RL_GO <- enrichGO(gene = names(geneList_entrez_RL_MAST),
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.2)
                
# Extract results
res_mast_RL_GO <- mast_RL_GO@result

# Filter results to include only those with p.adjust < 0.05
filtered_res_mast_RL_GO <- res_mast_RL_GO %>%
  filter(p.adjust < 0.05)

# Display filtered results
display(filtered_res_mast_RL_GO)

# COMMAND ----------

options(repr.plot.width=1200, repr.plot.height=1200,  echo= F)
barplot(mast_RL_GO, showCategory=50, label_format=100, font.size=9) + ggtitle("RL- MastC - L vs HC - GO")

# COMMAND ----------

# MAGIC %md
# MAGIC ###Save results

# COMMAND ----------

write.xlsx(res_mast_RL_GO, "/dbfs/mnt/sandbox/TFM_PAULA/common_markers/MastC/functional_enrich_mast_RL_GO.xlsx")
