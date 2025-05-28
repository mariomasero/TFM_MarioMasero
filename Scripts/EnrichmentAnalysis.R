## Load required libraries
library(readr)

## Read proteomics results file
ProteomicsResults <- read_delim("ResultadosProteomica.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

## GO ENRICHMENT ANALYSIS
library(dplyr)
library(tidyr)

# Extract Gene ID vector
gene_ids <- ProteomicsResults$`Gene ID`

# Create raw expression matrices
expr_A_raw <- as.matrix(ProteomicsResults[ ,(5:14)])
expr_H_raw <- as.matrix(ProteomicsResults[ ,(15:24)])

# Rename columns
colnames(expr_A_raw) <- c("C_1", "C_2", "C_3", "C_4", "C_5",
                          "T_1", "T_2", "T_3", "T_4", "T_5")
colnames(expr_H_raw) <- c("C_1", "C_2", "C_3", "C_4", "C_5",
                          "T_1", "T_2", "T_3", "T_4", "T_5")

# Convert matrices to data frames and add Gene ID column
expr_A_raw <- as.data.frame(expr_A_raw) %>%
  mutate(`Gene ID` = gene_ids) %>%
  relocate(`Gene ID`, .before = C_1)

expr_H_raw <- as.data.frame(expr_H_raw) %>%
  mutate(`Gene ID` = gene_ids) %>%
  relocate(`Gene ID`, .before = C_1)

# Remove rows where all LFQ values are NA
expr_A_filtered <- expr_A_raw[rowSums(!is.na(expr_A_raw[, -1])) > 0, ]
expr_H_filtered <- expr_H_raw[rowSums(!is.na(expr_H_raw[, -1])) > 0, ]

# Split multi-gene entries separated by ";" into separate rows
expr_A_expanded <- expr_A_filtered %>%
  separate_rows(`Gene ID`, sep = ";")

expr_H_expanded <- expr_H_filtered %>%
  separate_rows(`Gene ID`, sep = ";")

## GO Enrichment

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Prepare gene symbol vectors
gene_list_A <- expr_A_expanded$`Gene ID`
gene_list_H <- expr_H_expanded$`Gene ID`

# Run enrichGO for each list
ego_A <- enrichGO(
  gene          = gene_list_A,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "SYMBOL"
)

ego_H <- enrichGO(
  gene          = gene_list_H,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "SYMBOL"
)

## Plots
barplot(ego_A, showCategory = 20)
dotplot(ego_A, showCategory = 20)
treeplot(pairwise_termsim(ego_A), showCategory = 20)
emapplot(pairwise_termsim(ego_A), showCategory = 20)
cnetplot(ego_A, showCategory = 20)

barplot(ego_H, showCategory = 20)
dotplot(ego_H, showCategory = 20)
treeplot(pairwise_termsim(ego_H), showCategory = 20)
emapplot(pairwise_termsim(ego_H), showCategory = 20)
cnetplot(ego_H, showCategory = 20)

###############################################################################

## GSEA ANALYSIS

library(ggplot2)
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

# Sort gene list (required for GSEA)
proteinListA <- dep_A$logFC
names(proteinListA) <- dep_A$ID
proteinListA <- sort(proteinListA, decreasing = TRUE)

proteinListH <- dep_H$logFC
names(proteinListH) <- dep_H$ID
proteinListH <- sort(proteinListH, decreasing = TRUE)

gse_A <- gseGO(
  geneList       = proteinListA,
  ont            = "ALL",
  keyType        = "SYMBOL",
  OrgDb          = org.Hs.eg.db,
  minGSSize      = 3,
  maxGSSize      = 800,
  pvalueCutoff   = 0.05,
  pAdjustMethod  = "fdr",
  verbose        = TRUE
)

gse_H <- gseGO(
  geneList       = proteinListH,
  ont            = "ALL",
  keyType        = "SYMBOL",
  OrgDb          = org.Hs.eg.db,
  minGSSize      = 3,
  maxGSSize      = 800,
  pvalueCutoff   = 0.05,
  pAdjustMethod  = "fdr",
  verbose        = TRUE
)

## Visualization
gse_CC <- gse_A
gse_CC@result <- gse_CC@result[gse_CC@result$ONTOLOGY == "CC", ]

require(DOSE)  # Required for semantic similarity and GSEA plotting
dotplot(gse_CC, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign)

ema_gse_A <- pairwise_termsim(gse_A)
emapplot(ema_gse_A, showCategory = 10)

# Literature trends
terms <- c(
  "plasma membrane",
  "cell surface",
  "transmembrane protein",
  "cell membrane protein",
  "membrane receptor",
  "cell adhesion molecule",
  "membrane transporter"
)
pmcplot(terms, period = 2010:2023, proportion = TRUE)

# Filter membrane-related GO terms
term_pattern <- "membrane|cell surface|plasma membrane|transmembrane|junction|immune response-activating cell surface receptor signaling pathway"

membrane_terms_A <- grep(term_pattern, gse_A@result$Description, 
                         value = TRUE, ignore.case = TRUE)

membrane_terms_H <- grep(term_pattern, gse_H@result$Description, 
                         value = TRUE, ignore.case = TRUE)

# Subset GSEA objects
gse_membraneA <- gse_A
gse_membraneA@result <- gse_A@result[gse_A@result$Description %in% membrane_terms_A, ]

cnetplot(gse_membraneA,
         showCategory = length(membrane_terms_A),
         foldChange = proteinListA,
         layout = "kk",
         circular = FALSE,
         node_label = "all",
         cex_label_category = 0.9,
         cex_label_gene = 0.5
) +
  ggtitle("A549 Cell Line - Gene-Term Network: Membrane & Surface-Related GO Terms")

gse_membraneH <- gse_H
gse_membraneH@result <- gse_H@result[gse_H@result$Description %in% membrane_terms_H, ]

###############################################################################

## compareCluster – GSEA comparison

df_gsea <- rbind(
  data.frame(SYMBOL = dep_A$ID, logFC = dep_A$logFC, group = "A549"),
  data.frame(SYMBOL = dep_H$ID, logFC = dep_H$logFC, group = "H2009")
)

gsea_comparison <- compareCluster(
  SYMBOL | logFC ~ group,
  data = df_gsea,
  fun = "gseGO",
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05
)

png("Images/Gsea_comparison.png", width = 6500, height = 5000, res = 600)
dotplot(gsea_comparison, showCategory = 20, split = ".sign", font.size = 8) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA Comparison – A549 vs H2009")
dev.off()

## compareCluster – enrichGO comparison

df_go <- list(
  A549 = protein_ids_A,
  H2009 = protein_ids_H
)

go_comparison <- compareCluster(
  geneCluster = df_go,
  fun = "enrichGO",
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pvalueCutoff = 0.05
)

png("Images/FigureS4.png", width = 6500, height = 5000, res = 600)
dotplot(go_comparison, showCategory = 20, font.size = 8) +
  ggtitle("GOCC Comparison – A549 vs H2009")
dev.off()
