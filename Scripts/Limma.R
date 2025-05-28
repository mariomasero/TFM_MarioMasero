## Load required packages
library(readr)

## Read proteomics results file
ProteomicsResults <- read_delim("ResultadosProteomica.csv", 
                                delim = ";", escape_double = FALSE, 
                                trim_ws = TRUE)

## Extract protein names
protein_ids <- ProteomicsResults$`Gene ID`

## Create expression matrix for A cell line
expr_A <- as.matrix(ProteomicsResults[ ,(5:14)])  # Extract LFQ data
rownames(expr_A) <- protein_ids  # Assign protein names

## Create expression matrix for H cell line
expr_H <- as.matrix(ProteomicsResults[ ,(15:24)])  # Extract LFQ data
rownames(expr_H) <- protein_ids  # Assign protein names

## Rename columns
colnames(expr_A) <- c("C_1","C_2","C_3", "C_4", "C_5",
                      "T_1","T_2","T_3", "T_4", "T_5")

colnames(expr_H) <- c("C_1","C_2","C_3", "C_4", "C_5",
                      "T_1","T_2","T_3", "T_4", "T_5")

##### DATA PREPROCESSING BEFORE DIFFERENTIAL EXPRESSION ANALYSIS #####

## Remove proteins with zero variance
# Cell line A
variances_A <- apply(expr_A, 1, var, na.rm = TRUE)
expr_A <- expr_A[!is.na(variances_A) & variances_A != 0, ]

# Cell line H
variances_H <- apply(expr_H, 1, var, na.rm = TRUE)
expr_H <- expr_H[!is.na(variances_H) & variances_H != 0, ]

## Remove proteins with more than 70% missing values
## Helper function to calculate % of NAs per row
na_fraction <- function(x) mean(is.na(x))

# Cell line A
na_frac_A <- apply(expr_A, 1, na_fraction)
expr_A <- expr_A[na_frac_A <= 0.7, ]

# Cell line H
na_frac_H <- apply(expr_H, 1, na_fraction)
expr_H <- expr_H[na_frac_H <= 0.7, ]

## NA imputation using LoD (1/5 of the minimum positive value per row)
impute_LoD <- function(matrix) {
  t(apply(matrix, 1, function(row) {
    if (all(is.na(row))) {
      return(rep(NA, length(row)))
    }
    min_val <- min(row[row > 0], na.rm = TRUE)
    row[is.na(row)] <- min_val / 5
    return(row)
  }))
}

## Apply imputation to both cell lines
expr_A <- impute_LoD(expr_A)
expr_H <- impute_LoD(expr_H)

## ADDITIONAL VARIANCE FILTER: IQR (BOTTOM 30%)

# Function to compute IQR for each row
row_iqr <- function(matrix) {
  apply(matrix, 1, function(x) IQR(x, na.rm = TRUE))
}

# Cell line A
iqr_A <- row_iqr(expr_A)
threshold_A <- quantile(iqr_A, probs = 0.30, na.rm = TRUE)  # 30th percentile
expr_A <- expr_A[iqr_A > threshold_A, ]

# Cell line H
iqr_H <- row_iqr(expr_H)
threshold_H <- quantile(iqr_H, probs = 0.30, na.rm = TRUE)
expr_H <- expr_H[iqr_H > threshold_H, ]

#------------------------------------------------------------------------------
## Apply NormalyzerDE
# Load library
# library(NormalyzerDE)

# Save LFQ data in a suitable format for NormalyzerDE
# write.table(expr_A, file = "LFQ_A_quant.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Run NormalyzerDE
# normalyzer(jobName = "A_LFQ", 
#           designPath = "normalyzer_design.tsv",
#           dataPath = "LFQ_A_quant.tsv",
#           outputDir = "Normalyzer_output/",
#           noLogTransform = TRUE) 
#------------------------------------------------------------------------------

## Are the data normalized?
## Boxplots to visually inspect raw values
ylim_raw <- range(log2(c(expr_A, expr_H) + 1), na.rm = TRUE)

png("Images/boxplotraw.png", width = 4800, height = 2400, res = 600)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

boxplot(log2(expr_A + 1), 
        main = "A549 - Log2(LFQ + 1) (Raw)", 
        ylab = "Protein Expression [log2(LFQ + 1)]", 
        col = rep(c("steelblue", "firebrick"), each = 5),
        outline = FALSE, las = 2,
        cex.main = 0.85, cex.lab = 0.85, cex.axis = 0.7,
        ylim = ylim_raw)

boxplot(log2(expr_H + 1), 
        main = "H2009 - Log2(LFQ + 1) (Raw)", 
        ylab = "Protein Expression [log2(LFQ + 1)]", 
        col = rep(c("steelblue", "firebrick"), each = 5),
        outline = FALSE, las = 2,
        cex.main = 0.85, cex.lab = 0.85, cex.axis = 0.7,
        ylim = ylim_raw)

par(mfrow = c(1, 1))
dev.off()

## Apply quantile normalization

library(preprocessCore)  # provides quantile normalization function

# Normalize expression matrices
expr_A_qn <- normalize.quantiles(as.matrix(expr_A))
expr_H_qn <- normalize.quantiles(as.matrix(expr_H))

# Restore row and column names
rownames(expr_A_qn) <- rownames(expr_A)
colnames(expr_A_qn) <- colnames(expr_A)

rownames(expr_H_qn) <- rownames(expr_H)
colnames(expr_H_qn) <- colnames(expr_H)

## Boxplots to verify normalization
ylim_qn <- range(log2(c(expr_A_qn, expr_H_qn) + 1), na.rm = TRUE)

png("Images/boxplotqn.png", width = 4800, height = 2400, res = 600)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

boxplot(log2(expr_A_qn + 1), 
        main = "A549 - Log2(LFQ + 1) (Quantile Normalized)", 
        ylab = "Protein Expression [log2(LFQ + 1)]", 
        col = rep(c("steelblue", "firebrick"), each = 5),
        outline = FALSE, las = 2,
        cex.main = 0.85, cex.lab = 0.85, cex.axis = 0.7,
        ylim = ylim_qn)

boxplot(log2(expr_H_qn + 1), 
        main = "H2009 - Log2(LFQ + 1) (Quantile Normalized)", 
        ylab = "Protein Expression [log2(LFQ + 1)]", 
        col = rep(c("steelblue", "firebrick"), each = 5),
        outline = FALSE, las = 2,
        cex.main = 0.85, cex.lab = 0.85, cex.axis = 0.7,
        ylim = ylim_qn)

par(mfrow = c(1, 1))
dev.off()

## Apply log2 transformation
expr_A_qn <- log2(expr_A_qn + 1)
expr_H_qn <- log2(expr_H_qn + 1)

## ---------------- DIFFERENTIAL EXPRESSION ANALYSIS WITH LIMMA ---------------

library(limma)

# Experimental design (no intercept: -1)
experimental.design <- model.matrix(~ -1 + factor(c(rep("Control", 5), rep("Treatment", 5))))
colnames(experimental.design) <- c("Control", "Treatment")

# Fit linear model
fit_A <- lmFit(expr_A_qn, experimental.design)
fit_H <- lmFit(expr_H_qn, experimental.design)

# Define contrast: Treatment - Control
contrast.matrix <- makeContrasts(Treatment - Control, levels = experimental.design)

# Apply contrast fit
contrast.fit_A <- contrasts.fit(fit_A, contrast.matrix)
results_A <- eBayes(contrast.fit_A)

contrast.fit_H <- contrasts.fit(fit_H, contrast.matrix)
results_H <- eBayes(contrast.fit_H)

# Extract results for all proteins
dep_A <- topTable(results_A, 
                  number = nrow(expr_A_qn), 
                  coef = 1, 
                  sort.by = "logFC")

dep_H <- topTable(results_H, 
                  number = nrow(expr_H_qn), 
                  coef = 1, 
                  sort.by = "logFC")

# Extract key statistics
log2fc_A <- dep_A$logFC
pval_A <- dep_A$P.Value
protein_ids_A <- dep_A$ID
names(log2fc_A) <- protein_ids_A
names(pval_A) <- protein_ids_A

log2fc_H <- dep_H$logFC
pval_H <- dep_H$P.Value
protein_ids_H <- dep_H$ID
names(log2fc_H) <- protein_ids_H
names(pval_H) <- protein_ids_H

# Define up- and down-regulated proteins (|log2FC| > 1 and p < 0.05)
activated_protein_A <- protein_ids_A[log2fc_A > 1 & pval_A < 0.05]
repressed_protein_A <- protein_ids_A[log2fc_A < -1 & pval_A < 0.05]

length(activated_protein_A)
length(repressed_protein_A)

activated_protein_H <- protein_ids_H[log2fc_H > 1 & pval_H < 0.05]
repressed_protein_H <- protein_ids_H[log2fc_H < -1 & pval_H < 0.05]

length(activated_protein_H)
length(repressed_protein_H)

## ---------------- CREATE TABLE WITH DEP RESULTS ----------------

# Add FDR-adjusted p-values
adj.pval_A <- dep_A$adj.P.Val
names(adj.pval_A) <- rownames(dep_A)

adj.pval_H <- dep_H$adj.P.Val
names(adj.pval_H) <- rownames(dep_H)

# Label DEPs
dep.labels_A <- ifelse(log2fc_A > 1 & pval_A < 0.05, "Upregulated in Treatment",
                       ifelse(log2fc_A < -1 & pval_A < 0.05, "Downregulated in Treatment", 
                              "Not significant"))

dep.labels_H <- ifelse(log2fc_H > 1 & pval_H < 0.05, "Upregulated in Treatment",
                       ifelse(log2fc_H < -1 & pval_H < 0.05, "Downregulated in Treatment", 
                              "Not significant"))

# Create annotated DEP tables
dep.annotated_A <- data.frame(
  Protein = protein_ids_A,
  log2FC = log2fc_A,
  P.Value = pval_A,
  Adj.P.Val = adj.pval_A,
  DEP.Status = dep.labels_A
)

dep.annotated_H <- data.frame(
  Protein = protein_ids_H,
  log2FC = log2fc_H,
  P.Value = pval_H,
  Adj.P.Val = adj.pval_H,
  DEP.Status = dep.labels_H
)

# Filter to retain only significant DEPs
dep.annotated_A <- dep.annotated_A[dep.annotated_A$DEP.Status != "Not significant", ]
dep.annotated_H <- dep.annotated_H[dep.annotated_H$DEP.Status != "Not significant", ]

# Write DEP tables to files
write.table(dep.annotated_A, file = "DEP_A549_annotated.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(dep.annotated_H, file = "DEP_H2009_annotated.csv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---------------- VOLCANO PLOTS ----------------
library(ggplot2)
library(ggrepel)

# A549 Volcano plot
volcano.df_A <- data.frame(
  Protein = dep.annotated_A$Protein,
  log2FC = dep.annotated_A$log2FC,
  negLog10P = -log10(dep.annotated_A$P.Value),
  DEP.Status = factor(dep.annotated_A$DEP.Status,
                      levels = c("Upregulated in Treatment", "Downregulated in Treatment", "Not significant"))
)

ggsave("Images/VolcanoPlot_A549.png",
       width = 6000, height = 4500, dpi = 600, units = "px") 

ggplot(volcano.df_A, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = DEP.Status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated in Treatment" = "red",
                                "Downregulated in Treatment" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(volcano.df_A, DEP.Status != "Not significant"),
                  aes(label = Protein),
                  size = 2.5, max.overlaps = Inf, box.padding = 0.5, segment.size = 0.2) +
  labs(title = "Volcano Plot – A549 (Treatment vs Control)",
       x = expression(log[2]*"(Fold-change)"),
       y = expression(-log[10]*"(p-value)"),
       color = "DEP Status") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

# Extract DEP IDs
deps_A <- dep_A[(log2fc_A > 1 & pval_A < 0.05) | (log2fc_A < -1 & pval_A < 0.05), ]
deps_A <- deps_A$ID

# H2009 Volcano plot
volcano.df_H <- data.frame(
  Protein = dep.annotated_H$Protein,
  log2FC = dep.annotated_H$log2FC,
  negLog10P = -log10(dep.annotated_H$P.Value),
  DEP.Status = factor(dep.annotated_H$DEP.Status,
                      levels = c("Upregulated in Treatment", "Downregulated in Treatment", "Not significant"))
)

ggplot(volcano.df_H, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = DEP.Status), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated in Treatment" = "red",
                                "Downregulated in Treatment" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = subset(volcano.df_H, DEP.Status != "Not significant"),
                  aes(label = Protein),
                  size = 2.5, max.overlaps = Inf, box.padding = 0.5, segment.size = 0.2) +
  labs(title = "Volcano Plot – H2009 (Treatment vs Control)",
       x = expression(log[2]*"(Fold-change)"),
       y = expression(-log[10]*"(p-value)"),
       color = "DEP Status") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")


# Extract DEP IDs
deps_H <- dep_H[(log2fc_H > 1 & pval_H < 0.05) | (log2fc_H < -1 & pval_H < 0.05), ]
deps_H <- deps_H$ID


## ---------------- VENN DIAGRAMS ----------------
library(VennDiagram)
library(grid)

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    A549 = deps_A,
    H2009 = deps_H
  ),
  filename = NULL,
  fill = c("#66a3ff", "#ffc04c"),
  col = "white",
  alpha = 0.7,
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20),
  cat.col = c("#1f78b4", "#d95f02"),
  main = "Overlap of DEPs between A549 and H2009",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plot)

## ---------------- COMPARISON WITH METABOANALYST ----------------

meta_DEP_A <- read_delim("MetaboAnalyst/DEP_LineaA_70Missings_FC2.csv")
meta_DEP_H <- read_delim("MetaboAnalyst/DEP_LineaH_70Missings_FC2.csv")

meta_DEP_A <- as.vector(meta_DEP_A[,1])$...1
meta_DEP_H <- as.vector(meta_DEP_H[,1])$...1

# A549
grid.newpage()
venn.plot_A <- venn.diagram(
  x = list(
    Limma = deps_A,
    MetaboAnalyst = meta_DEP_A
  ),
  filename = NULL,
  fill = c("#66a3ff", "#ffc04c"),
  col = "white",
  alpha = 0.7,
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20),
  cat.col = c("#1f78b4", "#d95f02"),
  main = "A549 – Overlap of DEPs (Limma vs MetaboAnalyst)",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plot_A)

# H2009
grid.newpage()
venn.plot_H <- venn.diagram(
  x = list(
    Limma = deps_H,
    MetaboAnalyst = meta_DEP_H
  ),
  filename = NULL,
  fill = c("#66a3ff", "#ffc04c"),
  col = "white",
  alpha = 0.7,
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 20),
  cat.col = c("#1f78b4", "#d95f02"),
  main = "H2009 – Overlap of DEPs (Limma vs MetaboAnalyst)",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plot_H)

## ---------------- HEATMAPS ----------------
library(pheatmap)

# Subset normalized expression matrix by DEPs
df_heatmap_A <- expr_A_qn[rownames(expr_A_qn) %in% deps_A, ]
df_heatmap_H <- expr_H_qn[rownames(expr_H_qn) %in% deps_H, ]

# Ensure unique row names
rownames(df_heatmap_A) <- make.unique(rownames(df_heatmap_A))
rownames(df_heatmap_H) <- make.unique(rownames(df_heatmap_H))

# Column annotation: treatment condition
conditions_A <- ifelse(grepl("^C_", colnames(df_heatmap_A)), "Control", "Irradiation")
annotations_A <- data.frame(Condition = factor(conditions_A))
rownames(annotations_A) <- colnames(df_heatmap_A)

conditions_H <- ifelse(grepl("^C_", colnames(df_heatmap_H)), "Control", "Irradiation")
annotations_H <- data.frame(Condition = factor(conditions_H))
rownames(annotations_H) <- colnames(df_heatmap_H)

# Heatmap A549
heatmap_A <- pheatmap(df_heatmap_A,
                      annotation_col = annotations_A,
                      clustering_distance_rows = "euclidean",
                      cluster_cols = FALSE,
                      clustering_method = "complete",
                      scale = "row",
                      fontsize_row = 10,
                      fontsize_col = 10,
                      main = "A549 Cell Line - DEPs between conditions")

png("Images/figure2c.png", width = 6000, height = 4500, res = 600)
grid.newpage()
grid.draw(heatmap_A$gtable)
dev.off()

# Heatmap H2009
heatmap_H <- pheatmap(df_heatmap_H,
                      annotation_col = annotations_H,
                      clustering_distance_rows = "euclidean",
                      cluster_cols = FALSE,
                      clustering_method = "complete",
                      scale = "row",
                      fontsize_row = 10,
                      fontsize_col = 10,
                      main = "H2009 Cell Line - DEPs between conditions")

png("Images/figure2d.png", width = 6000, height = 4500, res = 600)
grid.newpage()
grid.draw(heatmap_H$gtable)
dev.off()
