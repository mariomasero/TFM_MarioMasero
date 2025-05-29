## Load required libraries
library(readr)

## Read proteomics results file
ProteomicsResults <- read_delim("ResultadosProteomica.csv", 
                                delim = ";", escape_double = FALSE, 
                                trim_ws = TRUE)

## Extract protein accession IDs
protein_ids <- ProteomicsResults$`Protein Accession`

## Create expression matrix for A549 cell line
expr_A <- as.matrix(ProteomicsResults[ ,(5:14)])  # Extract LFQ values
rownames(expr_A) <- protein_ids  # Assign protein names

## Create expression matrix for H2009 cell line
expr_H <- as.matrix(ProteomicsResults[ ,(15:24)])  # Extract LFQ values
rownames(expr_H) <- protein_ids  # Assign protein names

## Rename columns
colnames(expr_A) <- c("C_1","C_2","C_3", "C_4", "C_5",
                      "T_1","T_2","T_3", "T_4", "T_5")

colnames(expr_H) <- c("C_1","C_2","C_3", "C_4", "C_5",
                      "T_1","T_2","T_3", "T_4", "T_5")

## Preview first rows of each matrix
head(expr_A)
head(expr_H)

## Define Y-axis limits for histograms
ymax_raw <- max(hist(as.numeric(expr_A), plot = FALSE)$counts,
                hist(as.numeric(expr_H), plot = FALSE)$counts)

xlim_raw <- range(c(as.numeric(expr_A), as.numeric(expr_H)), na.rm = TRUE)

## Visual inspection of LFQ value distribution (to check log scale)
png("Images/histogram.png", width = 6000, height = 3700, res = 600)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(as.numeric(expr_A), 
     main = "A549 - Raw LFQ values", 
     xlab = "LFQ", col = "#1f77b4", border = "white",
     ylim = c(0, ymax_raw), xlim = xlim_raw,
     cex.main = 0.95, cex.lab = 0.95)

hist(as.numeric(expr_H), 
     main = "H2009 - Raw LFQ values", 
     xlab = "LFQ", col = "#ff7f0e", border = "white",
     ylim = c(0, ymax_raw), xlim = xlim_raw,
     cex.main = 0.95, cex.lab = 0.95)

## Summary statistics of raw expression values
summary(as.numeric(expr_A))
summary(as.numeric(expr_H))

## Apply log2 transformation (add 1 to avoid log2(0))
expr_A_log2 <- log2(expr_A + 1)
expr_H_log2 <- log2(expr_H + 1)

## Define new Y-axis limits for log2 histograms
ymax_log2 <- max(hist(as.numeric(expr_A_log2), plot = FALSE)$counts,
                 hist(as.numeric(expr_H_log2), plot = FALSE)$counts)

xlim_log2 <- range(c(as.numeric(expr_A_log2), as.numeric(expr_H_log2)), na.rm = TRUE)

## Visual inspection of log2-transformed distributions
hist(as.numeric(expr_A_log2), 
     main = "A549 - Log2(LFQ + 1)", 
     xlab = "log2(LFQ + 1)", col = "#1f77b4", border = "white",
     ylim = c(0, ymax_log2), xlim = xlim_log2,
     cex.main = 0.95, cex.lab = 0.95)

hist(as.numeric(expr_H_log2), 
     main = "H2009 - Log2(LFQ + 1)", 
     xlab = "log2(LFQ + 1)", col = "#ff7f0e", border = "white",
     ylim = c(0, ymax_log2), xlim = xlim_log2,
     cex.main = 0.95, cex.lab = 0.95)

par(mfrow = c(1, 1))
dev.off()

## Compute total protein abundance per sample (sum of LFQ intensities)
abundance_A <- apply(expr_A, 2, sum)
abundance_H <- apply(expr_H, 2, sum)

## Create metadata tables for plotting
df_A <- data.frame(
  Sample = colnames(expr_A),
  Abundance = abundance_A,
  Condition = ifelse(grepl("C_", colnames(expr_A)), "Control", "Irradiated"),
  CellLine = "A549"
)

df_H <- data.frame(
  Sample = colnames(expr_H),
  Abundance = abundance_H,
  Condition = ifelse(grepl("C_", colnames(expr_H)), "Control", "Irradiated"),
  CellLine = "H2009"
)
