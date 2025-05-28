## Accessing subcellular location data from HPA and GeneCards

## Note: The default `hpaDownload()` function from the HPAanalyze package was manually adjusted
## due to outdated download URLs in the official release at the time of analysis.

hpaDownload_fixed <- function (downloadList = "histology", version = "latest") 
{
  allDatasets <- tibble(datasetnames = c("Normal tissue", "Pathology", 
                                         "Subcellular location", "RNA consensus tissue", "RNA HPA tissue", 
                                         "RNA GTEx tissue", "RNA FANTOM tissue", "RNA single cell type", 
                                         "RNA single cell type tissue cluster", "RNA GTEx brain region", 
                                         "RNA FANTOM brain region", "RNA pig brain region", "RNA pig brain subregion sample", 
                                         "RNA mouse brain region", "RNA mouse brain subregion sample", 
                                         "RNA Allen mouse brain region", "RNA HPA immune cell", 
                                         "RNA HPA immune cell sample", "RNA Monaco immune cell", 
                                         "RNA Schmiedel immune cell", "RNA HPA blood cell", "RNA HPA blood cell sample", 
                                         "RNA Monaco blood cell", "RNA Schmiedel blood cell", 
                                         "RNA HPA cell line cancer", "RNA HPA cell line", "RNA TCGA cancer sample", 
                                         "RNA transcript tissue", "RNA transcript GTEx retina", 
                                         "RNA transcript immune cells", "RNA transcript cell line", 
                                         "RNA transcript pig brain", "RNA transcript mouse brain"), 
                        tidycols = list(normal_tissue = c("ensembl", "gene", 
                                                          "tissue", "cell_type", "level", "reliability"), pathology = c("ensembl", 
                                                                                                                        "gene", "cancer", "high", "medium", "low", "not_detected", 
                                                                                                                        "prognostic_favorable", "unprognostic_favorable", 
                                                                                                                        "prognostic_unfavorable", "unprognostic_unfavorable"), 
                                        subcellular_location = c("ensembl", "gene", "reliability", 
                                                                 "main_location", "additional_location", "extracellular_location", 
                                                                 "enhanced", "supported", "approved", "uncertain", 
                                                                 "single_cell_var_intensity", "single_cell_var_spatial", 
                                                                 "cell_cycle_dependency", "go_id"), rna_tissue_consensus = c("ensembl", 
                                                                                                                             "gene", "tissue", "nx"), rna_tissue_hpa = c("ensembl", 
                                                                                                                                                                         "gene", "tissue", "tpm", "ptpm", "nx"), rna_tissue_gtex = c("ensembl", 
                                                                                                                                                                                                                                     "gene", "tissue", "tpm", "ptpm", "nx"), rna_tissue_fantom = c("ensembl", 
                                                                                                                                                                                                                                                                                                   "gene", "tissue", "tags_per_million", "scaled_tags_per_million", 
                                                                                                                                                                                                                                                                                                   "nx"), rna_single_cell_type = c("ensembl", "gene", 
                                                                                                                                                                                                                                                                                                                                   "cell_type", "nx"), rna_single_cell_type_tissue = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                       "gene", "tissue", "cluster", "cell_type", "read_count", 
                                                                                                                                                                                                                                                                                                                                                                                       "ptpm"), rna_brain_gtex = c("ensembl", "gene", 
                                                                                                                                                                                                                                                                                                                                                                                                                   "brain_region", "tpm", "ptpm", "nx"), rna_brain_fantom = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "gene", "brain_region", "tags_per_million", "scaled_tags_per_million", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "nx"), rna_pig_brain_hpa = c("ensembl", "gene", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           "brain_region", "tpm", "ptpm", "nx"), rna_pig_brain_sample_hpa = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "main_region", "subregion", "animal", "tpm", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "ptpm"), rna_mouse_brain_hpa = c("ensembl", "gene", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "brain_region", "tpm", "ptpm", "nx"), rna_mouse_brain_sample_hpa = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "main_region", "subregion", "animal", "tpm", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "ptpm"), rna_mouse_brain_allen = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "gene", "brain_region", "expression_energy"), 
                                        rna_immune_cell = c("ensembl", "gene", "immune_cell", 
                                                            "tpm", "ptpm", "ntpm"), rna_immune_cell_sample = c("sample_id", 
                                                                                                               "donor", "immune_cell", "ensembl", "gene", "tpm", 
                                                                                                               "ptpm", "ntpm"), rna_immune_cell_monaco = c("ensembl", 
                                                                                                                                                           "gene", "immune_cell", "tpm", "ptpm"), rna_immune_cell_schmiedel = c("ensembl", 
                                                                                                                                                                                                                                "gene", "immune_cell", "tpm"), rna_blood_cell = c("ensembl", 
                                                                                                                                                                                                                                                                                  "gene", "blood_cell", "tpm", "ptpm", "nx"), rna_blood_cell_sample = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                        "blood_cell_type", "donor", "ptpm", "nx"), rna_blood_cell_monaco = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                             "gene", "blood_cell", "tpm", "ptpm"), rna_blood_cell_schmiedel = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "gene", "blood_cell", "tpm"), rna_celline_cancer = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "gene", "cancer", "tpm", "ptpm", "ntpm"), rna_celline = c("ensembl", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "gene", "cell_line", "tpm", "ptpm", "ntpm"), 
                                        rna_cancer_sample = c("ensembl", "sample", "cancer", 
                                                              "fpkm"), transcript_rna_tissue = c("ensgid", 
                                                                                                 "enstid", "sample", "tpm"), transcript_rna_gtexretina = c("ensgid", 
                                                                                                                                                           "enstid", "sample", "tpm"), transcript_rna_immunecells = c("ensgid", 
                                                                                                                                                                                                                      "enstid", "sample", "tpm"), transcript_rna_celline = c("ensgid", 
                                                                                                                                                                                                                                                                             "enstid", "sample", "tpm"), transcript_rna_pigbrain = c("ensgid", 
                                                                                                                                                                                                                                                                                                                                     "enstid", "sample", "tpm"), transcript_rna_mousebrain = c("ensgid", 
                                                                                                                                                                                                                                                                                                                                                                                               "enstid", "sample", "tpm")), urls = c(normal_tissue = "https://www.proteinatlas.org/download/tsv/normal_tissue.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     pathology = "https://www.proteinatlas.org/download/tsv/pathology.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     subcellular_location = "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_tissue_consensus = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_tissue_hpa = "https://www.proteinatlas.org/download/tsv/rna_tissue_hpa.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_tissue_gtex = "https://www.proteinatlas.org/download/tsv/rna_tissue_gtex.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_tissue_fantom = "https://www.proteinatlas.org/download/tsv/rna_tissue_fantom.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_single_cell_type = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_single_cell_type_tissue = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type_tissue.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_brain_gtex = "https://www.proteinatlas.org/download/tsv/rna_brain_gtex.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_brain_fantom = "https://www.proteinatlas.org/download/tsv/rna_brain_fantom.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_pig_brain_hpa = "https://www.proteinatlas.org/download/tsv/rna_pig_brain_hpa.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_pig_brain_sample_hpa = "https://www.proteinatlas.org/download/tsv/rna_pig_brain_sample_hpa.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_mouse_brain_hpa = "https://www.proteinatlas.org/download/tsv/rna_mouse_brain_hpa.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_mouse_brain_sample_hpa = "https://www.proteinatlas.org/download/tsv/rna_mouse_brain_sample_hpa.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_mouse_brain_allen = "https://www.proteinatlas.org/download/tsv/rna_mouse_brain_allen.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_immune_cell = "https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_immune_cell_sample = "https://www.proteinatlas.org/download/tsv/rna_immune_cell_sample.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_immune_cell_monaco = "https://www.proteinatlas.org/download/tsv/rna_immune_cell_monaco.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_immune_cell_schmiedel = "https://www.proteinatlas.org/download/tsv/rna_immune_cell_schmiedel.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_blood_cell = "https://v21.proteinatlas.org/download/tsv/rna_blood_cell.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_blood_cell_sample = "https://v21.proteinatlas.org/download/tsv/rna_blood_cell_sample.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_blood_cell_monaco = "https://v21.proteinatlas.org/download/tsv/rna_blood_cell_monaco.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_blood_cell_schmiedel = "https://v21.proteinatlas.org/download/tsv/rna_blood_cell_schmiedel.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_celline_cancer = "https://www.proteinatlas.org/download/tsv/rna_celline_cancer.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_celline = "https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     rna_cancer_sample = "https://www.proteinatlas.org/download/tsv/rna_cancer_sample.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_tissue = "https://www.proteinatlas.org/download/tsv/transcript_rna_tissue.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_gtexretina = "https://www.proteinatlas.org/download/tsv/transcript_rna_gtexretina.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_immunecells = "https://www.proteinatlas.org/download/tsv/transcript_rna_immunecells.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_celline = "https://v21.proteinatlas.org/download/tsv/transcript_rna_celline.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_pigbrain = "https://www.proteinatlas.org/download/tsv/transcript_rna_pigbrain.tsv.zip", 
                                                                                                                                                                                                                                                                                                                                                                                                                                     transcript_rna_mousebrain = "https://www.proteinatlas.org/download/tsv/transcript_rna_mousebrain.tsv.zip"))
  replace_shortcut <- function(x, shortcut, with) {
    x <- rep(x, 1 + (length(with) - 1) * (x == shortcut))
    x[x == shortcut] <- with
    return(x)
  }
  downloadList <- downloadList %>% replace_shortcut("all", 
                                                    allDatasets$datasetnames) %>% replace_shortcut("histology", 
                                                                                                   allDatasets$datasetnames[1:3]) %>% replace_shortcut("rna tissue", 
                                                                                                                                                       allDatasets$datasetnames[4:7]) %>% replace_shortcut("rna cell type", 
                                                                                                                                                                                                           allDatasets$datasetnames[8:9]) %>% replace_shortcut("rna brain region", 
                                                                                                                                                                                                                                                               allDatasets$datasetnames[10:16]) %>% replace_shortcut("rna immune cell", 
                                                                                                                                                                                                                                                                                                                     allDatasets$datasetnames[17:20]) %>% replace_shortcut("rna blood cell", 
                                                                                                                                                                                                                                                                                                                                                                           allDatasets$datasetnames[21:24]) %>% replace_shortcut("isoform", 
                                                                                                                                                                                                                                                                                                                                                                                                                                 allDatasets$datasetnames[28:33])
  downloadDatasets <- filter(allDatasets, datasetnames %in% 
                               downloadList)
  loadedData <- list()
  if (version %in% c("example", "built-in")) {
    message("Only the followings are example/built-in datasets: \n - Normal tissue \n - Pathology \n - Subcellular location \nOther datasets will not be loaded")
    downloadDatasets <- filter(downloadDatasets, datasetnames %in% 
                                 c("Normal tissue", "Pathology", "Subcellular location"))
    for (i in names(downloadDatasets$urls)) {
      loadedData[[i]] <- hpa_histology_data[[i]]
    }
  }
  else if (version == "latest") {
    for (i in seq_along(downloadDatasets$urls)) {
      temp <- tempfile()
      download.file(url = downloadDatasets$urls[[i]], destfile = temp)
      loadedData[[i]] <- read.delim2(unz(temp, unzip(temp, 
                                                     list = TRUE)$Name[1]), stringsAsFactors = FALSE, 
                                     check.names = FALSE, strip.white = TRUE, sep = "\t", 
                                     na.strings = c("", " "))
      unlink(temp)
      if (downloadDatasets$datasetnames[[i]] %in% c("RNA transcript tissue", 
                                                    "RNA transcript GTEx retina", "RNA transcript immune cells", 
                                                    "RNA transcript cell line", "RNA transcript pig brain", 
                                                    "RNA transcript mouse brain")) {
        loadedData[[i]] <- stats::reshape(loadedData[[i]], 
                                          direction = "long", varying = list(3:ncol(loadedData[[i]])), 
                                          v.names = "tpm", timevar = "sample", times = c(colnames(loadedData[[i]][, 
                                                                                                                  3:ncol(loadedData[[i]])]))) %>% subset(select = -id)
      }
      colnames(loadedData[[i]]) <- downloadDatasets$tidycols[[i]]
    }
    loadedData <- lapply(loadedData, as_tibble)
    names(loadedData) <- downloadDatasets$urls %>% gsub(".tsv.zip|https://www.proteinatlas.org/download/tsv/|https://v21.proteinatlas.org/download/tsv/", 
                                                        "", .)
  }
  return(loadedData)
}


## Access subcellular location data from HPA and GeneCards

## Note: The default hpaDownload() function was modified (hpaDownload_fixed())
## due to outdated URLs in the original HPAanalyze package version

library(BiocStyle)
library(HPAanalyze)
library(dplyr)

## Download HPA subcellular location data
subcell_locations <- hpaDownload_fixed(downloadList = 'Subcellular location')

## Visualize subcellular locations for DEPs
localization_DEPA <- hpaVisSubcell(
  data = subcell_locations, 
  targetGene = deps_A,
  reliability = c("enhanced", "supported", "approved")
)

localization_DEPH <- hpaVisSubcell(
  data = subcell_locations, 
  targetGene = deps_H,
  reliability = c("enhanced", "supported", "approved")
)

localization_DEPA
localization_DEPH

###############################################################################
## Load subcellular localization data from GeneCards

library(readxl)

# Load data from different GeneCards sheets
DEP_uniprot_locationA <- read_excel("GeneCards/All_DEP_A.xlsx", sheet = "UniProtSubcellularLocations")
DEP_genecard_locationA <- read_excel("GeneCards/All_DEP_A.xlsx", sheet = "CompartmentsSubcellularLocation")
DEP_HPA_locationA <- read_excel("GeneCards/All_DEP_A.xlsx", sheet = "HPASubcellularLocations")
DEP_CellularComponentsA <- read_excel("GeneCards/All_DEP_A.xlsx", sheet = "CellularComponents")

DEP_uniprot_locationH <- read_excel("GeneCards/All_DEP_H.xlsx", sheet = "UniProtSubcellularLocations")
DEP_genecard_locationH <- read_excel("GeneCards/All_DEP_H.xlsx", sheet = "CompartmentsSubcellularLocation")
DEP_HPA_locationH <- read_excel("GeneCards/All_DEP_H.xlsx", sheet = "HPASubcellularLocations")
DEP_CellularComponentsH <- read_excel("GeneCards/All_DEP_H.xlsx", sheet = "CellularComponents")

## Filter entries related to membrane localization

library(stringr)

# UniProt – A549
filt_locations_uniprotA <- DEP_uniprot_locationA |>
  filter(str_detect(tolower(Compartments), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_uniprotA$InputTerm))

# UniProt – H2009
filt_locations_uniprotH <- DEP_uniprot_locationH |>
  filter(str_detect(tolower(Compartments), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_uniprotH$InputTerm))

# Cellular Components – A549
filt_locations_CCA <- DEP_CellularComponentsA |>
  filter(str_detect(tolower(Term), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_CCA$InputTerm))

# Cellular Components – H2009
filt_locations_CCH <- DEP_CellularComponentsH |>
  filter(str_detect(tolower(Term), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_CCH$InputTerm))

# GeneCards Compartments – A549
filt_locations_GCA <- DEP_genecard_locationA |>
  filter(Confidence >= 3) |>
  filter(str_detect(tolower(Name), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_GCA$InputTerm))

# GeneCards Compartments – H2009
filt_locations_GCH <- DEP_genecard_locationH |>
  filter(Confidence >= 3) |>
  filter(str_detect(tolower(Name), "membrane|extracellular|apical|cell"))
length(unique(filt_locations_GCH$InputTerm))

# HPA – A549
filt_locations_HPAA <- DEP_HPA_locationA |>
  filter(str_detect(tolower(Text), "membrane|cell")) |>
  filter(str_detect(tolower(Description), "approved|enhanced|supported")) |>
  filter(str_detect(tolower(Score), "4|5"))
length(unique(filt_locations_HPAA$InputTerm))

# HPA – H2009
filt_locations_HPAH <- DEP_HPA_locationH |>
  filter(str_detect(tolower(Text), "membrane|cell")) |>
  filter(str_detect(tolower(Description), "approved|enhanced|supported")) |>
  filter(str_detect(tolower(Score), "4|5"))
length(unique(filt_locations_HPAH$InputTerm))

###############################################################################
## Venn Diagrams of subcellular annotations
library(VennDiagram)
library(grid)

# A549
protein_uniprotA <- unique(filt_locations_uniprotA$InputTerm)
protein_CCA <- unique(filt_locations_CCA$InputTerm)
protein_GCA <- unique(filt_locations_GCA$InputTerm)
protein_HPAA <- unique(filt_locations_HPAA)

grid.newpage()
venn.plotA <- venn.diagram(
  x = list(
    UniProt = protein_uniprotA,
    CellularComponents = protein_CCA,
    GeneCards = protein_GCA
  ),
  filename = NULL,
  fill = c("#66a3ff", "#ffc04c", "#a1d99b"),
  col = "white",
  alpha = 0.7,
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03, 0.03),
  cat.pos = c(-20, 20, 0),
  cat.col = c("#1f78b4", "#d95f02", "#238b45"),
  main = "A549 – Overlap of subcellular annotations",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plotA)

# H2009
protein_uniprotH <- unique(filt_locations_uniprotH$InputTerm)
protein_CCH <- unique(filt_locations_CCH$InputTerm)
protein_GCH <- unique(filt_locations_GCH$InputTerm)
protein_HPAH <- unique(filt_locations_HPAH)

grid.newpage()
venn.plotH <- venn.diagram(
  x = list(
    UniProt = protein_uniprotH,
    CellularComponents = protein_CCH,
    GeneCards = protein_GCH
  ),
  filename = NULL,
  fill = c("#66a3ff", "#ffc04c", "#a1d99b"),
  col = "white",
  alpha = 0.7,
  cex = 1.4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03, 0.03),
  cat.pos = c(-20, 20, 0),
  cat.col = c("#1f78b4", "#d95f02", "#238b45"),
  main = "H2009 – Overlap of subcellular annotations",
  main.cex = 1.5,
  main.fontface = "bold",
  margin = 0.05
)
grid.draw(venn.plotH)

###############################################################################
## Identify proteins annotated in at least 3 databases

# Create list of sources
db_list_A <- list(UniProt = protein_uniprotA, CCA = protein_CCA, GCA = protein_GCA, HPA = protein_HPAA)
db_list_H <- list(UniProt = protein_uniprotH, CCA = protein_CCH, GCA = protein_GCH, HPA = protein_HPAH)

# Count number of databases per gene
all_protein_A <- unique(unlist(db_list_A))
all_protein_H <- unique(unlist(db_list_H))

counts_A <- sapply(all_protein_A, function(g) sum(sapply(db_list_A, function(x) g %in% x)))
counts_H <- sapply(all_protein_H, function(g) sum(sapply(db_list_H, function(x) g %in% x)))

protein_in_3plus_A <- names(counts_A)[counts_A >= 3]
protein_in_3plus_H <- names(counts_H)[counts_H >= 3]

###############################################################################
## Cross-reference with upregulated DEPs

upregulated_A <- dep_A[(log2fc_A > 1 & pval_A < 0.05), ]$ID
upregulated_H <- dep_H[(log2fc_H > 1 & pval_H < 0.05), ]$ID

candidates_A <- intersect(upregulated_A, protein_in_3plus_A)
candidates_H <- intersect(upregulated_H, protein_in_3plus_H)

## Create Table 2 – Candidate membrane-associated DEPs
filt_total <- bind_rows(
  filt_locations_uniprotA |> mutate(Source = "UniProt-A549"),
  filt_locations_uniprotH |> mutate(Source = "UniProt-H2009"),
  filt_locations_CCA       |> mutate(Source = "CellComp-A549"),
  filt_locations_CCH       |> mutate(Source = "CellComp-H2009"),
  filt_locations_GCA       |> mutate(Source = "GeneCards-A549"),
  filt_locations_GCH       |> mutate(Source = "GeneCards-H2009"),
  filt_locations_HPAA      |> mutate(Source = "HPA-A549"),
  filt_locations_HPAH      |> mutate(Source = "HPA-H2009")
)

candidates <- unique(c(candidates_A, candidates_H))
table2 <- filt_total |> filter(InputTerm %in% candidates)

## Create Table S1 – All genes present in 3+ databases
all_3bases <- unique(c(protein_in_3plus_A, protein_in_3plus_H))
table1 <- filt_total |> filter(InputTerm %in% all_3bases)
