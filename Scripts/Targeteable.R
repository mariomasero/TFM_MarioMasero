## Load drug annotations associated with DEPs
library(readxl)

DEP_DrugsA <- read_excel("GeneCards/All_DEP_A.xlsx", sheet = "UnifiedDrugs")
DEP_DrugsH <- read_excel("GeneCards/All_DEP_H.xlsx", sheet = "UnifiedDrugs")

## Filter drugs targeting candidate proteins
library(dplyr)

filt_drugsA <- DEP_DrugsA |>
  filter(InputTerm %in% candidates_A) |>
  filter(str_detect(tolower(Status), "approved|investigational"))

filt_drugsH <- DEP_DrugsH |>
  filter(InputTerm %in% candidates_H) |>
  filter(str_detect(tolower(Status), "approved|investigational"))

## Export filtered tables
write.csv(filt_drugsA, "drugs_targetsA.csv")
write.csv(filt_drugsH, "drugs_targetsH.csv")

## Create Table 3 – Targetable candidates by approved/investigational drugs
filt_drugs <- rbind(
  filt_drugsA |> mutate(Cell_line = "A549"),
  filt_drugsH |> mutate(Cell_line = "H2009")
)
