## ============================================================================
## 10_external_resources.R — Section 2.2
##
## Generates Supplementary Table S1: every external resource used in the study
## with its version, release date and file type.
##
## This table is the single authoritative record of the resource versions. Keep
## it in step with the versions cited in Section 2.2 and with the ExplorEnz
## snapshot in 07_ec_class_distribution.R: if one is updated, update the others.
## ============================================================================

source("R/00_config.R")

resources <- tribble(
  ~Resource,                                     ~Version,               ~`Release date`, ~`File type`, ~`Used for`,
  "Yeast-GEM",                                   "v9.0.2",               "2024-11-23",    "XLSX",       "Model source (Section 2.1)",
  "UniProtKB / Swiss-Prot",                      "2026_02",              "2026-06-10",    "TSV",        "EC numbers, GO terms, Pfam, localisation (2.2, 2.5, 2.6, 2.7)",
  "SGD, S288C reference genome annotation",      "R64.5.1",              "2024-05-29",    "GFF3",       "Per-chromosome gene counts (2.3)",
  "SGD / GO Consortium, gene association file",  "R64.5.1",              "2024-05-29",    "GAF",        "GO Cellular Component evidence codes (2.6)",
  "KEGG",                                        "118.2",                "2026-06-01",    "flat file",  "Reaction and pathway context (2.2)",
  "MetaCyc",                                     "29.6",                 "2026-02-26",    "flat file",  "Reaction definitions (2.2)",
  "BioCyc",                                      "29.6",                 "2026-02-26",    "flat file",  "Organism-specific pathway context (2.2)",
  "ExplorEnz (IUBMB Enzyme List)",               "6,972 active entries", "accessed 2026-07", "HTML",     "Global EC-class statistics (2.2, 2.7)",
  "R",                                           "4.3.2",                "-",             "-",          "All analyses (2.1)"
)

banner("Supplementary Table S1 - external resources")
print(as.data.frame(resources %>% select(Resource, Version, `Release date`, `File type`)))

write_out(resources, "TableS1_external_resources")

## ---- consistency check against the model file ------------------------------
if (file.exists(MODEL_XLSX)) {
  mi <- read_excel(MODEL_XLSX, sheet = "MODEL")
  cat("\nModel sheet metadata:\n")
  print(as.data.frame(mi %>% select(any_of(c("ID", "NAME", "TAXONOMY")))))
}

cat("\nNote: MetaCyc and BioCyc share a common version numbering, so the two\n")
cat("entries must always carry the same value.\n")
cat("ExplorEnz has no numbered releases; the entry count and the access date\n")
cat("are reported instead.\n")

message("\n10_external_resources.R done.")
