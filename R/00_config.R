## ============================================================================
## 00_config.R — shared configuration and helper functions
##
## Sourced by every analysis script. Defines paths, loads and pre-processes the
## model, and implements the definitions used throughout the manuscript:
##   - EC/gene confidence categories 1-4            (Methods 2.4)
##   - EC specificity sub-categories a / b          (Methods 2.4)
##   - the three counting units: entries, distinct annotated transformations,
##     distinct stoichiometries                     (Sections 2.4 and 3.1)
##
## Paths are RELATIVE to the repository root: set the working directory to the
## repository root before sourcing, or open the .Rproj file.
## ============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(tidyr)
  library(purrr);  library(readr); library(ggplot2); library(forcats)
  library(scales); library(patchwork)
})

## ---- paths -----------------------------------------------------------------
DATA_DIR   <- Sys.getenv("YGEM_DATA",   unset = "data")
OUTPUT_DIR <- Sys.getenv("YGEM_OUTPUT", unset = "output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MODEL_XLSX <- file.path(DATA_DIR, "yeast-GEM.xlsx")

## Optional external files. Scripts that need them fail gracefully if absent.
find_data <- function(pattern) {
  f <- list.files(DATA_DIR, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(f) == 0) NA_character_ else f[1]
}
UNIPROT_TSV <- file.path(DATA_DIR, "uniprot_yeast_reviewed.tsv")
SGD_GAF     <- find_data("\\.gaf(\\.gz)?$")
SGD_GFF     <- find_data("\\.gff(\\.gz)?$")
GO_OBO      <- file.path(DATA_DIR, "go-basic.obo")

## ---- plotting defaults -----------------------------------------------------
## No titles or sub-titles inside plots: captions carry that information.
COL <- list(
  primary   = "#1F4E79",
  secondary = "#E69F00",
  accent    = "#C0392B",
  grey      = "#7F8C8D",
  light     = "#D9D9D9",
  score     = c("0" = "#EFF3FA", "1" = "#C6DBEF", "2" = "#6BAED6",
                "3" = "#2171B5", "Unscored" = "#9E9E9E")
)

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_blank(),
          plot.subtitle    = element_blank(),
          plot.tag         = element_text(face = "bold", size = base_size + 2))
}

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
  invisible(plot)
}

## ---- field normalisation ---------------------------------------------------
norm_field <- function(x) {
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "na")] <- NA_character_
  x
}

## Systematic ORF identifiers: YxLnnnW / YxRnnnC (+ optional -A), or Qnnnn (mt)
GENE_RE <- "Y[A-P][LR][0-9]{3}[WC](-[A-Z])?|Q[0-9]{4}"

extract_genes <- function(gpr) {
  if (is.na(gpr)) return(character(0))
  unique(toupper(str_extract_all(gpr, GENE_RE)[[1]]))
}

split_ec <- function(ec) {
  if (is.na(ec)) return(character(0))
  p <- str_trim(str_split(ec, "[;,\\s]+")[[1]])
  unique(p[str_detect(p, "^[0-9]")])
}

## An EC number is FULLY SPECIFIED when all four levels are resolved (1.1.1.4);
## it is PARTIAL when one or more sub-class levels are undetermined (3.1.3.-).
is_full_ec <- function(x) str_detect(x, "^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$")

## ---- model loading ---------------------------------------------------------
load_model <- function(path = MODEL_XLSX) {

  if (!file.exists(path))
    stop("Model file not found: ", normalizePath(path, mustWork = FALSE),
         "\nSee data/README.md for the download link.")

  rxns  <- read_excel(path, sheet = "RXNS")
  comps <- read_excel(path, sheet = "COMPS")

  comp_map     <- setNames(comps$NAME, comps$ABBREVIATION)
  valid_comps  <- comps$ABBREVIATION[order(-nchar(comps$ABBREVIATION))]
  comp_pattern <- paste0("\\[(", paste(valid_comps, collapse = "|"), ")\\]")

  get_comps  <- function(eq) {
    if (is.na(eq) || str_trim(eq) == "") return(character(0))
    sort(unique(str_remove_all(str_extract_all(eq, comp_pattern)[[1]], "\\[|\\]")))
  }
  strip_comp <- function(eq) if (is.na(eq)) NA_character_ else
    str_squish(str_remove_all(eq, comp_pattern))

  out <- rxns %>%
    filter(!is.na(EQUATION), str_trim(EQUATION) != "") %>%
    mutate(
      EC   = norm_field(`EC-NUMBER`),
      GPR  = norm_field(`GENE ASSOCIATION`),

      ## --- EC/gene confidence categories (Methods 2.4) ---------------------
      category = case_when(
        !is.na(EC) & !is.na(GPR) ~ "Cat 1 (EC+gene)",
        !is.na(EC) &  is.na(GPR) ~ "Cat 2 (EC only)",
         is.na(EC) & !is.na(GPR) ~ "Cat 3 (gene only)",
        TRUE                      ~ "Cat 4 (neither)"),
      cat_num = as.integer(str_extract(category, "[0-9]")),

      ## --- EC specificity: sub-categories a / b (Methods 2.4) --------------
      ec_list  = map(EC, split_ec),
      n_full   = map_int(ec_list, ~ sum(is_full_ec(.x))),
      n_part   = map_int(ec_list, ~ sum(!is_full_ec(.x))),
      ec_spec  = case_when(n_full > 0 ~ "fully specified",
                           n_part > 0 ~ "partial only",
                           TRUE       ~ "none"),
      sub_category = case_when(
        cat_num %in% c(1, 2) & ec_spec == "fully specified" ~ paste0(cat_num, "a"),
        cat_num %in% c(1, 2)                                 ~ paste0(cat_num, "b"),
        TRUE                                                 ~ as.character(cat_num)),

      ## --- compartments ----------------------------------------------------
      comp_set   = map(EQUATION, get_comps),
      n_comps    = lengths(comp_set),
      compartment = map_chr(comp_set, ~ if (length(.x) == 1) .x else NA_character_),
      comp_pair   = map_chr(comp_set, ~ paste(.x, collapse = "+")),
      span        = ifelse(n_comps == 1, "intra-compartmental", "inter-compartmental"),

      ## --- deduplication keys (Section 3.1) --------------------------------
      ## eq_nocomp        : stoichiometry with compartment tags removed.
      ##                    Reaction directionality is part of the identity:
      ##                    "A => B" and "A <=> B" are distinct.
      ## annot_key        : stoichiometry + gene set + EC set. Gene and EC order
      ##                    is normalised so that "X or Y" == "Y or X".
      eq_nocomp = map_chr(EQUATION, strip_comp),
      gene_set  = map_chr(GPR, ~ paste(sort(extract_genes(.x)), collapse = ";")),
      ec_set    = map_chr(ec_list, ~ paste(sort(.x), collapse = ";")),
      annot_key = paste(eq_nocomp, gene_set, ec_set, sep = "||"),

      ## --- native SBML/COBRA evidence score --------------------------------
      native_score = ifelse(is.na(`CONFIDENCE SCORE`), "Unscored",
                            as.character(`CONFIDENCE SCORE`))
    )

  attr(out, "comp_map")     <- comp_map
  attr(out, "comp_pattern") <- comp_pattern
  out
}

## ---- shared summaries ------------------------------------------------------
CAT_LEVELS   <- c("Cat 1 (EC+gene)", "Cat 2 (EC only)",
                  "Cat 3 (gene only)", "Cat 4 (neither)")
SCORE_LEVELS <- c("0", "1", "2", "3", "Unscored")
CHR_MAP      <- c(A="I", B="II", C="III", D="IV", E="V", `F`="VI", G="VII", H="VIII",
                  I="IX", J="X", K="XI", L="XII", M="XIII", N="XIV", O="XV", P="XVI")
CHR_LEVELS   <- c(as.character(as.roman(1:16)), "Mitochondrial")

chromosome_of <- function(gene_id) {
  case_when(str_starts(gene_id, "Q")           ~ "Mitochondrial",
            str_detect(gene_id, "^Y[A-P][LR]") ~ unname(CHR_MAP[str_sub(gene_id, 2, 2)]),
            TRUE                                ~ NA_character_)
}

write_out <- function(x, name) {
  write_csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")))
  invisible(x)
}

banner <- function(txt) cat("\n", strrep("=", 74), "\n", txt, "\n",
                            strrep("=", 74), "\n", sep = "")
