## ============================================================================
## 05_category3_recovery.R — Section 3.2
##
## Category-3 reactions (gene association but no EC number): mapping to reviewed
## UniProtKB entries and identification of the reactions that can be upgraded to
## Category 1 by importing an EC number that already exists in UniProtKB.
## Produces Supplementary Table S3.
##
## Requires data/uniprot_yeast_reviewed.tsv. If absent, the script downloads it
## from the UniProtKB REST API and caches it; if the network is unreachable it
## stops with an explanatory message rather than producing partial results.
## ============================================================================

source("R/00_config.R")
library(httr)

m     <- load_model()
genes <- read_excel(MODEL_XLSX, sheet = "GENES")
comp_map <- attr(m, "comp_map")

## ---- UniProtKB reviewed S. cerevisiae proteome -----------------------------
get_uniprot <- function(path = UNIPROT_TSV) {
  if (file.exists(path)) {
    message("Using cached UniProtKB export: ", path)
    return(read_tsv(path, show_col_types = FALSE))
  }
  message("Querying UniProtKB REST API ...")
  r <- tryCatch(
    GET("https://rest.uniprot.org/uniprotkb/stream",
        query = list(query  = "organism_id:559292 AND reviewed:true",
                     format = "tsv",
                     fields = paste("accession", "gene_primary", "go_c", "go_f",
                                    "ec", "xref_pfam", "cc_subcellular_location",
                                    "reviewed", sep = ",")),
        timeout(180)),
    error = function(e) NULL)
  if (is.null(r) || http_error(r))
    stop("UniProtKB is unreachable and no cached file is present.\n",
         "Download the reviewed S. cerevisiae proteome to ", path,
         " and re-run (see data/README.md).")
  x <- read_tsv(content(r, "text", encoding = "UTF-8"), show_col_types = FALSE)
  write_tsv(x, path)
  message("Cached to ", path)
  x
}

uniprot <- get_uniprot() %>% distinct(Entry, .keep_all = TRUE)
banner("Category 3: EC recovery from UniProtKB")
cat(sprintf("UniProtKB reviewed entries: %d\n", nrow(uniprot)))

gene_to_acc <- setNames(str_extract(genes$MIRIAM, "(?<=uniprot/)[A-Za-z0-9]+"),
                        genes$NAME)

## ---- Category-3 reactions and their genes ----------------------------------
c3 <- m %>% filter(cat_num == 3)

c3_map <- c3 %>%
  transmute(rxn_id = ID, gene_id = map(GPR, extract_genes)) %>%
  unnest(gene_id)

c3_genes <- distinct(c3_map, gene_id) %>%
  mutate(accession = unname(gene_to_acc[gene_id]))

cat(sprintf("Category-3 reactions: %d  (%d distinct stoichiometries, %d distinct annotated transformations)\n",
            nrow(c3), n_distinct(c3$eq_nocomp), n_distinct(c3$annot_key)))
cat(sprintf("Distinct genes: %d, of which %d (%.1f%%) carry a UniProt accession\n",
            nrow(c3_genes), sum(!is.na(c3_genes$accession)),
            100 * mean(!is.na(c3_genes$accession))))

in_uniprot <- sum(c3_genes$accession %in% uniprot$Entry, na.rm = TRUE)
cat(sprintf("Confirmed as reviewed (Swiss-Prot) entries: %d\n", in_uniprot))

## ---- is the gene association uniform across compartments? ------------------
c3_eq  <- unique(c3$eq_nocomp)
in_c3  <- m %>% filter(eq_nocomp %in% c3_eq) %>%
  group_by(eq_nocomp) %>%
  summarise(all_c3 = all(cat_num == 3),
            others = paste(sort(unique(cat_num[cat_num != 3])), collapse = "/"),
            .groups = "drop")
cat(sprintf("\nOf the %d stoichiometries involved, %d carry a gene in every representation,\n",
            nrow(in_c3), sum(in_c3$all_c3)))
cat(sprintf("  whereas %d include at least one representation without a gene:\n",
            sum(!in_c3$all_c3)))
print(count(filter(in_c3, !all_c3), others))

## ---- which genes already have an EC in UniProtKB? --------------------------
up_ec <- uniprot %>% select(accession = Entry, ec_uniprot = `EC number`,
                            gene_name = `Gene Names (primary)`)

c3_genes <- c3_genes %>%
  left_join(up_ec, by = "accession") %>%
  mutate(has_ec  = !is.na(ec_uniprot) & str_trim(ec_uniprot) != "",
         is_ec7  = map_lgl(ec_uniprot, ~ !is.na(.x) &&
                             any(str_starts(str_trim(str_split(.x, ";")[[1]]), "7"))))

cat(sprintf("\nGenes with a curated EC in UniProtKB: %d / %d (%.1f%%)\n",
            sum(c3_genes$has_ec), nrow(c3_genes), 100 * mean(c3_genes$has_ec)))
cat(sprintf("  of which carrying an EC 7 (translocase) assignment: %d\n",
            sum(c3_genes$is_ec7, na.rm = TRUE)))
if (any(c3_genes$is_ec7, na.rm = TRUE))
  print(c3_genes %>% filter(is_ec7) %>% select(gene_id, gene_name, accession, ec_uniprot))

## ---- recoverable reactions -------------------------------------------------
rec <- c3_map %>%
  left_join(c3_genes %>% select(gene_id, has_ec, is_ec7, ec_uniprot, gene_name,
                                accession),
            by = "gene_id") %>%
  group_by(rxn_id) %>%
  summarise(recoverable = any(has_ec, na.rm = TRUE),
            ec7         = any(is_ec7, na.rm = TRUE),
            genes_ec    = paste(sort(unique(gene_id[has_ec])), collapse = "; "),
            names_ec    = paste(sort(unique(na.omit(gene_name[has_ec]))), collapse = "; "),
            accs        = paste(sort(unique(na.omit(accession[has_ec]))), collapse = "; "),
            ecs         = paste(sort(unique(na.omit(ec_uniprot[has_ec]))), collapse = "; "),
            .groups = "drop")

up <- rec %>% filter(recoverable)
cat(sprintf("\nCategory-3 reactions upgradable to Category 1: %d / %d (%.1f%%)\n",
            nrow(up), nrow(rec), 100 * nrow(up) / nrow(rec)))
cat(sprintf("  reactions whose genes carry an EC 7 assignment: %d\n", sum(up$ec7)))

## do the 79 entries collapse? does the upgrade fix all compartments?
sub <- m %>% filter(ID %in% up$rxn_id)
cat(sprintf("  corresponding to %d distinct annotated transformations and %d distinct stoichiometries\n",
            n_distinct(sub$annot_key), n_distinct(sub$eq_nocomp)))

tot_per_eq <- m   %>% count(eq_nocomp, name = "n_total")
rec_per_eq <- sub %>% count(eq_nocomp, name = "n_recoverable")
complete <- rec_per_eq %>% left_join(tot_per_eq, by = "eq_nocomp") %>%
  mutate(complete = n_recoverable == n_total)
cat(sprintf("  upgrade complete for %d of %d stoichiometries; partial for %d\n",
            sum(complete$complete), nrow(complete), sum(!complete$complete)))

## ---- Supplementary Table S3 ------------------------------------------------
S2 <- sub %>%
  left_join(up, by = c("ID" = "rxn_id")) %>%
  left_join(complete %>% select(eq_nocomp, n_total, n_recoverable, complete),
            by = "eq_nocomp") %>%
  transmute(`Reaction ID` = ID, `Reaction name` = NAME, Equation = EQUATION,
            Subsystem = SUBSYSTEM,
            `Compartment(s)` = map_chr(comp_set, ~ paste(comp_map[.x], collapse = "; ")),
            `Model gene association` = GPR,
            `Gene(s) with EC in UniProtKB` = genes_ec,
            `Standard gene name` = names_ec,
            `UniProt accession` = accs,
            `EC number from UniProtKB` = ecs,
            `EC 7 (translocase)` = ifelse(ec7, "yes", "no"),
            `Native confidence score` = native_score,
            `Representations in model` = n_total,
            `Representations recoverable` = n_recoverable,
            `Upgrade complete` = ifelse(complete, "yes", "no")) %>%
  arrange(`Reaction ID`)

write_out(S2, "TableS3_category3_recoverable")
write_out(c3_genes, "table_category3_genes_uniprot")

message("\n05_category3_recovery.R done.")
