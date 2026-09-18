## ============================================================================
## 09_compartment_validation.R — Section 3.4
##
## Reaction compartment assignments are cross-checked against Gene Ontology
## Cellular Component evidence from UniProtKB and SGD.
## Produces Figure 7, Supplementary Table S5, and a breakdown of the "No
## external evidence" reactions into IEA-only versus no-gene-association
## (table_no_evidence_IEA_breakdown.csv), reported in Section 3.4.
##
## Evidence codes: all manually assigned codes are retained, i.e. every code
## except IEA (electronically inferred annotations, which are not reviewed by a
## curator). Note that this set includes computational inferences such as ISS
## and IBA, which are reviewed but not experimental.
##
## The generic "cell envelope" compartment of the model has no single
## corresponding GO term (it spans plasma membrane, cell wall and periplasm) and
## is therefore reported as a separate outcome rather than scored.
##
## Requires: uniprot_yeast_reviewed.tsv and the SGD GO annotation file (.gaf).
## ============================================================================

source("R/00_config.R")

m        <- load_model()
comp_map <- attr(m, "comp_map")
genes    <- read_excel(MODEL_XLSX, sheet = "GENES")

if (!file.exists(UNIPROT_TSV))
  stop("UniProtKB export not found. Run 05_category3_recovery.R first.")
if (is.na(SGD_GAF))
  stop("SGD GO annotation file (.gaf) not found in ", DATA_DIR,
       "\nSee data/README.md.")

## ---- compartment -> GO Cellular Component ----------------------------------
COMPARTMENT_TO_GO <- list(
  c   = c("GO:0005737", "GO:0005829"),                             # cytoplasm, cytosol
  e   = c("GO:0005576", "GO:0005615"),                             # extracellular
  er  = c("GO:0005783", "GO:0005788"),
  erm = c("GO:0005789", "GO:0030176"),
  g   = c("GO:0005794", "GO:0005796"),
  gm  = c("GO:0000139"),
  lp  = c("GO:0005811"),
  m   = c("GO:0005739", "GO:0005759"),
  mm  = c("GO:0031966", "GO:0005743", "GO:0005741", "GO:0005758"),
  n   = c("GO:0005634", "GO:0005654"),
  p   = c("GO:0005777", "GO:0005778"),
  v   = c("GO:0005773", "GO:0000324", "GO:0000323"),
  vm  = c("GO:0005774", "GO:0000329"),
  ce  = character(0)      # no single corresponding GO term
)

MANUALLY_ASSIGNED <- c("EXP","IDA","IPI","IMP","IGI","IEP",     # experimental
                       "HTP","HDA","HMP","HGI","HEP",           # high-throughput
                       "ISS","ISA","ISM","ISO","IBA","RCA",     # computational, reviewed
                       "TAS","NAS","IC")                        # author / curator

## ---- external evidence -----------------------------------------------------
uniprot <- read_tsv(UNIPROT_TSV, show_col_types = FALSE) %>%
  distinct(Entry, .keep_all = TRUE)
up_go <- setNames(uniprot$`Gene Ontology (cellular component)`, uniprot$Entry)
gene_to_acc <- setNames(str_extract(genes$MIRIAM, "(?<=uniprot/)[A-Za-z0-9]+"),
                        genes$NAME)

GAF_COLS <- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID",
              "DB_Reference","Evidence_Code","With_From","Aspect","DB_Object_Name",
              "DB_Object_Synonym","DB_Object_Type","Taxon","Date","Assigned_By",
              "Annotation_Extension","Gene_Product_Form_ID")

gaf <- read_tsv(SGD_GAF, comment = "!", col_names = GAF_COLS,
                col_types = cols(.default = "c"), progress = FALSE)

## SGD indexes annotations by standard gene symbol, so the systematic ORF name
## must be recovered from the symbol or the synonym field.
sgd_go <- gaf %>%
  filter(Aspect == "C", Evidence_Code %in% MANUALLY_ASSIGNED) %>%
  mutate(all_names = paste(DB_Object_Symbol, DB_Object_Synonym, sep = "|")) %>%
  select(all_names, GO_ID) %>%
  separate_rows(all_names, sep = "\\|") %>%
  filter(all_names %in% genes$NAME) %>%
  distinct(all_names, GO_ID) %>%
  group_by(all_names) %>% summarise(go = list(GO_ID), .groups = "drop")
sgd_lookup <- setNames(sgd_go$go, sgd_go$all_names)

banner("Compartmental validation")
cat(sprintf("Model genes with SGD GO-CC evidence: %d of %d\n",
            sum(genes$NAME %in% names(sgd_lookup)), nrow(genes)))

## ---- consistency check -----------------------------------------------------
check <- function(gpr, compartment) {
  expected <- COMPARTMENT_TO_GO[[compartment]]
  if (is.null(expected) || length(expected) == 0)
    return(list(status = "Not checkable (cell envelope)", src = "", obs = ""))

  gs <- extract_genes(gpr)
  found <- character(0); src <- character(0)
  for (g in gs) {
    acc <- gene_to_acc[g]
    if (!is.na(acc) && acc %in% names(up_go)) {
      terms <- str_extract_all(up_go[[acc]], "GO:\\d+")[[1]]
      if (length(terms)) { found <- c(found, terms); src <- c(src, "UniProt") }
    }
    if (!is.null(sgd_lookup[[g]])) {
      found <- c(found, sgd_lookup[[g]]); src <- c(src, "SGD")
    }
  }
  if (!length(found))
    return(list(status = "No external evidence", src = "", obs = ""))
  list(status = if (any(expected %in% found)) "Consistent" else "Inconsistent",
       src = paste(unique(src), collapse = ","),
       obs = paste(head(unique(found), 12), collapse = ";"))
}

single <- m %>% filter(n_comps == 1)
res <- map2(single$GPR, single$compartment, check)

single <- single %>%
  mutate(validation = map_chr(res, "status"),
         sources    = map_chr(res, "src"),
         go_found   = map_chr(res, "obs"),
         validation = factor(validation,
                             levels = c("Consistent", "Inconsistent",
                                        "No external evidence",
                                        "Not checkable (cell envelope)")))

vc <- single %>% count(validation) %>% mutate(pct = round(100 * n / sum(n), 1))
print(as.data.frame(vc))
write_out(single %>% transmute(ID, NAME, compartment,
                               Compartment_name = comp_map[compartment],
                               EC, GPR, category, native_score,
                               validation, sources, go_found),
          "compartment_validation_full")

## ---- do our scores anticipate the inconsistencies? -------------------------
banner("Validation status versus the two scoring schemes")
tab_cat <- single %>% count(validation, category) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)
print(as.data.frame(tab_cat)); write_out(tab_cat, "table_validation_x_category")

tab_sc <- single %>%
  mutate(native_score = factor(native_score, levels = SCORE_LEVELS)) %>%
  count(validation, native_score) %>%
  pivot_wider(names_from = native_score, values_from = n, values_fill = 0)
print(as.data.frame(tab_sc)); write_out(tab_sc, "table_validation_x_native_score")

inc <- single %>% filter(validation == "Inconsistent")
con <- single %>% filter(validation == "Consistent")
cat(sprintf("\nInconsistent: %d, of which %d (%.1f%%) in Category 1 and %d (%.1f%%) with native score 3\n",
            nrow(inc), sum(inc$cat_num == 1), 100 * mean(inc$cat_num == 1),
            sum(inc$native_score == "3"), 100 * mean(inc$native_score == "3")))
cat(sprintf("Consistent  : %d, of which %.1f%% in Category 1 and %.1f%% with native score 3\n",
            nrow(con), 100 * mean(con$cat_num == 1), 100 * mean(con$native_score == "3")))
cat("Compartmental misassignment is therefore orthogonal to both scores and is\n")
cat("detectable only through explicit comparison with external localisation data.\n")

## ---- per-compartment rate --------------------------------------------------
by_comp <- single %>% count(compartment, validation) %>%
  group_by(compartment) %>% mutate(total = sum(n)) %>% ungroup() %>%
  mutate(Label = paste0(comp_map[compartment], " (", compartment, ")"))

rate <- by_comp %>% group_by(Label, total) %>%
  summarise(inc = sum(n[validation == "Inconsistent"]), .groups = "drop") %>%
  mutate(pct_inc = 100 * inc / total)
print(as.data.frame(rate %>% arrange(desc(pct_inc))), digits = 3)

ord     <- rate %>% arrange(total) %>% pull(Label)
by_comp <- by_comp %>% mutate(Label = factor(Label, levels = ord))
rate    <- rate    %>% mutate(Label = factor(Label, levels = ord))

## ---- Figure 7 --------------------------------------------------------------
VCOL <- c("Consistent"                    = "#2C7BB6",
          "Inconsistent"                  = "#D7191C",
          "No external evidence"          = "#D9D9D9",
          "Not checkable (cell envelope)" = "#7B7B7B")

pA <- ggplot(by_comp, aes(x = n, y = Label, fill = validation)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = rate, aes(x = total, y = Label, label = comma(total)),
            inherit.aes = FALSE, hjust = -0.18, size = 3.1, colour = "grey25") +
  scale_fill_manual(values = VCOL, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.13))) +
  labs(x = "Number of reactions", y = NULL, tag = "A") +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(size = 9.5, colour = "black"),
        legend.position = "bottom")

pB <- ggplot(rate, aes(x = pct_inc, y = Label)) +
  geom_col(width = 0.68, fill = VCOL[["Inconsistent"]]) +
  geom_text(aes(label = ifelse(inc > 0, sprintf("%d  (%.1f%%)", inc, pct_inc), "\u2013")),
            hjust = -0.15, size = 3.0, colour = "grey25") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     breaks = c(0, 5, 10, 15),
                     limits = c(0, max(rate$pct_inc) * 2.1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Inconsistent reactions", y = NULL, tag = "B") +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

save_fig((pA | pB) + plot_layout(widths = c(2.4, 1), guides = "collect") &
           theme(legend.position = "bottom"),
         "Figure7_compartment_validation", 11, 5.6)

## ---- Supplementary Table S5 ------------------------------------------------
S5 <- single %>%
  filter(validation != "Consistent") %>%
  transmute(`Reaction ID` = ID, `Reaction name` = NAME, Equation = EQUATION,
            `Model compartment` = comp_map[compartment],
            `Gene association` = GPR, `EC number` = EC,
            `EC/gene category` = category,
            `Native confidence score` = native_score,
            `Validation status` = as.character(validation),
            `External evidence source` = sources,
            `GO terms found (truncated)` = go_found) %>%
  arrange(`Validation status`, `Reaction ID`)

write_out(S5, "TableS5_compartment_validation")
cat(sprintf("\nSupplementary Table S5: %d reactions not consistently annotated\n", nrow(S5)))

## ============================================================================
## IEA-only vs no-annotation breakdown of the "No external evidence" reactions
## (Section 3.4: distinguishes a genuine annotation gap from the deliberate
## exclusion of electronically inferred, IEA-coded evidence, Methods 2.7)
##
## If a reaction was classified "no_evidence", neither UniProtKB nor SGD with
## manually assigned codes produced a GO-CC term. Re-checking the SAME SGD GAF
## WITHOUT the evidence-code filter: any term found for that gene can only be
## IEA-coded, since a manually assigned one would already have been found above.
## ============================================================================
sgd_cc_any <- gaf %>%
  filter(Aspect == "C") %>%
  mutate(all_names = paste(DB_Object_Symbol, DB_Object_Synonym, sep = "|")) %>%
  select(all_names, GO_ID) %>%
  separate_rows(all_names, sep = "\\|") %>%
  filter(all_names %in% gene_ids) %>%
  distinct(all_names, GO_ID)

genes_with_any_cc <- unique(sgd_cc_any$all_names)

classify_no_evidence <- function(gpr) {
  gs <- extract_genes(gpr)
  if (length(gs) == 0) return("no_gene_association")
  if (any(gs %in% genes_with_any_cc)) return("iea_only")
  "no_ccm_at_all"
}

noev <- single %>%
  filter(validation == "no_evidence") %>%
  mutate(iea_breakdown = map_chr(GPR, classify_no_evidence))

banner("Breakdown of the 'No external evidence' reactions (Section 2.7 / 3.4)")
print(noev %>% count(iea_breakdown) %>% mutate(pct = round(100 * n / sum(n), 1)))

n_iea    <- sum(noev$iea_breakdown == "iea_only")
n_none   <- sum(noev$iea_breakdown == "no_ccm_at_all")
n_nogene <- sum(noev$iea_breakdown == "no_gene_association")

cat(sprintf("  IEA-coded GO-CC term present (excluded by criteria, Methods 2.7): %d (%.1f%%)\n",
            n_iea, 100 * n_iea / nrow(noev)))
cat(sprintf("  No GO Cellular Component annotation at all, any source/evidence: %d (%.1f%%)\n",
            n_none, 100 * n_none / nrow(noev)))
cat(sprintf("  No gene association (cannot be evaluated in principle): %d (%.1f%%)\n",
            n_nogene, 100 * n_nogene / nrow(noev)))

write_out(noev %>% select(ID, NAME, compartment, GPR, iea_breakdown),
          "table_no_evidence_IEA_breakdown")

message("\n09_compartment_validation.R done.")
