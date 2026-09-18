## ============================================================================
## 06_genome_coverage.R — Section 3.3
##
## Distribution of the model gene set across the yeast genome.
## Produces Figure 5: (A) model genes per chromosome; (B) total protein-coding
## gene complement of each chromosome with the covered fraction embedded.
##
## Chromosome assignment uses the systematic ORF nomenclature: the second
## character of the locus identifier encodes the chromosome (A-P for I-XVI),
## while mitochondrially-encoded loci carry a "Q" prefix (Methods 2.3).
##
## Requires the SGD reference genome annotation in GFF3 format; download from
## http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/
## ============================================================================

source("R/00_config.R")

genes <- read_excel(MODEL_XLSX, sheet = "GENES")

## ---- model side ------------------------------------------------------------
model_genes <- genes %>%
  transmute(gene_id = NAME, chr = chromosome_of(NAME))

banner("Genome coverage by chromosome")
cat(sprintf("Model genes: %d (unassigned: %d)\n",
            nrow(model_genes), sum(is.na(model_genes$chr))))

model_chr <- model_genes %>% filter(!is.na(chr)) %>% count(chr, name = "n_model")

## ---- reference genome side -------------------------------------------------
if (is.na(SGD_GFF))
  stop("SGD GFF3 annotation not found in ", DATA_DIR,
       "\nSee data/README.md for the download link.")

gff <- read_tsv(SGD_GFF, comment = "#", col_names = FALSE,
                col_types = cols(.default = "c"), progress = FALSE)
names(gff)[1:9] <- c("seqid","source","type","start","end",
                     "score","strand","phase","attributes")

## SGD uses separate feature types for non-coding RNA genes (ncRNA_gene,
## tRNA_gene, rRNA_gene, snoRNA_gene), so filtering type == "gene" already
## isolates protein-coding ORFs.
genome_genes <- gff %>%
  filter(type == "gene") %>%
  mutate(chr     = ifelse(seqid %in% c("chrmt", "chrMT", "chrM"),
                          "Mitochondrial", str_remove(seqid, "^chr")),
         gene_id = str_match(attributes, "ID=([^;]+)")[, 2])

genome_chr <- genome_genes %>% count(chr, name = "n_genome")

cov <- genome_chr %>%
  left_join(model_chr, by = "chr") %>%
  mutate(n_model = replace_na(n_model, 0L),
         n_uncov = n_genome - n_model,
         pct     = 100 * n_model / n_genome,
         chr     = factor(chr, levels = CHR_LEVELS)) %>%
  filter(!is.na(chr)) %>% arrange(chr)

print(as.data.frame(cov %>% select(chr, n_model, n_genome, pct)), digits = 3)
write_out(cov, "table_genome_coverage")

nuc <- cov %>% filter(chr != "Mitochondrial")
cat(sprintf("\nNuclear: %d of %d protein-coding genes = %.1f%%\n",
            sum(nuc$n_model), sum(nuc$n_genome),
            100 * sum(nuc$n_model) / sum(nuc$n_genome)))
cat(sprintf("Per-chromosome mean: %.1f%% +/- %.1f%%   (range %.1f%% - %.1f%%)\n",
            mean(nuc$pct), sd(nuc$pct), min(nuc$pct), max(nuc$pct)))
mt <- cov %>% filter(chr == "Mitochondrial")
if (nrow(mt) > 0)
  cat(sprintf("Mitochondrial: %d of %d = %.1f%%\n",
              mt$n_model, mt$n_genome, mt$pct))

cat("\nNote: this proportion is not a measure of model incompleteness. The\n")
cat("denominator is the entire protein-coding complement, which includes\n")
cat("ribosomal, transcriptional, cytoskeletal and DNA-repair functions that a\n")
cat("metabolic reconstruction is not expected to represent.\n")

## ---- Figure 5A: model genes per chromosome ---------------------------------
pA <- ggplot(cov, aes(x = chr, y = n_model)) +
  geom_col(fill = COL$primary, width = 0.72) +
  geom_text(aes(label = comma(n_model)), vjust = -0.40, size = 3.0,
            colour = "grey25") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Model genes", tag = "A") +
  theme_paper() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))

## ---- Figure 5B: total complement with the covered fraction embedded --------
covB <- cov %>%
  select(chr, Covered = n_model, `Not covered` = n_uncov) %>%
  pivot_longer(-chr, names_to = "status", values_to = "n") %>%
  mutate(status = factor(status, levels = c("Covered", "Not covered")))

pB <- ggplot(covB, aes(x = chr, y = n, fill = status)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = cov, aes(x = chr, y = n_genome, label = sprintf("%.1f%%", pct)),
            inherit.aes = FALSE, vjust = -0.45, size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c(Covered = COL$primary, `Not covered` = COL$light),
                    name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Chromosome", y = "Protein-coding genes", tag = "B") +
  theme_paper() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        legend.position = "bottom")

save_fig((pA / pB) + plot_layout(heights = c(1, 1.15)),
         "Figure5_genome_coverage", 10, 8)

message("\n06_genome_coverage.R done.")
