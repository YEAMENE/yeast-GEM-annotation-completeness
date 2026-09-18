## ============================================================================
## 03_category4_breakdown.R — Section 3.2
##
## Structural composition of Category 4 (no EC number, no gene association).
## Produces Figure 3.
##
## Note on labels: entries whose only model note refers to a repository pull
## request are labelled "curation-provenance note only" rather than "residual
## model note", because those notes carry no biological information.
## ============================================================================

source("R/00_config.R")
m <- load_model()

cat4 <- m %>%
  filter(category == "Cat 4 (neither)") %>%
  mutate(subclass = case_when(
    SUBSYSTEM == "Exchange reaction"                   ~ "Exchange / boundary\nreactions",
    SUBSYSTEM == "SLIME reaction"                      ~ "SLIME lipid\npseudo-reactions",
    str_starts(replace_na(SUBSYSTEM, ""), "Transport") ~ "Transport reactions",
    str_detect(str_to_lower(replace_na(NAME, "")),
               "biomass|maintenance|growth")           ~ "Biomass / maintenance\npseudo-reactions",
    !is.na(NOTE)                                       ~ "Curation-provenance\nnote only",
    TRUE                                               ~ "No model note"))

banner("Category 4 composition")
cat(sprintf("Category 4: %d entries\n\n", nrow(cat4)))

c4 <- cat4 %>% count(subclass) %>%
  mutate(pct = 100 * n / sum(n),
         lab = sprintf("%s (%.1f%%)", comma(n), pct),
         highlight = subclass == "Transport reactions")

print(as.data.frame(c4 %>% select(subclass, n, pct) %>% arrange(desc(n))), digits = 3)
write_out(c4 %>% arrange(desc(n)), "table_category4_breakdown")

tech <- c4 %>% filter(str_detect(subclass, "Transport|Exchange")) %>% pull(pct) %>% sum()
cat(sprintf("\nTransport + exchange/boundary: %.1f%% of the category\n", tech))
cat("These do not represent enzyme-catalysed transformations and are therefore\n")
cat("not expected to carry EC or gene annotation.\n")

## reactions with no informative note (excluding biomass/maintenance)
uninformative <- c4 %>%
  filter(str_detect(subclass, "No model note|Curation-provenance")) %>%
  summarise(n = sum(n), pct = sum(pct))
cat(sprintf("Entries with no informative model note: %d (%.1f%%)\n",
            uninformative$n, uninformative$pct))

## what the "curation-provenance" notes actually contain
prov <- cat4 %>% filter(str_detect(subclass, "Curation-provenance"))
if (nrow(prov) > 0) {
  cat("\nExamples of curation-provenance notes:\n")
  print(head(unique(prov$NOTE), 3))
}

## ---- Figure 3 --------------------------------------------------------------
p <- ggplot(c4, aes(x = n, y = fct_reorder(subclass, n), fill = highlight)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.28) +
  geom_text(aes(label = lab), hjust = -0.08, size = 3.4, colour = "grey20") +
  scale_fill_manual(values = c(`TRUE` = COL$accent, `FALSE` = COL$grey),
                    guide = "none") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Number of reactions", y = NULL) +
  theme_paper(12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(size = 10, lineheight = 0.95, colour = "black"))

save_fig(p, "Figure3_Category4_composition", 8.5, 5)

message("\n03_category4_breakdown.R done.")
