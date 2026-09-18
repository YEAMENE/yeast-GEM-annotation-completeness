## ============================================================================
## 04_confidence_crosstab.R — Section 3.2
##
## EC/gene categories versus the native SBML/COBRA confidence score.
## Produces Figure 4 and Supplementary Table S4 (cross-tabulation, the 61
## Category-2 reactions, and the Category-1 reactions carrying no native score).
##
## Note on colour: reactions with no native score are shown in neutral grey,
## outside the sequential blue scale, so that "unscored" is not read as the
## strongest level of evidence.
## ============================================================================

source("R/00_config.R")
m <- load_model()
comp_map <- attr(m, "comp_map")

m <- m %>% mutate(native_score = factor(native_score, levels = SCORE_LEVELS),
                  category     = factor(category,     levels = CAT_LEVELS))

## ---- cross-tabulation ------------------------------------------------------
banner("EC/gene category x native confidence score")

ct <- m %>% count(category, native_score) %>%
  pivot_wider(names_from = native_score, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))
ct_all <- bind_rows(ct, ct %>% summarise(across(where(is.numeric), sum)) %>%
                      mutate(category = "Total", .before = 1))
print(as.data.frame(ct_all))
write_out(ct_all, "TableS4a_category_vs_native_score")

pct <- m %>% count(category, native_score) %>%
  group_by(category) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
write_out(pct, "table_category_vs_native_score_pct")

## the observations reported in the manuscript
c1 <- m %>% filter(cat_num == 1)
cat(sprintf("\nCategory 1: score 3 = %d (%.1f%%), score 2 = %d (%.1f%%)\n",
            sum(c1$native_score == "3"), 100 * mean(c1$native_score == "3"),
            sum(c1$native_score == "2"), 100 * mean(c1$native_score == "2")))
c1u <- c1 %>% filter(native_score == "Unscored")
cat(sprintf("Category 1 unscored: %d (%.1f%%); score 0 or 1: %d\n",
            nrow(c1u), 100 * nrow(c1u) / nrow(c1),
            sum(c1$native_score %in% c("0", "1"))))

## the unscored Category-1 entries cluster in the highest identifier block,
## suggesting they were added in recent releases without a native score
if (nrow(c1u) > 0) {
  num <- as.integer(str_extract(c1u$ID, "[0-9]+"))
  cat(sprintf("  identifier range: r_%04d - r_%04d\n", min(num), max(num)))
  cat("  most frequent subsystems:\n")
  print(head(sort(table(c1u$SUBSYSTEM), decreasing = TRUE), 3))
}

c4 <- m %>% filter(cat_num == 4)
lo <- sum(c4$native_score %in% c("0", "1"))
cat(sprintf("\nCategory 4 with score 0 or 1: %d (%.1f%%); including unscored: %d (%.1f%%)\n",
            lo, 100 * lo / nrow(c4),
            lo + sum(c4$native_score == "Unscored"),
            100 * (lo + sum(c4$native_score == "Unscored")) / nrow(c4)))

## ---- Figure 4 --------------------------------------------------------------
plot_df <- m %>% count(category, native_score) %>%
  group_by(category) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

p <- ggplot(plot_df, aes(x = category, y = pct, fill = native_score)) +
  geom_col(width = 0.66, colour = "black", linewidth = 0.22,
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = COL$score, name = "Native confidence score",
                    drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Percentage of reactions") +
  theme_paper(12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.x = element_text(size = 10, colour = "black"),
        legend.position = "right")

save_fig(p, "Figure4_category_vs_native_score", 8.5, 5)

## ---- Supplementary Table S4b: the 61 Category-2 reactions ------------------
## For each, check whether the same EC is already associated with a gene
## elsewhere in the model, i.e. whether the missing GPR can be recovered
## internally without external evidence.
ec2gene <- m %>% filter(cat_num == 1) %>%
  select(ec_list, GPR) %>% unnest(ec_list) %>%
  group_by(ec_list) %>% summarise(genes = paste(unique(GPR), collapse = "; "),
                                  .groups = "drop")
ec_lookup <- setNames(ec2gene$genes, ec2gene$ec_list)

c2 <- m %>% filter(cat_num == 2) %>%
  mutate(hit_full = map_lgl(ec_list, ~ any(.x %in% names(ec_lookup) & is_full_ec(.x))),
         hit_any  = map_lgl(ec_list, ~ any(.x %in% names(ec_lookup))),
         recoverable = case_when(hit_full ~ "yes (fully specified EC)",
                                 hit_any  ~ "weak (partial EC only)",
                                 TRUE     ~ "no"),
         source_genes = map_chr(ec_list, ~ paste(unique(unlist(
           ec_lookup[.x[.x %in% names(ec_lookup)]])), collapse = "; ")))

banner("Category 2: internal gene recovery")
print(count(c2, recoverable))
orphan <- c2 %>% filter(recoverable == "no")
cat(sprintf("\nOrphan enzymatic activities: %d\n", nrow(orphan)))
top_ec <- sort(table(unlist(orphan$ec_list)), decreasing = TRUE)
cat("Most frequent EC numbers among them:\n"); print(head(top_ec, 3))

S3b <- c2 %>%
  transmute(`Reaction ID` = ID, `Reaction name` = NAME, Equation = EQUATION,
            `Compartment(s)` = map_chr(comp_set, ~ paste(comp_map[.x], collapse = "; ")),
            `EC number` = EC, `EC specificity` = ec_spec,
            `Native confidence score` = as.character(native_score),
            `Gene recoverable from within the model` = recoverable,
            `Gene(s) carrying the same EC elsewhere` = source_genes)
write_out(S3b, "TableS4b_category2_reactions")

## ---- Supplementary Table S4c: Category-1 reactions with no native score ----
S3c <- c1u %>%
  transmute(`Reaction ID` = ID, `Reaction name` = NAME, Equation = EQUATION,
            `Compartment(s)` = map_chr(comp_set, ~ paste(comp_map[.x], collapse = "; ")),
            `EC number` = EC, `Gene association` = GPR, Subsystem = SUBSYSTEM)
write_out(S3c, "TableS4c_category1_unscored")

message("\n04_confidence_crosstab.R done.")
