## ============================================================================
## 02_annotation_categories.R — Section 3.2
##
## EC/gene confidence categories, EC specificity sub-categories (1a/1b, 2a/2b),
## and the breakdown by number of compartmental representations.
## Produces Figure 2 (categories on entries and on distinct annotated
## transformations) and the full reaction dataset deposited in the repository.
## ============================================================================

source("R/00_config.R")
m <- load_model()
comp_map <- attr(m, "comp_map")

n_entries <- nrow(m)
uniq      <- m %>% distinct(annot_key, .keep_all = TRUE)
n_uniq    <- nrow(uniq)

## ---- categories ------------------------------------------------------------
banner("EC/gene confidence categories")

cat_tab <- bind_rows(
  m    %>% count(category) %>% mutate(unit = "Model entries",                      total = n_entries),
  uniq %>% count(category) %>% mutate(unit = "Distinct annotated transformations", total = n_uniq)
) %>%
  mutate(pct = round(100 * n / total, 1),
         category = factor(category, levels = CAT_LEVELS),
         unit = factor(unit, levels = c("Model entries",
                                        "Distinct annotated transformations")))

print(cat_tab %>% select(unit, category, n, pct) %>%
        pivot_wider(names_from = unit, values_from = c(n, pct)) %>% as.data.frame())
write_out(cat_tab, "table_categories")

## sanity check: on the 3,701 key every transformation has exactly one category
amb <- uniq %>% count(annot_key) %>% filter(n > 1)
stopifnot(nrow(amb) == 0)
amb_stoich <- m %>% group_by(eq_nocomp) %>% summarise(k = n_distinct(category)) %>%
  filter(k > 1)
cat(sprintf("\nStoichiometries with an ambiguous category: %d\n", nrow(amb_stoich)))
cat("  (this is why the 3,398 stoichiometries cannot be used as a denominator)\n")

## ---- categories by number of representations -------------------------------
rep_n <- m %>% count(eq_nocomp, name = "n_rep")
m2 <- m %>% left_join(rep_n, by = "eq_nocomp")

single_rep <- m2 %>% filter(n_rep == 1)
multi_rep  <- m2 %>% filter(n_rep >  1)

banner("Categories by number of compartmental representations")
cat(sprintf("Stoichiometries represented once: %d\n", n_distinct(single_rep$eq_nocomp)))
print(single_rep %>% count(category) %>%
        mutate(pct = round(100 * n / sum(n), 1)))

grp <- multi_rep %>% group_by(eq_nocomp) %>%
  summarise(cats = paste(sort(unique(cat_num)), collapse = "/"), .groups = "drop")
cat(sprintf("\nStoichiometries represented more than once: %d\n", nrow(grp)))
cat(sprintf("  consistent category: %d | different category: %d (%.1f%%)\n",
            sum(!str_detect(grp$cats, "/")), sum(str_detect(grp$cats, "/")),
            100 * mean(str_detect(grp$cats, "/"))))
print(count(grp, cats, sort = TRUE))

## ---- EC specificity: sub-categories a / b ----------------------------------
banner("EC specificity (sub-categories a / b)")
spec <- m %>% filter(cat_num %in% c(1, 2)) %>% count(category, sub_category) %>%
  group_by(category) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
print(as.data.frame(spec))
write_out(spec, "table_ec_specificity")

ec_bearing <- m %>% filter(!is.na(EC))
partial    <- ec_bearing %>% filter(ec_spec == "partial only")
cat(sprintf("\nEC-bearing reactions: %d, of which %d (%.1f%%) carry only partial ECs\n",
            nrow(ec_bearing), nrow(partial), 100 * nrow(partial) / nrow(ec_bearing)))

## truncation depth: how much of the EC hierarchy is actually resolved
depth <- function(e) sum(str_split(e, "\\.")[[1]] != "-")
class_only <- partial %>%
  mutate(maxd = map_int(ec_list, ~ max(map_int(.x, depth)))) %>%
  filter(maxd == 1)
cat(sprintf("  truncated at the class level alone: %d\n", nrow(class_only)))
only_ec2 <- partial %>%
  filter(map_lgl(ec_list, ~ all(.x == "2.-.-.-")))
n_ec2 <- m %>% filter(map_lgl(ec_list, ~ any(str_starts(.x, "2")))) %>% nrow()
cat(sprintf("  annotated only as \"2.-.-.-\": %d (%.1f%% of all EC-2 reactions)\n",
            nrow(only_ec2), 100 * nrow(only_ec2) / n_ec2))

## partial ECs are invisible to the native evidence score
cat1b <- m %>% filter(sub_category == "1b")
cat(sprintf("  of the %d sub-category 1b reactions, %d (%.1f%%) carry native score 3\n",
            nrow(cat1b), sum(cat1b$native_score == "3"),
            100 * mean(cat1b$native_score == "3")))

## ---- Figure 2 --------------------------------------------------------------
plot_df <- cat_tab %>%
  mutate(xlab = factor(str_replace(category, " \\(", "\n("),
                       levels = str_replace(CAT_LEVELS, " \\(", "\n(")),
         label = sprintf("%s\n(%.1f%%)", comma(n), pct))

p <- ggplot(plot_df, aes(x = xlab, y = n, fill = unit)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66,
           colour = "white", linewidth = 0.25) +
  geom_text(aes(label = label), position = position_dodge(width = 0.74),
            vjust = -0.30, size = 3.0, lineheight = 0.95, colour = "grey20") +
  scale_fill_manual(values = setNames(c(COL$primary, COL$secondary),
                                      levels(cat_tab$unit)), name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Number of reactions") +
  theme_paper(12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.x = element_text(size = 10, lineheight = 1.05, colour = "black"),
        legend.position = "bottom")

save_fig(p, "Figure2_EC_gene_categories", 8.5, 5.6)

## ---- full reaction dataset (deposited in the repository) -------------------
grp_id <- setNames(seq_along(unique(m$eq_nocomp)),  sort(unique(m$eq_nocomp)))
ann_id <- setNames(seq_along(unique(m$annot_key)), sort(unique(m$annot_key)))

full <- m %>%
  transmute(`Reaction ID` = ID, `Reaction name` = NAME, Equation = EQUATION,
            Subsystem = SUBSYSTEM,
            `EC number` = EC, `Gene association` = GPR,
            Category = category, `EC specificity` = ec_spec,
            `Sub-category` = sub_category,
            `N compartments` = n_comps,
            `Compartment codes` = comp_pair,
            `Compartment names` = map_chr(comp_set,
                                          ~ paste(comp_map[.x], collapse = "; ")),
            Compartmentalisation = span,
            `Labelled as transport` = ifelse(
              str_starts(replace_na(SUBSYSTEM, ""), "Transport"), "yes", "no"),
            `Native confidence score` = native_score,
            `Stoichiometry group` = grp_id[eq_nocomp],
            `Annotated transformation group` = ann_id[annot_key]) %>%
  group_by(`Stoichiometry group`) %>%
  mutate(`Entries in stoichiometry group` = n()) %>% ungroup() %>%
  arrange(`Reaction ID`)

write_out(full, "full_reaction_dataset")
cat(sprintf("\nFull dataset written: %d rows\n", nrow(full)))

message("\n02_annotation_categories.R done.")
