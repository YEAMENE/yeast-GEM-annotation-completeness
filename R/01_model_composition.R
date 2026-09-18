## ============================================================================
## 01_model_composition.R — Section 3.1
##
## Model composition, the three counting units, compartmental distribution.
## Produces Figure 1 (A: reactions per compartment split into exclusive and
## shared; B: number of compartments per stoichiometry) and Supplementary
## Table S2 (stoichiometries occurring in more than one compartment).
## ============================================================================

source("R/00_config.R")
m <- load_model()
comp_map <- attr(m, "comp_map")

banner("Model composition")
cat(sprintf("Reaction entries : %d\n", nrow(m)))

n_entries <- nrow(m)
n_annot   <- n_distinct(m$annot_key)
n_stoich  <- n_distinct(m$eq_nocomp)

cat(sprintf("Distinct annotated transformations : %d  (%d true duplicates removed)\n",
            n_annot, n_entries - n_annot))
cat(sprintf("Distinct stoichiometries           : %d\n", n_stoich))

## ---- how the two deduplication levels differ -------------------------------
rep_n <- m %>% count(eq_nocomp, name = "n_rep")
once  <- sum(rep_n$n_rep == 1); many <- sum(rep_n$n_rep > 1)
cat(sprintf("  represented once: %d   more than once: %d\n", once, many))

multi <- m %>% semi_join(filter(rep_n, n_rep > 1), by = "eq_nocomp")
disc <- multi %>% group_by(eq_nocomp) %>%
  summarise(gene_diff = n_distinct(gene_set) > 1,
            ec_diff   = n_distinct(ec_set)   > 1,
            cat_diff  = n_distinct(category) > 1, .groups = "drop")
cat(sprintf("  of these: gene differs %d | EC differs %d | category differs %d | identical %d\n",
            sum(disc$gene_diff), sum(disc$ec_diff), sum(disc$cat_diff),
            sum(!disc$gene_diff & !disc$ec_diff)))

## reversibility is part of the stoichiometric identity
rev_split <- m %>%
  mutate(arrow  = ifelse(str_detect(EQUATION, "<=>"), "<=>", "=>"),
         no_arr = str_replace(eq_nocomp, "<=>", "=>")) %>%
  group_by(no_arr) %>% summarise(n = n_distinct(arrow), .groups = "drop")
cat(sprintf("  stoichiometries with both a reversible and an irreversible form: %d\n",
            sum(rev_split$n > 1)))

## ---- intra vs inter compartmental ------------------------------------------
banner("Compartmental span")
cat(sprintf("Maximum compartments per entry: %d\n", max(m$n_comps)))
print(m %>% count(span) %>% mutate(pct = round(100 * n / sum(n), 1)))

inter <- m %>% filter(n_comps == 2)
pairs <- inter %>% count(comp_pair, sort = TRUE)
cat(sprintf("Distinct compartment pairs observed: %d of %d possible\n",
            nrow(pairs), choose(length(comp_map), 2)))
cat(sprintf("Pairs involving the cytoplasm: %d (%.1f%%)\n",
            sum(inter$n_comps == 2 & map_lgl(inter$comp_set, ~ "c" %in% .x)),
            100 * mean(map_lgl(inter$comp_set, ~ "c" %in% .x))))
print(head(pairs, 3))

## stoichiometries occurring in both forms
form <- m %>% group_by(eq_nocomp) %>%
  summarise(intra = any(n_comps == 1), inter = any(n_comps == 2), .groups = "drop")
cat(sprintf("\nOf %d stoichiometries: %d exclusively intra, %d exclusively inter, %d both\n",
            nrow(form), sum(form$intra & !form$inter),
            sum(!form$intra & form$inter), sum(form$intra & form$inter)))

both <- m %>% semi_join(filter(form, intra & inter), by = "eq_nocomp") %>%
  arrange(eq_nocomp, n_comps) %>%
  select(ID, NAME, EQUATION, n_comps, comp_pair)
print(as.data.frame(both))
write_out(both, "table_stoichiometries_both_forms")

## ---- intra-compartmental detail --------------------------------------------
single <- m %>% filter(n_comps == 1)
occ <- single %>% distinct(eq_nocomp, compartment) %>%
  count(eq_nocomp, name = "n_comp_for_rxn")
single <- single %>% left_join(occ, by = "eq_nocomp")

n_intra <- nrow(single)
n_uniq  <- n_distinct(single$eq_nocomp)

banner("Intra-compartmental reactions")
cat(sprintf("%d entries corresponding to %d distinct stoichiometries\n", n_intra, n_uniq))
cat(sprintf("  exclusive to one compartment: %d (%.1f%%)\n",
            sum(occ$n_comp_for_rxn == 1), 100 * mean(occ$n_comp_for_rxn == 1)))
cat(sprintf("  shared across 2+            : %d (%.1f%%)\n",
            sum(occ$n_comp_for_rxn > 1), 100 * mean(occ$n_comp_for_rxn > 1)))

## ---- Figure 1A -------------------------------------------------------------
byc <- single %>%
  mutate(Status = ifelse(n_comp_for_rxn == 1, "Exclusive to one compartment",
                                              "Shared across \u22652 compartments")) %>%
  count(compartment, Status) %>%
  group_by(compartment) %>% mutate(tot = sum(n)) %>% ungroup() %>%
  mutate(Label  = paste0(comp_map[compartment], " (", compartment, ")"),
         Status = factor(Status, levels = c("Exclusive to one compartment",
                                            "Shared across \u22652 compartments")))

tot_lab <- byc %>% distinct(Label, tot) %>%
  mutate(pct = 100 * tot / n_intra, lab = sprintf("%d (%.1f%%)", tot, pct))
ord <- tot_lab %>% arrange(tot) %>% pull(Label)
byc     <- byc     %>% mutate(Label = factor(Label, levels = ord))
tot_lab <- tot_lab %>% mutate(Label = factor(Label, levels = ord))

write_out(tot_lab %>% arrange(desc(tot)), "table_reactions_per_compartment")

pA <- ggplot(byc, aes(x = n, y = Label, fill = Status)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = tot_lab, aes(x = tot, y = Label, label = lab),
            inherit.aes = FALSE, hjust = -0.10, size = 3.1, colour = "grey25") +
  scale_fill_manual(values = setNames(c(COL$primary, COL$secondary),
                                      levels(byc$Status)), name = NULL) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Number of reactions", y = NULL, tag = "A") +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(size = 9.5, colour = "black"),
        legend.position = "bottom")

## ---- Figure 1B -------------------------------------------------------------
dist <- occ %>% count(n_comp_for_rxn, name = "n_rxn") %>%
  mutate(pct = 100 * n_rxn / n_uniq,
         lab = sprintf("%s (%.1f%%)", comma(n_rxn), pct),
         x   = factor(n_comp_for_rxn))
write_out(dist, "table_compartment_spread")

pB <- ggplot(dist, aes(x = x, y = n_rxn)) +
  geom_col(width = 0.62, fill = COL$primary) +
  geom_text(aes(label = lab), vjust = -0.55, size = 3.1, colour = "grey25") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.16))) +
  labs(x = "Number of compartments in which the stoichiometry occurs",
       y = "Number of stoichiometries", tag = "B") +
  theme_paper() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))

save_fig((pA / pB) + plot_layout(heights = c(1.55, 1)),
         "Figure1_compartmental_distribution", 9, 9.5)

## the maximum-spread stoichiometries are worth reporting explicitly
maxn <- max(occ$n_comp_for_rxn)
top <- single %>% filter(n_comp_for_rxn == maxn) %>%
  group_by(eq_nocomp) %>%
  summarise(compartments = paste(sort(comp_map[unique(compartment)]), collapse = "; "),
            reaction_ids = paste(sort(ID), collapse = "; "),
            genes = paste(sort(unique(GPR[!is.na(GPR)])), collapse = " | "),
            .groups = "drop")
cat(sprintf("\nStoichiometries in %d compartments:\n", maxn))
print(as.data.frame(top))
write_out(top, "table_max_compartment_spread")

## ---- Supplementary Table S2 ------------------------------------------------
shared_eq <- occ %>% filter(n_comp_for_rxn > 1) %>% pull(eq_nocomp)
sub <- single %>% filter(eq_nocomp %in% shared_eq)

grp_info <- sub %>% group_by(eq_nocomp) %>%
  summarise(n_gene_sets = n_distinct(gene_set),
            no_gene     = all(is.na(GPR)),
            ec_differs  = n_distinct(ec_set) > 1,
            cats        = paste(sort(unique(cat_num)), collapse = "/"),
            .groups = "drop") %>%
  mutate(recurrence = case_when(no_gene          ~ "no gene association",
                                n_gene_sets > 1  ~ "isoenzymes (different genes)",
                                TRUE             ~ "multi-localised (same gene)"),
         group = row_number())

cat("\nRecurrence type of the shared stoichiometries:\n")
print(count(grp_info, recurrence))
cat(sprintf("EC assignment differs between compartments in %d groups\n",
            sum(grp_info$ec_differs)))
cat(sprintf("Entries involved: %d, of which %d in Category 1 (%.1f%%)\n",
            nrow(sub), sum(sub$cat_num == 1), 100 * mean(sub$cat_num == 1)))

S4 <- sub %>%
  left_join(grp_info %>% select(eq_nocomp, group, recurrence, ec_differs, cats),
            by = "eq_nocomp") %>%
  transmute(Group = group,
            `Stoichiometry (compartment tags removed)` = eq_nocomp,
            `Reaction ID` = ID, `Reaction name` = NAME,
            Compartment = comp_map[compartment],
            `EC number` = EC, `Gene association` = GPR,
            Category = category, `Categories in group` = cats,
            `Recurrence type` = recurrence,
            `EC differs between compartments` = ifelse(ec_differs, "yes", "no"),
            `Native confidence score` = native_score) %>%
  arrange(Group, `Reaction ID`)
write_out(S4, "TableS2_shared_stoichiometries")

message("\n01_model_composition.R done.")
