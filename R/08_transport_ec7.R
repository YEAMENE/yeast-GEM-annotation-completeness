## ============================================================================
## 08_transport_ec7.R — Sections 2.5 and 3.3
##
## Identification of transport reactions and assessment of the EC 7 annotation
## gap.
##
## Two sets are distinguished throughout, because they answer different
## questions and the manuscript uses both:
##   - CANDIDATE transport set: every reaction spanning two compartments.
##     Deliberately over-inclusive; the union of the three criteria of
##     Section 2.6 coincides with the stoichiometric criterion.
##   - ANNOTATED transport set: reactions labelled as transport in the model
##     subsystem field. This is the set used for the EC-coverage statistics.
## ============================================================================

source("R/00_config.R")
m <- load_model()

## ---- the three criteria of Section 2.6 -------------------------------------
tr <- m %>%
  mutate(by_subsystem = str_detect(replace_na(SUBSYSTEM, ""), regex("transport", TRUE)),
         by_name      = str_detect(replace_na(NAME, ""),      regex("transport", TRUE)),
         by_stoich    = n_comps == 2,
         candidate    = by_subsystem | by_name | by_stoich)

banner("Transport reactions")
cat(sprintf("subsystem label : %d\n", sum(tr$by_subsystem)))
cat(sprintf("reaction name   : %d\n", sum(tr$by_name)))
cat(sprintf("two compartments: %d\n", sum(tr$by_stoich)))
cat(sprintf("union           : %d\n", sum(tr$candidate)))

if (sum(tr$candidate) == sum(tr$by_stoich))
  cat("\nNote: every reaction labelled as transport by subsystem or by name also\n",
      "spans two compartments, so the union coincides with the stoichiometric\n",
      "criterion. This should be stated explicitly in Section 2.6.\n", sep = "")

cand <- tr %>% filter(candidate)
ann  <- tr %>% filter(by_subsystem)

cat(sprintf("\nCandidate set: %d entries = %d distinct annotated transformations (%.1f%% of %d)\n",
            nrow(cand), n_distinct(cand$annot_key),
            100 * n_distinct(cand$annot_key) / n_distinct(m$annot_key),
            n_distinct(m$annot_key)))
cat(sprintf("               %d of %d model entries (%.1f%%)\n",
            nrow(cand), nrow(m), 100 * nrow(cand) / nrow(m)))
cat(sprintf("Annotated set: %d entries = %d distinct annotated transformations (%.1f%% of the model)\n",
            nrow(ann), n_distinct(ann$annot_key), 100 * nrow(ann) / nrow(m)))

## what the difference between the two sets contains
diff <- cand %>% filter(!by_subsystem)
cat(sprintf("\nCandidates not labelled as transport: %d\n", nrow(diff)))
print(head(sort(table(replace_na(diff$SUBSYSTEM, "(none)")), decreasing = TRUE), 3))
cat("These span two compartments but do not represent membrane translocation.\n")

## ---- EC annotation of the annotated transport set --------------------------
banner("EC annotation of transport reactions")
with_ec <- ann %>% filter(!is.na(EC))
cat(sprintf("With an EC number: %d of %d (%.1f%%)\n",
            nrow(with_ec), nrow(ann), 100 * nrow(with_ec) / nrow(ann)))
cat(sprintf("With a gene association: %d (%.1f%%)\n",
            sum(!is.na(ann$GPR)), 100 * mean(!is.na(ann$GPR))))

if (nrow(with_ec) > 0) {
  cls <- table(unlist(map(with_ec$ec_list, ~ unique(str_extract(.x, "^[0-9]+")))))
  cat("EC classes represented:\n"); print(cls)
  cat(sprintf("  of class 3 (hydrolases): %d of %d\n",
              as.integer(cls["3"]), nrow(with_ec)))
  cat("\nEC numbers assigned:\n")
  print(sort(table(unlist(with_ec$ec_list)), decreasing = TRUE))
}

cat("\nCategory distribution of the annotated transport set:\n")
print(ann %>% count(category) %>% mutate(pct = round(100 * n / sum(n), 1)))

## ---- EC 7 in the model -----------------------------------------------------
ec7 <- m %>% filter(map_lgl(ec_list, ~ any(str_starts(.x, "7"))))
cat(sprintf("\nReactions annotated with an EC 7 (translocase) number: %d\n", nrow(ec7)))

atp <- m %>% filter(map_lgl(ec_list, ~ any(str_starts(.x, "3.6.3"))))
cat(sprintf("Reactions annotated as ATP hydrolases (EC 3.6.3.-): %d\n", nrow(atp)))
cat("These correspond to the pre-2018 classification of what are now translocases.\n")

write_out(ann %>% transmute(`Reaction ID` = ID, `Reaction name` = NAME,
                            Equation = EQUATION, Subsystem = SUBSYSTEM,
                            `EC number` = EC, `Gene association` = GPR,
                            Category = category,
                            `Native confidence score` = native_score),
          "table_annotated_transport_reactions")

cat("\nThe EC 7 assignments recoverable from UniProtKB are identified in\n")
cat("05_category3_recovery.R (see Supplementary Table S2).\n")

message("\n08_transport_ec7.R done.")
