## ============================================================================
## 07_ec_class_distribution.R — Section 3.3
##
## Distribution of top-level EC classes in the model, in the ExplorEnz global
## catalogue and in the reviewed UniProtKB S. cerevisiae proteome.
## Produces Figure 6.
##
## Counting rule (Methods 2.8): each top-level EC class is counted ONCE per
## reaction (model) or per entry (UniProtKB), regardless of how many EC numbers
## are assigned. Partial ECs contribute their class.
## ============================================================================

source("R/00_config.R")
library(httr)

m <- load_model()

EC_NAMES <- c("1" = "Oxidoreductases", "2" = "Transferases", "3" = "Hydrolases",
              "4" = "Lyases", "5" = "Isomerases", "6" = "Ligases",
              "7" = "Translocases")

## unique top-level classes per record
classes_of <- function(ec_lists) {
  unlist(map(ec_lists, ~ unique(str_extract(.x, "^[0-9]+"))))
}

## ---- model -----------------------------------------------------------------
model_counts <- table(classes_of(m$ec_list))

banner("EC class distribution")
cat(sprintf("Model: %d class assignments over %d EC-bearing reactions\n",
            sum(model_counts), sum(!is.na(m$EC))))

## ---- UniProtKB -------------------------------------------------------------
if (!file.exists(UNIPROT_TSV))
  stop("UniProtKB export not found: ", UNIPROT_TSV,
       "\nRun 05_category3_recovery.R first (it caches the file), ",
       "or see data/README.md.")

uniprot <- read_tsv(UNIPROT_TSV, show_col_types = FALSE) %>%
  distinct(Entry, .keep_all = TRUE)
up_lists  <- map(uniprot$`EC number`, split_ec)
up_counts <- table(classes_of(up_lists))

cat(sprintf("UniProtKB: %d reviewed entries, %d with an EC number, %d class assignments\n",
            nrow(uniprot), sum(lengths(up_lists) > 0), sum(up_counts)))

## ---- ExplorEnz -------------------------------------------------------------
## Live counts when reachable, otherwise a dated snapshot. Update the snapshot
## and its date together, and keep them consistent with Supplementary Table S1.
EXPLORENZ_SNAPSHOT      <- c("1"=2041, "2"=2096, "3"=1369, "4"=786,
                             "5"=325,  "6"=257,  "7"=98)   # 6,972 active entries
EXPLORENZ_SNAPSHOT_DATE <- "2026-07"

get_explorenz <- function() {
  out <- tryCatch({
    r <- GET("https://www.enzyme-database.org/stats.php", timeout(20))
    stop_for_status(r)
    txt <- content(r, "text", encoding = "UTF-8")
    hit <- regmatches(txt, regexpr("Current:\\s*[0-9 ,]+", txt))
    if (length(hit) == 0) stop("pattern not found")
    v <- as.integer(str_remove_all(str_extract_all(hit, "[0-9,]+")[[1]], ","))
    if (length(v) < 7) stop("unexpected format")
    message("ExplorEnz: live counts retrieved on ", Sys.Date())
    setNames(v[1:7], as.character(1:7))
  }, error = function(e) {
    message("ExplorEnz live fetch failed (", conditionMessage(e),
            "); using snapshot of ", EXPLORENZ_SNAPSHOT_DATE)
    EXPLORENZ_SNAPSHOT
  })
  out
}
explorenz <- get_explorenz()
cat(sprintf("ExplorEnz: %d active entries\n\n", sum(explorenz)))

## ---- assemble --------------------------------------------------------------
get_n <- function(tbl, k) { v <- tbl[as.character(k)]; if (is.na(v)) 0L else as.integer(v) }

comp <- tibble(
  ec_class  = paste0("EC ", 1:7),
  enzyme    = unname(EC_NAMES[as.character(1:7)]),
  explorenz = as.integer(explorenz[as.character(1:7)]),
  uniprot   = map_int(1:7, ~ get_n(up_counts,    .x)),
  model     = map_int(1:7, ~ get_n(model_counts, .x))) %>%
  mutate(explorenz_pct = round(100 * explorenz / sum(explorenz), 1),
         uniprot_pct   = round(100 * uniprot   / sum(uniprot),   1),
         model_pct     = round(100 * model     / sum(model),     1))

print(as.data.frame(comp))
write_out(comp, "table_ec_class_distribution")

cat(sprintf("\nEC 2 + EC 3 in the model: %.1f%%\n",
            comp$model_pct[2] + comp$model_pct[3]))
cat(sprintf("EC 7: ExplorEnz %.1f%%, UniProtKB %.1f%%, model %.1f%%\n",
            comp$explorenz_pct[7], comp$uniprot_pct[7], comp$model_pct[7]))

## the Category-2 reactions show the opposite profile
c2 <- m %>% filter(cat_num == 2)
c2_counts <- table(classes_of(c2$ec_list))
cat("\nEC classes among the Category-2 reactions:\n")
for (k in names(sort(c2_counts, decreasing = TRUE)))
  cat(sprintf("  EC %s %-16s %3d (%.1f%%)\n", k, EC_NAMES[k],
              c2_counts[[k]], 100 * c2_counts[[k]] / nrow(c2)))

## ---- Figure 6 --------------------------------------------------------------
LEVELS <- c("ExplorEnz (all organisms)", "UniProtKB (S. cerevisiae)",
            "Yeast-GEM 9.0.2")

plot_df <- comp %>%
  select(ec_class, enzyme, explorenz_pct, uniprot_pct, model_pct) %>%
  pivot_longer(ends_with("_pct"), names_to = "source", values_to = "pct") %>%
  mutate(source = recode(source, explorenz_pct = LEVELS[1],
                                 uniprot_pct   = LEVELS[2],
                                 model_pct     = LEVELS[3]),
         source = factor(source, levels = LEVELS),
         xlab   = factor(paste0(ec_class, "\n", enzyme),
                         levels = paste0(comp$ec_class, "\n", comp$enzyme)))

p <- ggplot(plot_df, aes(x = xlab, y = pct, fill = source)) +
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = 0, ymax = Inf,
           fill = COL$accent, alpha = 0.07) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.78),
            vjust = -0.45, size = 2.7, colour = "grey25") +
  scale_fill_manual(values = setNames(c("#8FA8C8", COL$secondary, "#E8734A"), LEVELS),
                    name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Proportion of annotated reactions or entries") +
  theme_paper() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        axis.text.x = element_text(size = 9, lineheight = 1.1),
        legend.position = "bottom")

save_fig(p, "Figure6_EC_class_distribution", 9.5, 5.2)

message("\n07_ec_class_distribution.R done.")
