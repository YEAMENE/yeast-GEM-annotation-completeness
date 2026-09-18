## ============================================================================
## run_all.R — reproduce the complete analysis
##
## Working directory must be the repository root.
##   source("R/run_all.R")
## or, from a terminal:
##   Rscript R/run_all.R
##
## Scripts 01-04 and 08 need only data/yeast-GEM.xlsx.
## Scripts 05, 07 and 09 additionally need external data (see data/README.md);
## if a required file is missing they stop with an explanatory message and the
## remaining scripts still run.
## ============================================================================

scripts <- c(
  "R/01_model_composition.R",       # Section 3.1 · Figure 1 · Table S4
  "R/02_annotation_categories.R",   # Section 3.2 · Figure 2 · full dataset
  "R/03_category4_breakdown.R",     # Section 3.2 · Figure 3
  "R/04_confidence_crosstab.R",     # Section 3.2 · Figure 4 · Table S3
  "R/05_category3_recovery.R",      # Section 3.2 · Table S2      [UniProtKB]
  "R/06_genome_coverage.R",         # Section 3.3 · Figure 5      [SGD GFF3]
  "R/07_ec_class_distribution.R",   # Section 3.3 · Figure 6      [UniProtKB, ExplorEnz]
  "R/08_transport_ec7.R",           # Sections 2.5, 3.3
  "R/09_compartment_validation.R",  # Section 3.4 · Figure 7 · Table S5  [UniProtKB, SGD GAF]
  "R/10_external_resources.R"       # Section 2.2 · Table S1
)

t0 <- Sys.time()
status <- character(0)

for (s in scripts) {
  cat("\n\n", strrep("#", 74), "\n## ", s, "\n", strrep("#", 74), "\n", sep = "")
  ok <- tryCatch({ source(s, echo = FALSE); TRUE },
                 error = function(e) { message("FAILED: ", conditionMessage(e)); FALSE })
  status[s] <- if (ok) "ok" else "failed"
}

cat("\n\n", strrep("=", 74), "\n", sep = "")
cat("SUMMARY\n")
for (s in names(status)) cat(sprintf("  %-34s %s\n", basename(s), status[[s]]))
cat(sprintf("\nElapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Outputs written to: ", normalizePath(OUTPUT_DIR, mustWork = FALSE), "\n", sep = "")

if (any(status == "failed"))
  cat("\nSome scripts failed. Missing external data is the usual cause:\n",
      "see data/README.md for the download links.\n", sep = "")
