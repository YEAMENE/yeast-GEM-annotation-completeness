# Yeast-GEM v9.0.2 — Reaction-level annotation completeness

Supporting code for:

> **Reaction-level evaluation of annotation completeness in the consensus yeast
> genome-scale model: a framework to guide curation priorities**
> De Martino L., Pesole G., Picardi E., Musicco C., Giannattasio S., Klapa M.I.

This repository contains the R code used to evaluate annotation completeness
(EC numbers, gene associations and compartment assignments) in the community
consensus yeast genome-scale metabolic model **Yeast-GEM v9.0.2**,
cross-validated against UniProtKB, SGD and ExplorEnz.

---

## Counting units

Three units are used, and the distinction matters for every reported percentage:

| Unit | n | Definition |
|---|---|---|
| **Reaction entries** | 4,131 | rows of the model. **Denominator of all percentages in the manuscript**, because EC and gene fields are attached to individual entries and each requires separate curation |
| Distinct annotated transformations | 3,701 | entries collapsed when stoichiometry, gene association and EC assignment all coincide |
| Distinct stoichiometries | 3,398 | entries collapsed on stoichiometry alone |

Distinct stoichiometries are **not** used to classify reactions: 75 of them fall
into different EC/gene categories depending on the compartment, so the category
is not well defined at that level. Reaction directionality is part of the
stoichiometric identity (`A => B` and `A <=> B` are distinct).

## Annotation categories

Each entry scores one point for at least one EC number and one for at least one
gene association, giving four categories. The two EC-bearing categories are
further split by EC specificity: **a** = at least one fully specified four-level
EC (e.g. 1.1.1.4); **b** = partial ECs only (e.g. 3.1.3.-, 2.-.-.-).

---

## How to run

```r
# 1. Dependencies (one-off)
install.packages(c("readxl", "dplyr", "stringr", "tidyr", "purrr", "readr",
                   "ggplot2", "forcats", "scales", "patchwork", "httr"))

# 2. Download the external data into data/ — see data/README.md

# 3. Working directory = repository root, then
source("R/run_all.R")
```

or from a terminal:

```bash
Rscript R/run_all.R
```

Paths are **relative** to the repository root. To use different locations set
the environment variables `YGEM_DATA` and `YGEM_OUTPUT`.

All outputs (figures as `.png` and `.pdf`, tables as `.csv`) are written to
`output/`, which is not version-controlled.

## Repository structure

```
.
├── R/
│   ├── 00_config.R                  # paths, model loading, shared definitions
│   ├── 01_model_composition.R       # Section 3.1 · Figure 1 · Table S2
│   ├── 02_annotation_categories.R   # Section 3.2 · Figure 2 · full dataset
│   ├── 03_category4_breakdown.R     # Section 3.2 · Figure 3
│   ├── 04_confidence_crosstab.R     # Section 3.2 · Figure 4 · Table S4
│   ├── 05_category3_recovery.R      # Section 3.2 · Table S3
│   ├── 06_genome_coverage.R         # Section 3.3 · Figure 5
│   ├── 07_ec_class_distribution.R   # Section 3.3 · Figure 6
│   ├── 08_transport_ec7.R           # Sections 2.5 and 3.3
│   ├── 09_compartment_validation.R  # Section 3.4 · Figure 7 · Table S5
│   ├── 10_external_resources.R      # Section 2.2 · Table S1
│   └── run_all.R
├── data/
│   └── README.md                    # download links for the external data
├── output/                          # created automatically, not version-controlled
├── .gitignore
├── LICENSE
└── README.md
```

Every script is self-contained: it sources `R/00_config.R` and loads the model
on its own, so any script can be re-run individually.

## Required data

`data/README.md` gives the exact download links for:

| File | Needed by |
|---|---|
| `yeast-GEM.xlsx` | all scripts |
| `uniprot_yeast_reviewed.tsv` | 05, 07, 09 (script 05 downloads and caches it) |
| `sgd.gaf` (or `.gaf.gz`) | 09 |
| `saccharomyces_cerevisiae_*.gff` | 06 |

Scripts that require a missing file stop with an explanatory message; the rest of
the pipeline still runs.

External datasets are not redistributed here because they remain subject to
their original licences.

## Outputs

**Figures** (no titles or sub-titles inside the plot: captions carry that
information)

| | |
|---|---|
| Figure 1 | compartmental distribution (A: reactions per compartment, exclusive vs shared; B: number of compartments per stoichiometry) |
| Figure 2 | EC/gene categories, on entries and on distinct annotated transformations |
| Figure 3 | composition of Category 4 |
| Figure 4 | categories versus the native SBML/COBRA confidence score |
| Figure 5 | genome coverage by chromosome (A: model genes; B: total complement with covered fraction) |
| Figure 6 | EC class distribution: model, ExplorEnz, UniProtKB |
| Figure 7 | compartmental validation (A: counts; B: inconsistency rate) |

**Supplementary tables**

| | |
|---|---|
| S1 | external resources with versions, release dates and file types |
| S2 | stoichiometries occurring in more than one compartment |
| S3 | Category-3 reactions recoverable to Category 1 from UniProtKB |
| S4 | categories versus native confidence score (a: cross-tab; b: the 61 Category-2 reactions; c: Category-1 reactions with no native score) |
| S5 | reactions whose compartment assignment is not consistently annotated |

`full_reaction_dataset.csv` contains all 4,131 entries with the category,
sub-category, EC specificity, compartment assignment, native score and
deduplication group identifiers.

## Reproducibility and external data versions

External resources change over time, so re-running the pipeline today may
produce slightly different numbers. This is expected, not a bug. The versions
used for the manuscript are listed in `output/TableS1_external_resources.csv`,
generated by `R/10_external_resources.R`:

- Yeast-GEM v9.0.2 (23 November 2024)
- UniProtKB release 2026_02 (10 June 2026)
- SGD R64.5.1 (29 May 2024)
- KEGG release 118.2 (1 June 2026)
- MetaCyc and BioCyc release 29.6 (26 February 2026)
- ExplorEnz, 6,972 active entries (accessed July 2026)

ExplorEnz has no numbered releases. `07_ec_class_distribution.R` queries it live
and falls back to a dated snapshot when the site is unreachable; the snapshot and
its date are declared at the top of that script and must be kept in step with
Supplementary Table S1.

## License

See `LICENSE`. Data from Yeast-GEM, UniProtKB, SGD, KEGG, MetaCyc, BioCyc and
ExplorEnz remain subject to their original licences and are not redistributed
here.

## Citation

De Martino L., Pesole G., Picardi E., Musicco C., Giannattasio S., Klapa M.I.
*Reaction-level evaluation of annotation completeness in the consensus yeast
genome-scale model: a framework to guide curation priorities.*

(Journal information and DOI will be added after publication.)
