# External data

None of these files is redistributed in this repository: they remain subject to
their original licences. Download them into this directory before running the
pipeline.

The scripts locate the GFF3 and GAF files by extension, so the exact file name
does not matter for those two.

---

## 1. Yeast-GEM v9.0.2 — required by every script

Save as `data/yeast-GEM.xlsx`.

- Repository: https://github.com/SysBioChalmers/yeast-GEM
- Release v9.0.2: https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v9.0.2
- Archived: https://doi.org/10.5281/zenodo.14210050

---

## 2. UniProtKB reviewed *S. cerevisiae* proteome — scripts 05, 07, 09

Save as `data/uniprot_yeast_reviewed.tsv`.

`R/05_category3_recovery.R` downloads and caches this file automatically when
the network is available. To fetch it manually, use the REST API:

```
https://rest.uniprot.org/uniprotkb/stream?query=organism_id:559292%20AND%20reviewed:true&format=tsv&fields=accession,gene_primary,go_c,go_f,ec,xref_pfam,cc_subcellular_location,reviewed
```

Version used in the manuscript: release 2026_02 (10 June 2026), 6,733 entries.

---

## 3. SGD GO annotation file — script 09

Save as `data/sgd.gaf.gz` (or `.gaf`, uncompressed; either works).

- http://current.geneontology.org/annotations/sgd.gaf.gz

Contains the Gene Ontology annotations curated by SGD. Only Cellular Component
annotations with manually assigned evidence codes are used, i.e. everything
except IEA.

---

## 4. SGD S288C reference genome annotation — script 06

Save as `data/saccharomyces_cerevisiae_R64-5-1.gff` (any `.gff` name works).

- http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/

Pick the release folder matching the version cited in the manuscript
(R64.5.1, 29 May 2024), download the `.gff.gz` and decompress it.

Only features of type `gene` are used: SGD assigns separate types to non-coding
RNA genes (`ncRNA_gene`, `tRNA_gene`, `rRNA_gene`, `snoRNA_gene`), so this filter
already isolates protein-coding ORFs.

---

## 5. ExplorEnz — script 07

No download required. `R/07_ec_class_distribution.R` queries
https://www.enzyme-database.org/stats.php directly and falls back to a dated
snapshot if the site is unreachable.

If you update the snapshot, update its date at the same time, and keep both in
step with Supplementary Table S1 and with Section 2.2 of the manuscript.

---

## Checking the setup

```r
source("R/00_config.R")
file.exists(MODEL_XLSX)
file.exists(UNIPROT_TSV)
SGD_GAF   # NA if not found
SGD_GFF   # NA if not found
```
