# Host Phylogeny and Fungal Community Structure in Tropical Bark and Ambrosia Beetles

R analysis code associated with the MSc thesis:

**Host Phylogeny and Fungal Community Structure in Tropical Bark and Ambrosia Beetles**  
MSc Taxonomy, Biodiversity and Evolution, Imperial College London, 2026.

## Repository contents

The scripts are numbered in approximate workflow order:

1. `01_phylogeny_pruning.R` — prune the supplied COX1 host phylogenies to metabarcoded beetle tips.
2. `02_H1_alpha_diversity.R` — prepare H1 data and fit regional fungal alpha-diversity PGLMMs.
3. `03_H1_phylogeny_figures.R` — generate the phylogenetic alpha-diversity figures. Run after script 02 in the same R session.
4. `04_H2_beta_diversity.R` — fungal OTU filtering, binary Jaccard dissimilarities, sequential PERMANOVA, PERMDISP and PCoA.
5. `05_H3_phylogenetic_analysis.R` — pool replicate fungal samples by beetle-tree tip, calculate host patristic and fungal Jaccard distances, run Mantel tests, and generate the analysis objects used for the tanglegrams.
6. `06_H3_tanglegram_figure.R` — generate the final H3 tanglegram figure from the objects produced by script 05.
7. `07_focal_order_associations.R` — descriptive Microascales and Ophiostomatales host-fungus association analyses and figures.

## Data availability

The input datasets and phylogenies are **not included in this repository** because they were supplied by the research group and are not redistributed here.

Input file names and their analytical roles are documented in Table S1 of the thesis. To run the scripts locally, place the authorised input files in a data directory and set the environment variable `BEETLE_DATA_DIR` to that directory. If the environment variable is not set, the scripts assume a local folder named `data`.

Example in R:

```r
Sys.setenv(BEETLE_DATA_DIR = "/path/to/authorised/data")
```

Scripts 05 and 06 additionally use `BEETLE_RESULTS_DIR`; if this is not set, outputs are written to a local `results` directory.

## Note on the reconstructed H3 script

The original working script used for the H3 Mantel analysis was not retained. `05_H3_phylogenetic_analysis.R` was reconstructed from the documented analytical workflow and validated against the saved final H3 analysis object.

The reconstructed workflow reproduced the same retained host tips, fungal OTU sets, host patristic-distance matrices, fungal binary-Jaccard distance matrices, and Mantel statistics reported in the thesis:

- Scolytinae: 135 host tips, 2,463 fungal OTUs, Mantel r = 0.166, P = 0.0002.
- Platypodinae: 52 host tips, 1,630 fungal OTUs, Mantel r = 0.149, P = 0.0139.

## Software

Analyses were conducted in R. Package versions used for the submitted thesis are reported in the Methods section of the thesis.

## Scope

This repository contains analysis code only. Exploratory, debugging and superseded code has been omitted.
