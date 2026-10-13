# Multiplex Imaging Analysis (WCMC) — PDGFRA+ stroma / epithelial / NE cancer

R pipeline for Lunaphore COMET sequential-IF prostate cancer sections
("JJ_full panel": CEBPD, FOLH1, AR, DLL3, pan-CK, PDGFRA + secondary-only
negative controls). It segments nuclei, measures control-subtracted marker
signal, calls cell populations, and analyses abundance, proximity and
neighbourhoods. The COMET files originally named the PDGFRA channel "PDGFa"
and the CEBPD channel "CEBPg"; both were corrected on 2026-10-12.

## Cell populations

| Call | Rule | Role |
|---|---|---|
| PDGFRA | PDGFRA+ | main stromal population |
| PDGFRA_CEBPDpos | PDGFRA+ CEBPD+ | exploratory split of PDGFRA+ |
| PDGFRA_CEBPDneg | PDGFRA+ CEBPD− | exploratory split of PDGFRA+ |
| EPI_CA1 | AR+ FOLH1+ pan-CK+ | epithelial cancer |
| EPI_CA2 | AR+ FOLH1− pan-CK+ | epithelial cancer |
| NE_CA | DLL3+ | neuroendocrine cancer |
| Multi | meets more than one of PDGFRA / EPI_CA1 / EPI_CA2 / NE_CA | |

Blank false-positive controls apply the same rules to the matched
secondary-only channels.

## Requirements

- R ≥ 4.3 with EBImage, dplyr, tidyr, purrr, tibble, ggplot2, MASS, RANN, xml2, reticulate
- Python virtualenv named `r-reticulate` with tifffile, zarr, numpy, imagecodecs
- Raw data: one OME-TIFF per slide under `ROOT` (default `/Volumes/T9/Case Studies`)

## Run

```bash
Rscript -e 'source("00_config.R"); source("01_channels.R")'   # slide table, channel map
./run_02_parallel.sh                                           # segmentation (hours; 3 workers)
Rscript -e 'source("03_classify.R")'                           # marker gates + population calls
Rscript -e 'source("04_abundance_proximity.R")'                # abundance, proximity, neighbourhood
```

Results go to `OUT` (default `/Volumes/T9/Case Studies/multiplex_analysis_20261013_PDGFRA_v2`).
Set `ROOT` / `OUT` in R (or `OUT` in the environment for `run_02_parallel.sh`) to change locations.
Population rules only affect 03 and 04, so they can be re-run on an existing
`cells_raw.rds` in a few minutes.

## Outputs of 04

- Main analysis: `population_*`, `stroma_*`, `spatial_bin_*`, `hypothesis_summary.csv`;
  figures `plots/04_*.pdf`.
- Exploratory CEBPD split: the same analyses with the `exploratory_` prefix; figures
  `plots/04x_exploratory_*.pdf`.
- `plots/04_figures.pdf` collects every figure. Box plots use square panels of the same size.

## Notes

- `05_visualize.R` has not been updated. It expects the earlier ASPC outputs
  (`cells_spatial.rds`, `aspc_neighbourhood_composition.csv`) and does not run on
  the outputs of this version of 04.
- Slides: NEPC PM0T32 (2026-10-02 re-stain), PM677_A8, PM12_B3; ADENO PM1024_B, PM423, PM159T5.
- With three slides per group, group comparisons are descriptive only (minimum
  two-sided Wilcoxon p = 0.1).
