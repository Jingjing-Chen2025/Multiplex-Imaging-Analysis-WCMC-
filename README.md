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

## Segmentation update, 2026-10-22 (02_segment_measure.R v11.3)

v11.2 built the per-tile tissue mask from Otsu on 15-px-blurred DAPI and thresholded nuclei with a global per-tile
Otsu. In tiles that also held bright, densely packed tumour nuclei both thresholds rose and the loose stroma between
nests, with its dimmer nuclei, was treated as background: on test tiles of PM12_B3 and PM677_A8, 15-36% of nuclei,
nearly all stromal, were never segmented. v11.3 adds a local adaptive candidate-nucleus mask above an absolute floor
(3x the tile dark-background median), folds its dilation into the tissue mask and ORs it into the nuclear mask.
Parameters: ADAPT_W = 25 px, ADAPT_OFFSET = 0.02 of the tile maximum, FLOOR_K = 3, MASK_DILATE_UM = 8. Everything
downstream (opening, hull fill, watershed, 8-400 um2 area filter, ring measurement, gates) is unchanged.
`run_all.sh` is the three-worker from-scratch driver used for the v11.3 run (`run_02_parallel.sh` is the earlier one).

## Control-pairing fix, 2026-10-22 (02_segment_measure.R v11.4)

The matched secondary-only blank is now chosen by species and fluorophore. From the COMET run logs,
cycle 3 is anti-rabbit 647 + anti-mouse 555 and cycle 5 is anti-mouse 647 + anti-rabbit 555, so the
blanks are: CEBPD, PDGFRA -> NC3 Cy5; pan-CK -> NC3 TRITC; FOLH1, AR -> NC5 Cy5; DLL3 -> NC5 TRITC.
v11.3 paired FOLH1 with NC3 and PDGFRA / pan-CK with NC5 (fluorophore only, wrong species). The
`Negative-control pairing` comment in `00_config.R` already described the correct pairing for cycles 6
and 7. Segmentation (02), gates (03) and analyses (04) must be re-run; results before this fix used the
wrong blank for FOLH1, PDGFRA and pan-CK.
