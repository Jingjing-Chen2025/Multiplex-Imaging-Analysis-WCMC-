# Multiplex Imaging Analysis (WCMC) — ASPC / epithelial / NEPC

R pipeline for Lunaphore COMET sequential-IF prostate cancer sections
("JJ_full panel": CEBPg, FOLH1, AR, DLL3, pan-CK, PDGFa + secondary-only
negative controls). It segments nuclei, measures control-subtracted marker
signal, calls cell populations, analyses the ASPC neighbourhood and draws
per-biopsy figures.

## Cell populations

| Call | Rule |
|---|---|
| ASPC | PDGFa+ (split into CEBPg+ / CEBPg− subsets) |
| EPI_CA | pan-CK+ and (FOLH1+ or AR+) |
| NE_CA | DLL3+ |
| Multi | more than one of the above |

## Requirements

- R ≥ 4.3 with EBImage, dplyr, tidyr, purrr, tibble, ggplot2, MASS, RANN, xml2, reticulate
- Python virtualenv named `r-reticulate` with tifffile, zarr (≥ 3), numpy, imagecodecs
- Raw data: one OME-TIFF per slide under `ROOT` (default `/Volumes/T9/Case Studies`)

## Run

```bash
Rscript -e 'source("00_config.R"); source("01_channels.R")'   # slide table, channel map
./run_02_parallel.sh                                           # segmentation (hours; 3 workers)
Rscript -e 'source("03_classify.R")'                           # marker gates + population calls
Rscript -e 'source("04_abundance_proximity.R")'                # abundance, proximity, neighbourhood
Rscript -e 'source("05_visualize.R")'                          # per-slide figures
```

Results go to `OUT` (default `/Volumes/T9/results_ASPC_NEPC`). Set `ROOT` / `OUT`
before sourcing to change locations. Step 02 writes one checkpoint per slide
(`cells_raw_<slide>.rds`) and skips finished slides unless `FORCE_RERUN <- TRUE`.

| Script | What it does |
|---|---|
| `00_config.R` | slide table (single source of truth), dependency checks |
| `01_channels.R` | channel map from OME-XML, pixel size, DAPI registration check |
| `02_segment_measure.R` | slide tissue mask, tiled nuclear segmentation with halo, nuclear / 4-6-8 µm ring features after matched-control subtraction |
| `03_classify.R` | per-slide robust gates, blank false-positive audit, population calls |
| `04_abundance_proximity.R` | densities per mm² tissue, proximity vs permutation null, 50 µm ASPC neighbourhood |
| `05_visualize.R` | marker overview, cell-type map, ASPC neighbourhood map, hotspot crops, composition charts, per-cell CSV for QuPath/napari |

## Outputs (per slide, `OUT/plots/per_slide/<slide>/`)

`_A_marker_overview.png`, `_B_celltype_map.png`, `_C_aspc_neighbourhood.png`,
`_D_hotspot_<k>.png`, `_E_nbhd_composition.png`, `cells_<slide>.csv`.
All figures are also collected in `OUT/plots/05_all_slides.pdf`.

## Known caveats

- The panel stains PDGFa (ligand), not PDGFRA.
- PDGFa is called positive at a rate close to its matched negative-control
  channel on most slides (`qc_marker_usability_primary_6um.csv`); treat ASPC
  counts as provisional.
- PM0 T32 was stained without blocking and with different cycle 6/7
  secondaries than the other slides; its PDGFa/CEBPg positivity is not
  comparable.
- Some DLL3 calls in adenocarcinomas fall on cells autofluorescent in both
  TRITC channels.
- With three slides per group, group comparisons are descriptive only.
