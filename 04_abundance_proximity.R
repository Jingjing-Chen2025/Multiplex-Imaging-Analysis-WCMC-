###############################################################################
# 04_abundance_proximity.R  v11.5
#
# PURPOSE
#   Quantify ASPC-like stromal-cell abundance and spatial relationships to
#   epithelial and neuroendocrine tumour populations using the QC-aware outputs
#   from 03_classify.R v11.3+.
#
# v11.5 CHANGE
#   This script no longer resurrects ASPC calls by recomputing ASPC as plain
#   PDGF_pos. It respects `ASPC_callable` / `PDGF_usable` from 03_classify.R.
#   Raw PDGF-derived calls on non-callable slides are tracked as
#   `ASPC_not_callable` and shown separately in spatial maps, not as ASPC.
#
# INPUT
#   results_ASPC_NEPC/cells_classified.rds
#   results_ASPC_NEPC/channels.rds
#
# EXPECTED POPULATION GATES FROM 03_classify.R v11.3+
#   ASPC_raw      = PDGF_pos, with CEBPg strata retained
#   ASPC_callable = PDGF/PDGFRα slide-level marker usability flag
#   ASPC          = ASPC_raw & ASPC_callable
#   EPI_CA        = panCK_pos & (FOLH1_pos | AR_pos)
#   NE_CA         = DLL3_pos
#
# OUTPUT
#   aspc_abundance.csv
#   aspc_abundance_tests.csv
#   aspc_blank_false_positive_controls.csv
#   aspc_blank_adjusted_enrichment.csv
#   aspc_blank_adjusted_enrichment_tests.csv
#   spatial_bin_counts.csv
#   spatial_bin_correlations.csv
#   spatial_bin_correlation_tests.csv
#   aspc_proximity_perslide.csv
#   aspc_proximity_group_tests.csv
#   aspc_neighbourhood.csv
#   hypothesis_summary.csv
#   slide_qc_summary_v11_5_all_slides.csv
#   plots/04_figures.pdf and individual PDFs
#
# STATISTICAL NOTE
#   Unit of replication = slide/patient. Cell-level p-values are deliberately
#   not computed because they would be pseudoreplication.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(RANN)
})

###############################################################################
## USER SETTINGS
###############################################################################

OUT     <- "results_ASPC_NEPC"
PLOTDIR <- file.path(OUT, "plots")

## v11.5 default: include all slides present in cells_classified.rds, but keep
## ASPC non-callability explicit instead of treating non-callable raw calls as
## true ASPC.
USE_ALL_SLIDES_FOR_ANALYSIS <- TRUE
MAP_INCLUDE_OTHER <- TRUE
MAP_SHOW_ASPC_NOT_CALLABLE <- TRUE

MIN_N  <- 20
N_PERM <- 199

BIN_UM <- 250
MIN_BINS_FOR_RHO <- 10
MIN_TOTAL_CELLS_FOR_RHO <- 20

RADIUS_UM <- 50

set.seed(1)
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

PDF_DEV <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

save_pdf <- function(p, file, width = 10, height = 6) {
  ok <- tryCatch({
    ggsave(
      file.path(PLOTDIR, file),
      p,
      device = PDF_DEV,
      width = width,
      height = height,
      units = "in"
    )
    TRUE
  }, error = function(e) {
    warning("figure '", file, "' failed: ", conditionMessage(e), call. = FALSE)
    FALSE
  })

  if (ok) message("  wrote ", file)
  invisible(ok)
}

.figs <- list()

keep_fig <- function(p, w = 10, h = 6) {
  .figs[[length(.figs) + 1]] <<- list(plot = p, w = w, h = h)
  invisible(p)
}

###############################################################################
## SLIDE MANIFEST
###############################################################################

ROOT <- "/Volumes/T9/Case Studies"

MANIFEST <- data.frame(
  slide_id = c(
    "PM0T32", "PM677_A8", "PM12_B3",
    "PM1024_B", "PM423", "PM159T5"
  ),
  group = c(
    "NEPC", "NEPC", "NEPC",
    "ADENO", "ADENO", "ADENO"
  ),
  tiff = file.path(ROOT, c(
    "20260826_155725_1_ybucNp_JJ_full panel_PM0 T32/20260826_155725_1_ybucNp_JJ_full panel_PM0 T32.tiff",
    "20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8/20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8.tiff",
    "20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3/20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3.tiff",
    "20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B/20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B.tiff",
    "20260625_131559_2_ZxKUqB_JJ_full panel_PM423/20260625_131559_2_ZxKUqB_JJ_full panel_PM423.tiff",
    "20260625_131600_3_VHgSnt_JJ_full panel_PM159T5/20260625_131600_3_VHgSnt_JJ_full panel_PM159T5.tiff"
  )),
  stringsAsFactors = FALSE
)

GROUPS <- c("ADENO", "NEPC")

###############################################################################
## LOAD INPUTS
###############################################################################

cells_file <- file.path(OUT, "cells_classified.rds")
channels_file <- file.path(OUT, "channels.rds")

if (!file.exists(cells_file)) {
  stop("Run 03_classify.R v11.3+ first: missing cells_classified.rds.", call. = FALSE)
}

if (!file.exists(channels_file)) {
  stop("Missing channels.rds. Re-run 01_channels.R.", call. = FALSE)
}

cells <- readRDS(cells_file)

if (!nrow(cells)) {
  stop("cells_classified.rds has zero rows.", call. = FALSE)
}

z <- readRDS(channels_file)
PX_UM <- if (!is.null(z$px_um)) z$px_um else 12500 / 44643

present <- intersect(MANIFEST$slide_id, unique(cells$slide_id))
absent  <- setdiff(MANIFEST$slide_id, unique(cells$slide_id))

if (length(absent)) {
  warning(
    "No cells for: ", paste(absent, collapse = ", "),
    " -- excluded because they are absent from cells_classified.rds.",
    call. = FALSE
  )
}

if ("has_negctrl" %in% names(cells)) {
  .no_nc <- unique(cells$slide_id[!cells$has_negctrl])
  if (length(.no_nc)) {
    warning(
      "Slide(s) without negative controls in cells_classified.rds: ",
      paste(.no_nc, collapse = ", "),
      " -- blank-adjusted metrics are not interpretable for these. ",
      "Re-run 01-03 against the current acquisitions.",
      call. = FALSE
    )
  }
}

cells <- cells |>
  filter(slide_id %in% present) |>
  select(-any_of("group")) |>
  left_join(MANIFEST[, c("slide_id", "group")], by = "slide_id") |>
  mutate(group = factor(group, levels = GROUPS))

if (!"tier2" %in% names(cells)) cells$tier2 <- TRUE
if (!"tier3" %in% names(cells)) cells$tier3 <- TRUE

###############################################################################
## POPULATION DEFINITIONS -- QC-aware, do not resurrect raw ASPC calls
###############################################################################

required_gates <- paste0(
  c("CEBPg", "PDGF", "AR", "FOLH1", "panCK", "DLL3"),
  "_pos"
)

missing_gates <- setdiff(required_gates, names(cells))
if (length(missing_gates)) {
  stop(
    "Missing gates required for v11.5 population logic: ",
    paste(missing_gates, collapse = ", "),
    call. = FALSE
  )
}

if (!"CEBPg_neg" %in% names(cells)) {
  cells$CEBPg_neg <- !cells$CEBPg_pos
}

if (!"ASPC_raw" %in% names(cells)) {
  warning(
    "ASPC_raw is missing; constructing legacy raw ASPC as PDGF_pos. ",
    "Re-run 03_classify.R v11.3+ to retain full ASPC audit columns.",
    call. = FALSE
  )
  cells$ASPC_raw <- (cells$CEBPg_pos & cells$PDGF_pos) |
                    (cells$CEBPg_neg & cells$PDGF_pos)
}

if (!"ASPC_callable" %in% names(cells)) {
  warning(
    "ASPC_callable is missing; assuming TRUE for compatibility. ",
    "This can resurrect old PDGF-derived ASPC artefacts. Re-run 03_classify.R v11.3+.",
    call. = FALSE
  )
  cells$ASPC_callable <- TRUE
}

if (!"PDGF_usable" %in% names(cells)) {
  cells$PDGF_usable <- cells$ASPC_callable
}

cells <- cells |>
  mutate(
    ASPC_callable = ASPC_callable %in% TRUE,
    PDGF_usable = PDGF_usable %in% TRUE,

    ## Final ASPC used for abundance/spatial tests.
    ASPC = ASPC_raw & ASPC_callable,

    ## Raw non-callable calls are QC events, not biological ASPC.
    ASPC_not_callable = ASPC_raw & !ASPC_callable,

    ASPC_CEBPg_pos = ASPC & CEBPg_pos,
    ASPC_CEBPg_neg = ASPC & !CEBPg_pos,
    ASPC_raw_CEBPg_pos = ASPC_raw & CEBPg_pos,
    ASPC_raw_CEBPg_neg = ASPC_raw & !CEBPg_pos,

    EPI_CA = panCK_pos & (FOLH1_pos | AR_pos),
    NE_CA  = DLL3_pos
  )

bad_callable <- cells |>
  group_by(slide_id, group) |>
  summarise(
    ASPC_callable = all(ASPC_callable %in% TRUE),
    PDGF_usable = all(PDGF_usable %in% TRUE),
    ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
    ASPC = sum(ASPC, na.rm = TRUE),
    ASPC_not_callable = sum(ASPC_not_callable, na.rm = TRUE),
    .groups = "drop"
  )

if (any(!bad_callable$ASPC_callable)) {
  warning(
    "ASPC is not callable on slide(s): ",
    paste(bad_callable$slide_id[!bad_callable$ASPC_callable], collapse = ", "),
    ". Raw PDGF-derived calls are counted as ASPC_not_callable and are not used as ASPC.",
    call. = FALSE
  )
}

write.csv(
  bad_callable,
  file.path(OUT, "aspc_callability_used_by_04.csv"),
  row.names = FALSE
)

## Blank false-positive population controls if blank gates exist.
if (all(c(
  "PDGF_blank_pos",
  "panCK_blank_pos",
  "FOLH1_blank_pos",
  "AR_blank_pos",
  "DLL3_blank_pos"
) %in% names(cells))) {
  cells <- cells |>
    mutate(
      ASPC_blank = PDGF_blank_pos,
      EPI_CA_blank = panCK_blank_pos & (FOLH1_blank_pos | AR_blank_pos),
      NE_CA_blank = DLL3_blank_pos
    )
} else {
  warning(
    "Blank false-positive gate columns not all present; blank controls set to NA. ",
    "Re-run 02_segment_measure.R and 03_classify.R.",
    call. = FALSE
  )
  cells$ASPC_blank   <- NA
  cells$EPI_CA_blank <- NA
  cells$NE_CA_blank  <- NA
}

## Plotting label. Non-callable raw ASPC is explicit and never coloured as ASPC.
cells <- cells |>
  mutate(
    n_pop_flags = rowSums(cbind(ASPC, EPI_CA, NE_CA), na.rm = TRUE),
    multi_pop = n_pop_flags > 1,
    pop = case_when(
      multi_pop ~ "Multi",
      ASPC ~ "ASPC",
      NE_CA ~ "NE_CA",
      EPI_CA ~ "EPI_CA",
      ASPC_not_callable & MAP_SHOW_ASPC_NOT_CALLABLE ~ "ASPC_not_callable",
      TRUE ~ "Other"
    )
  )

###############################################################################
## ANALYSIS INCLUSION
###############################################################################

all_present_slides <- sort(unique(cells$slide_id))

tier2_flagged <- unique(cells$slide_id[cells$tier2 %in% c(TRUE, "TRUE", 1)])
tier3_flagged <- unique(cells$slide_id[cells$tier3 %in% c(TRUE, "TRUE", 1)])

if (USE_ALL_SLIDES_FOR_ANALYSIS) {
  tier2_ok <- all_present_slides
  tier3_ok <- all_present_slides
} else {
  tier2_ok <- tier2_flagged
  tier3_ok <- tier3_flagged
}

nslide <- cells |>
  distinct(slide_id, group) |>
  count(group, .drop = FALSE)

slide_qc_summary <- cells |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    tier2 = first(tier2),
    tier3 = first(tier3),
    ASPC_callable = all(ASPC_callable %in% TRUE),
    PDGF_usable = all(PDGF_usable %in% TRUE),
    n_ASPC = sum(ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(ASPC_not_callable, na.rm = TRUE),
    n_EPI_CA = sum(EPI_CA, na.rm = TRUE),
    n_NE_CA = sum(NE_CA, na.rm = TRUE),
    n_Multi = sum(multi_pop, na.rm = TRUE),
    n_Other = sum(pop == "Other", na.rm = TRUE),
    n_non_other = sum(pop != "Other", na.rm = TRUE),
    .groups = "drop"
  )

cat("\nSlides analysed:\n")
print(as.data.frame(nslide))
cat("\nSlide inclusion mode:\n")
cat("  USE_ALL_SLIDES_FOR_ANALYSIS = ", USE_ALL_SLIDES_FOR_ANALYSIS, "\n", sep = "")
cat("  slides used for abundance/tier2 analyses: ", paste(tier2_ok, collapse = ", "), "\n")
cat("  slides used for spatial/tier3 analyses:   ", paste(tier3_ok, collapse = ", "), "\n")
cat("\nOriginal tier2-flagged slides: ", paste(tier2_flagged, collapse = ", "), "\n")
cat("Original tier3-flagged slides: ", paste(tier3_flagged, collapse = ", "), "\n")
cat("\nSlide QC summary:\n")
print(as.data.frame(slide_qc_summary))

write.csv(
  slide_qc_summary,
  file.path(OUT, "slide_qc_summary_v11_5_all_slides.csv"),
  row.names = FALSE
)

if (any(nslide$n < 2)) {
  warning("A group has <2 slides -- no test is possible; effect sizes only.", call. = FALSE)
}

###############################################################################
## SLIDE-LEVEL COMPARISON HELPER
###############################################################################

cmp <- function(d, valcol, label,
                direction = c("lower_in_NEPC", "higher_in_NEPC")) {
  direction <- match.arg(direction)

  blank <- data.frame(
    metric = label,
    n_ADENO = 0L,
    n_NEPC = 0L,
    med_ADENO = NA_real_,
    med_NEPC = NA_real_,
    diff = NA_real_,
    log2_ratio = NA_real_,
    hodges_lehmann = NA_real_,
    p = NA_real_,
    expected = direction,
    supports = NA,
    row.names = NULL
  )

  if (!valcol %in% names(d) || !nrow(d)) return(blank)

  gs <- split(d[[valcol]], factor(d$group, levels = GROUPS))
  A <- gs$ADENO[is.finite(gs$ADENO)]
  N <- gs$NEPC[is.finite(gs$NEPC)]

  nA <- length(A)
  nN <- length(N)
  mA <- if (nA) median(A) else NA_real_
  mN <- if (nN) median(N) else NA_real_
  hl <- if (nA && nN) median(outer(N, A, "-")) else NA_real_

  p <- if (nA >= 2 && nN >= 2) {
    signif(suppressWarnings(wilcox.test(A, N))$p.value, 3)
  } else {
    NA_real_
  }

  sup <- if (is.na(hl)) {
    NA
  } else if (direction == "higher_in_NEPC") {
    hl > 0
  } else {
    hl < 0
  }

  data.frame(
    metric = label,
    n_ADENO = nA,
    n_NEPC = nN,
    med_ADENO = round(mA, 4),
    med_NEPC = round(mN, 4),
    diff = round(mN - mA, 4),
    log2_ratio = round(log2((mN + 1e-6) / (mA + 1e-6)), 3),
    hodges_lehmann = round(hl, 4),
    p = p,
    expected = direction,
    supports = sup,
    row.names = NULL
  )
}

###############################################################################
## AREA / TILE HELPERS
###############################################################################

infer_tile_size <- function(d, default = 2048L) {
  tx <- sort(unique(d$tile_x))
  ty <- sort(unique(d$tile_y))
  step <- c(diff(tx)[diff(tx) > 0], diff(ty)[diff(ty) > 0])
  step <- step[is.finite(step)]
  if (!length(step)) return(default)
  as.integer(round(median(step)))
}

tile_area_mm2_by_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id) |>
  summarise(
    inferred_tile_size_px = infer_tile_size(cur_data_all()),
    n_cellbearing_tiles = n_distinct(paste(tile_x, tile_y)),
    cellbearing_tile_mm2 =
      n_cellbearing_tiles * inferred_tile_size_px^2 * PX_UM^2 / 1e6,
    nuclear_area_mm2 = sum(area, na.rm = TRUE) * PX_UM^2 / 1e6,
    .groups = "drop"
  )

###############################################################################
## 1. ASPC ABUNDANCE
###############################################################################

ab_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    ASPC_callable = all(ASPC_callable %in% TRUE),
    PDGF_usable = all(PDGF_usable %in% TRUE),

    n_ASPC = sum(ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(ASPC_not_callable, na.rm = TRUE),
    n_ASPC_CEBPg_pos = sum(ASPC_CEBPg_pos, na.rm = TRUE),
    n_ASPC_CEBPg_neg = sum(ASPC_CEBPg_neg, na.rm = TRUE),

    n_EPI_CA = sum(EPI_CA, na.rm = TRUE),
    n_NE_CA = sum(NE_CA, na.rm = TRUE),
    n_Multi = sum(multi_pop, na.rm = TRUE),

    frac_ASPC = mean(ASPC, na.rm = TRUE),
    frac_ASPC_raw = mean(ASPC_raw, na.rm = TRUE),
    frac_ASPC_not_callable = mean(ASPC_not_callable, na.rm = TRUE),
    frac_ASPC_CEBPg_pos = mean(ASPC_CEBPg_pos, na.rm = TRUE),
    frac_ASPC_CEBPg_neg = mean(ASPC_CEBPg_neg, na.rm = TRUE),

    frac_EPI_CA = mean(EPI_CA, na.rm = TRUE),
    frac_NE_CA = mean(NE_CA, na.rm = TRUE),
    frac_Multi = mean(multi_pop, na.rm = TRUE),

    aspc_per_1000_cells = 1000 * mean(ASPC, na.rm = TRUE),
    aspc_raw_per_1000_cells = 1000 * mean(ASPC_raw, na.rm = TRUE),
    aspc_per_epi_ca = sum(ASPC, na.rm = TRUE) / pmax(sum(EPI_CA, na.rm = TRUE), 1),
    aspc_per_tumour = sum(ASPC, na.rm = TRUE) / pmax(sum(EPI_CA | NE_CA, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  left_join(tile_area_mm2_by_slide, by = "slide_id") |>
  mutate(
    aspc_per_mm2 = n_ASPC / pmax(cellbearing_tile_mm2, 1e-9),
    aspc_raw_per_mm2 = n_ASPC_raw / pmax(cellbearing_tile_mm2, 1e-9)
  )

cat("\nAbundance by slide:\n")
print(as.data.frame(ab_slide))

write.csv(ab_slide, file.path(OUT, "aspc_abundance.csv"), row.names = FALSE)

ab_tests <- rbind(
  cmp(ab_slide, "frac_ASPC", "ASPC fraction of all cells", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_1000_cells", "ASPC per 1000 cells", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_mm2", "ASPC per mm2 cell-bearing tile area", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_epi_ca", "ASPC per epithelial cancer cell", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_tumour", "ASPC per tumour cell", "higher_in_NEPC"),
  cmp(ab_slide, "frac_ASPC_CEBPg_pos", "ASPC CEBPg+ fraction of all cells", "higher_in_NEPC"),
  cmp(ab_slide, "frac_ASPC_CEBPg_neg", "ASPC CEBPg- fraction of all cells", "higher_in_NEPC"),
  cmp(ab_slide, "frac_EPI_CA", "Epithelial cancer fraction", "lower_in_NEPC"),
  cmp(ab_slide, "frac_NE_CA", "NE cancer fraction", "higher_in_NEPC")
)

ab_tests$p_adj <- signif(p.adjust(ab_tests$p, "BH"), 3)

cat("\nAbundance: ADENO vs NEPC, slide-level\n")
print(ab_tests, row.names = FALSE)

write.csv(ab_tests, file.path(OUT, "aspc_abundance_tests.csv"), row.names = FALSE)

###############################################################################
## 1B. BLANK FALSE-POSITIVE POPULATION CONTROLS
###############################################################################

blank_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    ASPC_callable = all(ASPC_callable %in% TRUE),

    n_ASPC = sum(ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(ASPC_not_callable, na.rm = TRUE),
    n_ASPC_blank = sum(ASPC_blank, na.rm = TRUE),

    n_EPI_CA = sum(EPI_CA, na.rm = TRUE),
    n_EPI_CA_blank = sum(EPI_CA_blank, na.rm = TRUE),
    n_NE_CA = sum(NE_CA, na.rm = TRUE),
    n_NE_CA_blank = sum(NE_CA_blank, na.rm = TRUE),

    frac_ASPC = mean(ASPC, na.rm = TRUE),
    frac_ASPC_raw = mean(ASPC_raw, na.rm = TRUE),
    frac_ASPC_not_callable = mean(ASPC_not_callable, na.rm = TRUE),
    frac_ASPC_blank = mean(ASPC_blank, na.rm = TRUE),

    frac_EPI_CA = mean(EPI_CA, na.rm = TRUE),
    frac_EPI_CA_blank = mean(EPI_CA_blank, na.rm = TRUE),
    frac_NE_CA = mean(NE_CA, na.rm = TRUE),
    frac_NE_CA_blank = mean(NE_CA_blank, na.rm = TRUE),

    ASPC_enrichment_over_blank = mean(ASPC, na.rm = TRUE) / pmax(mean(ASPC_blank, na.rm = TRUE), 1e-9),
    ASPC_raw_enrichment_over_blank = mean(ASPC_raw, na.rm = TRUE) / pmax(mean(ASPC_blank, na.rm = TRUE), 1e-9),
    EPI_CA_enrichment_over_blank = mean(EPI_CA, na.rm = TRUE) / pmax(mean(EPI_CA_blank, na.rm = TRUE), 1e-9),
    NE_CA_enrichment_over_blank = mean(NE_CA, na.rm = TRUE) / pmax(mean(NE_CA_blank, na.rm = TRUE), 1e-9),
    .groups = "drop"
  )

cat("\nBlank false-positive controls by slide:\n")
print(as.data.frame(blank_slide))

write.csv(blank_slide, file.path(OUT, "aspc_blank_false_positive_controls.csv"), row.names = FALSE)

###############################################################################
## 1C. BLANK-ADJUSTED MARKER ENRICHMENT
###############################################################################

aspc_enrich_slide <- blank_slide |>
  mutate(
    frac_ASPC_net = pmax(frac_ASPC - frac_ASPC_blank, 0),
    frac_ASPC_raw_net = pmax(frac_ASPC_raw - frac_ASPC_blank, 0),
    frac_EPI_CA_net = pmax(frac_EPI_CA - frac_EPI_CA_blank, 0),
    frac_NE_CA_net = pmax(frac_NE_CA - frac_NE_CA_blank, 0),
    ASPC_marker_to_blank = frac_ASPC / pmax(frac_ASPC_blank, 1e-9),
    ASPC_raw_marker_to_blank = frac_ASPC_raw / pmax(frac_ASPC_blank, 1e-9),
    EPI_CA_marker_to_blank = frac_EPI_CA / pmax(frac_EPI_CA_blank, 1e-9),
    NE_CA_marker_to_blank = frac_NE_CA / pmax(frac_NE_CA_blank, 1e-9)
  )

cat("\nBlank-adjusted marker enrichment by slide:\n")
print(as.data.frame(aspc_enrich_slide))

write.csv(aspc_enrich_slide, file.path(OUT, "aspc_blank_adjusted_enrichment.csv"), row.names = FALSE)

aspc_enrich_tests <- rbind(
  cmp(aspc_enrich_slide, "frac_ASPC_net", "ASPC net fraction above blank", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "ASPC_marker_to_blank", "ASPC marker-to-blank ratio", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "frac_EPI_CA_net", "EPI_CA net fraction above blank", "lower_in_NEPC"),
  cmp(aspc_enrich_slide, "EPI_CA_marker_to_blank", "EPI_CA marker-to-blank ratio", "lower_in_NEPC"),
  cmp(aspc_enrich_slide, "frac_NE_CA_net", "NE_CA net fraction above blank", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "NE_CA_marker_to_blank", "NE_CA marker-to-blank ratio", "higher_in_NEPC")
)

aspc_enrich_tests$p_adj <- signif(p.adjust(aspc_enrich_tests$p, "BH"), 3)

cat("\nBlank-adjusted enrichment tests, slide-level:\n")
print(aspc_enrich_tests, row.names = FALSE)

write.csv(aspc_enrich_tests, file.path(OUT, "aspc_blank_adjusted_enrichment_tests.csv"), row.names = FALSE)

p1c <- aspc_enrich_slide |>
  select(slide_id, group, frac_ASPC_net, ASPC_marker_to_blank,
         frac_EPI_CA_net, EPI_CA_marker_to_blank,
         frac_NE_CA_net, NE_CA_marker_to_blank) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Blank-adjusted marker enrichment",
    subtitle = "ASPC uses 03_classify.R callable gate; non-callable raw PDGF calls are excluded",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_pdf(p1c, "04_aspc_blank_adjusted_enrichment.pdf", 12, 6)
keep_fig(p1c, 12, 6)

###############################################################################
## 2. SPATIAL-BIN SPEARMAN RHO ANALYSIS
###############################################################################

make_bin_counts_one_slide <- function(d, sid, bin_um = BIN_UM) {
  if (!nrow(d)) return(NULL)

  bin_px <- bin_um / PX_UM

  d <- d |>
    mutate(
      x_bin = floor((x_px - min(x_px, na.rm = TRUE)) / bin_px),
      y_bin = floor((y_px - min(y_px, na.rm = TRUE)) / bin_px),
      STROMA_OTHER = pop == "Other"
    )

  d |>
    group_by(slide_id, group, x_bin, y_bin) |>
    summarise(
      n_total = n(),
      ASPC = sum(ASPC, na.rm = TRUE),
      ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
      ASPC_not_callable = sum(ASPC_not_callable, na.rm = TRUE),
      ASPC_CEBPg_pos = sum(ASPC_CEBPg_pos, na.rm = TRUE),
      ASPC_CEBPg_neg = sum(ASPC_CEBPg_neg, na.rm = TRUE),
      EPI_CA = sum(EPI_CA, na.rm = TRUE),
      NE_CA = sum(NE_CA, na.rm = TRUE),
      Multi = sum(multi_pop, na.rm = TRUE),
      STROMA_OTHER = sum(STROMA_OTHER, na.rm = TRUE),
      .groups = "drop"
    )
}

spatial_bin_counts <- do.call(rbind, lapply(tier3_ok, function(sid) {
  d <- cells |> filter(slide_id == sid)
  make_bin_counts_one_slide(d, sid, bin_um = BIN_UM)
}))

if (is.null(spatial_bin_counts) || !nrow(spatial_bin_counts)) {
  warning("No spatial bins available for rho analysis.", call. = FALSE)
  spatial_rho <- NULL
  spatial_rho_tests <- NULL
} else {
  write.csv(spatial_bin_counts, file.path(OUT, "spatial_bin_counts.csv"), row.names = FALSE)

  rho_one_pair <- function(b, xcol, ycol, pair_label) {
    x <- b[[xcol]]
    y <- b[[ycol]]
    ok <- is.finite(x) & is.finite(y)
    n_bins <- sum(ok)
    x_total <- sum(x[ok], na.rm = TRUE)
    y_total <- sum(y[ok], na.rm = TRUE)

    if (
      n_bins < MIN_BINS_FOR_RHO ||
      x_total < MIN_TOTAL_CELLS_FOR_RHO ||
      y_total < MIN_TOTAL_CELLS_FOR_RHO ||
      length(unique(x[ok])) < 2 ||
      length(unique(y[ok])) < 2
    ) {
      return(data.frame(
        slide_id = unique(b$slide_id)[1],
        group = unique(b$group)[1],
        bin_um = BIN_UM,
        pair = pair_label,
        x_population = xcol,
        y_population = ycol,
        n_bins = n_bins,
        x_total = x_total,
        y_total = y_total,
        rho = NA_real_,
        p = NA_real_,
        skipped = TRUE,
        row.names = NULL
      ))
    }

    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))

    data.frame(
      slide_id = unique(b$slide_id)[1],
      group = unique(b$group)[1],
      bin_um = BIN_UM,
      pair = pair_label,
      x_population = xcol,
      y_population = ycol,
      n_bins = n_bins,
      x_total = x_total,
      y_total = y_total,
      rho = round(unname(ct$estimate), 3),
      p = signif(ct$p.value, 3),
      skipped = FALSE,
      row.names = NULL
    )
  }

  spatial_rho <- do.call(rbind, lapply(split(spatial_bin_counts, spatial_bin_counts$slide_id), function(b) {
    rbind(
      rho_one_pair(b, "ASPC", "EPI_CA", "ASPC vs EPI_CA"),
      rho_one_pair(b, "ASPC", "NE_CA", "ASPC vs NE_CA"),
      rho_one_pair(b, "ASPC", "Multi", "ASPC vs Multi"),
      rho_one_pair(b, "STROMA_OTHER", "EPI_CA", "STROMA_OTHER vs EPI_CA"),
      rho_one_pair(b, "STROMA_OTHER", "NE_CA", "STROMA_OTHER vs NE_CA")
    )
  }))

  spatial_rho$group <- factor(spatial_rho$group, levels = GROUPS)
  cat("\nSpatial-bin Spearman rho by slide:\n")
  print(spatial_rho, row.names = FALSE)
  write.csv(spatial_rho, file.path(OUT, "spatial_bin_correlations.csv"), row.names = FALSE)

  spatial_rho_tests <- do.call(rbind, lapply(unique(spatial_rho$pair), function(pp) {
    d <- spatial_rho |> filter(pair == pp, !skipped)
    if (!nrow(d)) return(NULL)
    expected <- if (pp %in% c("ASPC vs EPI_CA", "ASPC vs NE_CA", "ASPC vs Multi")) {
      "higher_in_NEPC"
    } else {
      "lower_in_NEPC"
    }
    cmp(d, "rho", paste0("spatial rho: ", pp), expected)
  }))

  if (!is.null(spatial_rho_tests)) {
    spatial_rho_tests$p_adj <- signif(p.adjust(spatial_rho_tests$p, "BH"), 3)
    cat("\nSpatial-bin rho tests, slide-level:\n")
    print(spatial_rho_tests, row.names = FALSE)
    write.csv(spatial_rho_tests, file.path(OUT, "spatial_bin_correlation_tests.csv"), row.names = FALSE)
  }

  p_rho <- spatial_rho |>
    filter(!skipped) |>
    ggplot(aes(group, rho, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.5) +
    facet_wrap(~ pair, scales = "free_y") +
    labs(
      title = paste0("Spatial-bin Spearman correlation, ", BIN_UM, " µm bins"),
      subtitle = "ASPC excludes raw PDGF calls on non-callable slides",
      x = NULL,
      y = "Spearman rho"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  save_pdf(p_rho, "04_spatial_bin_rho.pdf", 12, 6)
  keep_fig(p_rho, 12, 6)

  scatter_dat <- bind_rows(
    spatial_bin_counts |> transmute(slide_id, group, pair = "ASPC vs EPI_CA", x = EPI_CA, y = ASPC),
    spatial_bin_counts |> transmute(slide_id, group, pair = "ASPC vs NE_CA", x = NE_CA, y = ASPC),
    spatial_bin_counts |> transmute(slide_id, group, pair = "STROMA_OTHER vs EPI_CA", x = EPI_CA, y = STROMA_OTHER),
    spatial_bin_counts |> transmute(slide_id, group, pair = "STROMA_OTHER vs NE_CA", x = NE_CA, y = STROMA_OTHER)
  )

  p_rho_scatter <- ggplot(scatter_dat, aes(log1p(x), log1p(y), colour = group)) +
    geom_point(alpha = 0.45, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_grid(pair ~ slide_id, scales = "free") +
    labs(
      title = paste0("Spatial-bin population co-occurrence, ", BIN_UM, " µm bins"),
      subtitle = "each point = one occupied spatial bin; axes are log1p cell counts",
      x = "log1p target population count per bin",
      y = "log1p ASPC or STROMA_OTHER count per bin"
    ) +
    theme_bw()

  save_pdf(p_rho_scatter, "04_spatial_bin_scatter.pdf", 14, 9)
  keep_fig(p_rho_scatter, 14, 9)
}

###############################################################################
## 3. SPATIAL PROXIMITY
###############################################################################

nn_med <- function(A, B) {
  if (!nrow(A) || !nrow(B)) return(NA_real_)
  median(RANN::nn2(B, A, k = 1)$nn.dists[, 1]) * PX_UM
}

prox_one <- function(d, from_lgl, to_lgl, from_lab, to_lab, sid, grp) {
  xy <- as.matrix(d[, c("x_px", "y_px")])
  from_lgl[is.na(from_lgl)] <- FALSE
  to_lgl[is.na(to_lgl)] <- FALSE

  nf <- sum(from_lgl)
  nt <- sum(to_lgl)

  if (nf < MIN_N || nt < MIN_N) {
    return(data.frame(
      slide_id = sid,
      group = grp,
      from = from_lab,
      to = to_lab,
      n_from = nf,
      n_to = nt,
      obs_um = NA_real_,
      null_um = NA_real_,
      z = NA_real_,
      frac_within_50um = NA_real_,
      frac_within_100um = NA_real_,
      p_perm = NA_real_,
      skipped = TRUE,
      row.names = NULL
    ))
  }

  dd <- RANN::nn2(
    xy[to_lgl, , drop = FALSE],
    xy[from_lgl, , drop = FALSE],
    k = 1
  )$nn.dists[, 1] * PX_UM

  obs <- median(dd)

  nullv <- replicate(N_PERM, {
    s <- sample.int(nrow(d))
    nn_med(
      xy[from_lgl[s], , drop = FALSE],
      xy[to_lgl[s], , drop = FALSE]
    )
  })

  mu  <- mean(nullv, na.rm = TRUE)
  sdv <- sd(nullv, na.rm = TRUE)

  data.frame(
    slide_id = sid,
    group = grp,
    from = from_lab,
    to = to_lab,
    n_from = nf,
    n_to = nt,
    obs_um = round(obs, 1),
    null_um = round(mu, 1),
    z = round((obs - mu) / max(sdv, 1e-9), 2),
    frac_within_50um = round(mean(dd <= 50), 3),
    frac_within_100um = round(mean(dd <= 100), 3),
    p_perm = signif((1 + sum(nullv <= obs, na.rm = TRUE)) / (N_PERM + 1), 3),
    skipped = FALSE,
    row.names = NULL
  )
}

prox <- do.call(rbind, lapply(tier3_ok, function(sid) {
  d <- cells |> filter(slide_id == sid)
  grp <- as.character(unique(d$group))[1]
  rbind(
    prox_one(d, d$ASPC, d$EPI_CA, "ASPC", "EPI_CA", sid, grp),
    prox_one(d, d$ASPC, d$NE_CA, "ASPC", "NE_CA", sid, grp)
  )
}))

prox$group <- factor(prox$group, levels = GROUPS)
cat("\nProximity by slide:\n")
print(prox, row.names = FALSE)
write.csv(prox, file.path(OUT, "aspc_proximity_perslide.csv"), row.names = FALSE)

skip <- prox |> filter(skipped)
if (nrow(skip)) {
  warning(
    "Proximity skipped (<", MIN_N, " cells) for ", nrow(skip),
    " slide-pair combinations -- see skipped column.",
    call. = FALSE
  )
}

prox_tests <- do.call(rbind, lapply(c("EPI_CA", "NE_CA"), function(tgt) {
  d <- prox |> filter(to == tgt, !skipped)
  if (!nrow(d)) return(NULL)
  rbind(
    cmp(d, "obs_um", paste0("ASPC->", tgt, ": median NN dist (um)"), "lower_in_NEPC"),
    cmp(d, "z", paste0("ASPC->", tgt, ": proximity z"), "lower_in_NEPC"),
    cmp(d, "frac_within_50um", paste0("ASPC->", tgt, ": frac within 50um"), "higher_in_NEPC"),
    cmp(d, "frac_within_100um", paste0("ASPC->", tgt, ": frac within 100um"), "higher_in_NEPC")
  )
}))

if (!is.null(prox_tests)) {
  prox_tests$p_adj <- signif(p.adjust(prox_tests$p, "BH"), 3)
  cat("\nProximity: ADENO vs NEPC, slide-level\n")
  print(prox_tests, row.names = FALSE)
  write.csv(prox_tests, file.path(OUT, "aspc_proximity_group_tests.csv"), row.names = FALSE)
}

###############################################################################
## 4. NEIGHBOURHOOD COMPOSITION AROUND ASPC
###############################################################################

rad_px <- RADIUS_UM / PX_UM

nbr <- do.call(rbind, lapply(tier3_ok, function(sid) {
  d <- cells |> filter(slide_id == sid)
  if (sum(d$ASPC, na.rm = TRUE) < MIN_N) return(NULL)

  xy <- as.matrix(d[, c("x_px", "y_px")])
  k <- min(100L, nrow(d))
  nn <- RANN::nn2(xy, xy[which(d$ASPC), , drop = FALSE], k = k)
  pops <- c("EPI_CA", "NE_CA", "ASPC", "Multi", "ASPC_not_callable", "Other")

  cnt <- vapply(pops, function(p) {
    sum(vapply(seq_len(nrow(nn$nn.idx)), function(i) {
      within <- nn$nn.dists[i, ] <= rad_px
      j <- nn$nn.idx[i, ][within]
      dist_j <- nn$nn.dists[i, ][within]
      j <- j[dist_j > 0]
      sum(d$pop[j] == p, na.rm = TRUE)
    }, numeric(1)))
  }, numeric(1))

  obs <- as.numeric(cnt / max(sum(cnt), 1))
  exp <- as.numeric(table(factor(d$pop, levels = pops))) / nrow(d)

  data.frame(
    slide_id = sid,
    group = as.character(unique(d$group))[1],
    neighbour = pops,
    obs_frac = round(obs, 3),
    exp_frac = round(exp, 3),
    log2_enrich = round(log2((obs + 1e-4) / (exp + 1e-4)), 2),
    row.names = NULL
  )
}))

if (!is.null(nbr) && nrow(nbr)) {
  nbr$group <- factor(nbr$group, levels = GROUPS)
  cat("\nNeighbourhood within ", RADIUS_UM, " um of ASPC:\n", sep = "")
  print(nbr, row.names = FALSE)
  write.csv(nbr, file.path(OUT, "aspc_neighbourhood.csv"), row.names = FALSE)
}

###############################################################################
## 5. HYPOTHESIS SUMMARY
###############################################################################

hyp <- rbind(
  ab_tests |> filter(metric == "ASPC per mm2 cell-bearing tile area") |> mutate(prediction = "(i) ASPC more abundant in NEPC"),
  if (exists("aspc_enrich_tests")) {
    aspc_enrich_tests |> filter(metric == "ASPC net fraction above blank") |> mutate(prediction = "(i-b) ASPC enriched above blank in NEPC")
  },
  if (!is.null(prox_tests)) {
    prox_tests |> filter(metric == "ASPC->EPI_CA: proximity z") |> mutate(prediction = "(ii) ASPC closer to epithelial cancer in NEPC")
  },
  if (!is.null(prox_tests)) {
    prox_tests |> filter(metric == "ASPC->NE_CA: proximity z") |> mutate(prediction = "(iii) ASPC closer to NE cancer in NEPC")
  },
  if (exists("spatial_rho_tests") && !is.null(spatial_rho_tests)) {
    spatial_rho_tests |> filter(metric == "spatial rho: ASPC vs EPI_CA") |> mutate(prediction = "(iv) ASPC spatial rho with epithelial cancer")
  },
  if (exists("spatial_rho_tests") && !is.null(spatial_rho_tests)) {
    spatial_rho_tests |> filter(metric == "spatial rho: ASPC vs NE_CA") |> mutate(prediction = "(v) ASPC spatial rho with NE cancer")
  },
  if (exists("spatial_rho_tests") && !is.null(spatial_rho_tests)) {
    spatial_rho_tests |> filter(metric == "spatial rho: STROMA_OTHER vs EPI_CA") |> mutate(prediction = "(vi) non-ASPC stromal spatial rho with epithelial cancer")
  }
) |>
  select(prediction, metric, n_ADENO, n_NEPC, med_ADENO, med_NEPC,
         hodges_lehmann, p, supports)

cat("\n================ HYPOTHESIS SUMMARY ================\n")
print(as.data.frame(hyp), row.names = FALSE)
cat("`supports` = direction of effect only. With n<=3 per arm, treat every row as descriptive.\n")

write.csv(hyp, file.path(OUT, "hypothesis_summary.csv"), row.names = FALSE)

###############################################################################
## 6. FIGURES
###############################################################################

message("\nWriting figures (PDF):")

p1 <- ab_slide |>
  select(slide_id, group, frac_ASPC, frac_ASPC_raw, frac_ASPC_not_callable,
         frac_ASPC_CEBPg_pos, frac_ASPC_CEBPg_neg,
         aspc_per_1000_cells, aspc_per_mm2, aspc_per_epi_ca, aspc_per_tumour) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Prediction (i): ASPC abundance",
    subtitle = "ASPC = raw PDGF-derived call only when slide-level PDGF/PDGFRα is callable",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_pdf(p1, "04_aspc_abundance.pdf", 12, 6)
keep_fig(p1, 12, 6)

p1b <- blank_slide |>
  select(slide_id, group, frac_ASPC, frac_ASPC_raw, frac_ASPC_not_callable,
         frac_ASPC_blank, frac_EPI_CA, frac_EPI_CA_blank, frac_NE_CA, frac_NE_CA_blank) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  mutate(
    control = ifelse(grepl("_blank$", metric), "blank_as_marker", "marker_or_qc"),
    population = gsub("^frac_|_blank$", "", metric)
  ) |>
  ggplot(aes(control, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.4) +
  facet_grid(population ~ group, scales = "free_y") +
  labs(
    title = "Blank false-positive controls and ASPC callability",
    subtitle = "negative-control channels processed as pseudo-markers; non-callable raw ASPC shown separately",
    x = NULL,
    y = "fraction of cells"
  ) +
  theme_bw()

save_pdf(p1b, "04_blank_false_positive_controls.pdf", 10, 7)
keep_fig(p1b, 10, 7)

if (any(!prox$skipped)) {
  pd <- prox |> filter(!skipped)

  p2 <- ggplot(pd, aes(group, obs_um, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.6) +
    facet_wrap(~ to, scales = "free_y") +
    labs(
      title = "ASPC -> target: median nearest-neighbour distance",
      subtitle = "ASPC excludes non-callable raw PDGF-derived calls",
      x = NULL,
      y = "median NN distance (µm)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  save_pdf(p2, "04_proximity_distance.pdf", 9, 4.5)
  keep_fig(p2, 9, 4.5)

  p3 <- ggplot(pd, aes(group, z, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.6) +
    facet_wrap(~ to) +
    labs(
      title = "Proximity z vs label-permuted null",
      subtitle = "negative = closer than chance; density-corrected within slide",
      x = NULL,
      y = "z"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  save_pdf(p3, "04_proximity_z.pdf", 9, 4.5)
  keep_fig(p3, 9, 4.5)

  p4 <- pd |>
    select(slide_id, group, to, frac_within_50um, frac_within_100um) |>
    pivot_longer(starts_with("frac_"), names_to = "radius", values_to = "frac") |>
    mutate(radius = sub("frac_within_", "", radius)) |>
    ggplot(aes(group, frac, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.4) +
    facet_grid(radius ~ to) +
    labs(title = "Fraction of ASPC with a target cell within radius", x = NULL, y = "fraction") +
    theme_bw() +
    theme(legend.position = "none")

  save_pdf(p4, "04_proximity_within_radius.pdf", 9, 6)
  keep_fig(p4, 9, 6)
}

if (!is.null(nbr) && nrow(nbr)) {
  p5 <- ggplot(nbr, aes(neighbour, slide_id, fill = log2_enrich)) +
    geom_tile() +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = paste0("Neighbourhood enrichment within ", RADIUS_UM, " µm of ASPC"), x = NULL, y = NULL) +
    theme_bw()

  save_pdf(p5, "04_neighbourhood.pdf", 8, 5)
  keep_fig(p5, 8, 5)
}

###############################################################################
## 7. SPATIAL MAPS
###############################################################################

map_dat <- cells |> filter(slide_id %in% tier3_ok)

if (!MAP_INCLUDE_OTHER) {
  map_dat <- map_dat |> filter(pop != "Other")
}

map_dat <- map_dat |>
  group_by(slide_id) |>
  mutate(
    x_mm = (x_px - min(x_px, na.rm = TRUE)) * PX_UM / 1000,
    y_mm = (y_px - min(y_px, na.rm = TRUE)) * PX_UM / 1000
  ) |>
  ungroup() |>
  mutate(
    panel = paste0(group, " | ", slide_id),
    pop = factor(pop, levels = c("ASPC", "EPI_CA", "NE_CA", "Multi", "ASPC_not_callable", "Other"))
  )

pop_cols <- c(
  ASPC = "#F8766D",
  EPI_CA = "#7CAE00",
  NE_CA = "#00BFC4",
  Multi = "#C77CFF",
  ASPC_not_callable = "grey35",
  Other = "grey85"
)

if (nrow(map_dat)) {
  pmap <- ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
    geom_point(data = dplyr::filter(map_dat, pop == "Other"), size = 0.08, alpha = 0.18) +
    geom_point(data = dplyr::filter(map_dat, pop == "ASPC_not_callable"), size = 0.10, alpha = 0.22) +
    geom_point(data = dplyr::filter(map_dat, !pop %in% c("Other", "ASPC_not_callable")), size = 0.25, alpha = 0.65) +
    facet_wrap(~ panel) +
    coord_equal() +
    scale_y_reverse() +
    scale_colour_manual(values = pop_cols, drop = FALSE) +
    labs(
      title = "Spatial distribution of ASPC and tumour populations",
      subtitle = "scaled blank correction + strict gate; raw PDGF calls on non-callable slides are grey, not ASPC",
      x = "mm",
      y = "mm"
    ) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  save_pdf(pmap, "04_spatial_maps.pdf", 12, 8)
  keep_fig(pmap, 12, 8)

  PDF_DEV(file.path(PLOTDIR, "04_spatial_maps_perslide.pdf"), width = 8, height = 8, onefile = TRUE)

  for (sid in unique(map_dat$slide_id)) {
    d <- map_dat |> filter(slide_id == sid)
    print(
      ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
        geom_point(data = dplyr::filter(d, pop == "Other"), size = 0.08, alpha = 0.18) +
        geom_point(data = dplyr::filter(d, pop == "ASPC_not_callable"), size = 0.10, alpha = 0.22) +
        geom_point(data = dplyr::filter(d, !pop %in% c("Other", "ASPC_not_callable")), size = 0.25, alpha = 0.65) +
        coord_equal() +
        scale_y_reverse() +
        scale_colour_manual(values = pop_cols, drop = FALSE) +
        labs(
          title = paste0(sid, "  [", unique(d$group), "]"),
          subtitle = paste0(
            "ASPC ", sum(d$pop == "ASPC", na.rm = TRUE),
            " | ASPC_not_callable ", sum(d$pop == "ASPC_not_callable", na.rm = TRUE),
            " | EPI_CA ", sum(d$pop == "EPI_CA", na.rm = TRUE),
            " | NE_CA ", sum(d$pop == "NE_CA", na.rm = TRUE),
            " | Multi ", sum(d$pop == "Multi", na.rm = TRUE),
            " | Other ", sum(d$pop == "Other", na.rm = TRUE)
          ),
          x = "mm",
          y = "mm"
        ) +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    )
  }

  invisible(dev.off())
  message("  wrote 04_spatial_maps_perslide.pdf")
} else {
  warning("No cells available for spatial maps after slide filtering.", call. = FALSE)
}

###############################################################################
## 8. MULTIPAGE SUMMARY PDF
###############################################################################

if (length(.figs)) {
  PDF_DEV(
    file.path(PLOTDIR, "04_figures.pdf"),
    width  = max(vapply(.figs, function(f) f$w, numeric(1))),
    height = max(vapply(.figs, function(f) f$h, numeric(1))),
    onefile = TRUE
  )

  for (f in .figs) print(f$plot)
  invisible(dev.off())
  message("  wrote 04_figures.pdf (", length(.figs), " pages)")
}

message("\n04_abundance_proximity.R complete -> ", PLOTDIR)
