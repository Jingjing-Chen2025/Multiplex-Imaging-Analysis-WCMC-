###############################################################################
# 04_abundance_proximity.R  v11.4-compatible
#
# PURPOSE
#   Quantify PDGFRA+ stromal-cell abundance and spatial relationships to AR+
#   epithelial (EPI_CA1 FOLH1+, EPI_CA2 FOLH1-) and neuroendocrine (NE_CA)
#   tumour-cell populations. The same analyses are repeated, as exploratory
#   outputs (exploratory_* files, 04x_* figures), for PDGFRA+ split by CEBPD.
#
# v11.4 CHANGE
#   PM0T32 re-acquired 20260826 on the standard 8-cycle protocol WITH
#   secondary-only negative controls (AR/DLL3 = Cycle 6). It is no longer a
#   protocol variant, so all six slides now carry the full blank-subtracted
#   feature set and the blank false-positive controls are defined for every
#   slide. MANIFEST path updated accordingly.
#
# v11.3 CHANGE (retained)
#   Relaxed filters so all slides present in cells_classified.rds are included
#   in abundance, spatial-bin rho, proximity, neighbourhood, and spatial-map
#   outputs. Original tier2/tier3 annotations are retained as QC labels but do
#   not exclude slides by default.
#
# INPUT
#   OUT/cells_classified.rds
#   OUT/channels.rds
#
# POPULATIONS (2026-10-13)
#   PDGFRA  = PDGFRA_pos                                 (main)
#   PDGFRA_CEBPDpos = PDGFRA_pos & CEBPD_pos             (exploratory)
#   PDGFRA_CEBPDneg = PDGFRA_pos & !CEBPD_pos            (exploratory)
#   EPI_CA1 = AR_pos & FOLH1_pos & panCK_pos
#   EPI_CA2 = AR_pos & !FOLH1_pos & panCK_pos
#   NE_CA   = DLL3_pos
#
# OUTPUT
#   (exploratory_* versions of the population/stroma/spatial/hypothesis files)
#   population_abundance.csv
#   population_abundance_tests.csv
#   population_blank_false_positive_controls.csv
#   population_blank_adjusted_enrichment.csv
#   population_blank_adjusted_enrichment_tests.csv
#   spatial_bin_counts.csv
#   spatial_bin_correlations.csv
#   spatial_bin_correlation_tests.csv
#   stroma_proximity_perslide.csv
#   stroma_proximity_group_tests.csv
#   stroma_neighbourhood.csv
#   hypothesis_summary.csv
#   slide_qc_summary_v11_3_all_slides.csv
#   plots/04_figures.pdf and individual PDFs
#
# STATISTICAL NOTE
#   Unit of replication = slide/patient.
#   With n = 3 vs 3 at best, the minimum two-sided Wilcoxon p-value is 0.1.
#   Cell-level p-values are deliberately not computed because they would be
#   pseudoreplication.
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

if (!exists("OUT")) OUT <- "/Volumes/T9/Case Studies/multiplex_analysis_20261013_PDGFRA_v2"
PLOTDIR <- file.path(OUT, "plots")

## v11.3 default: include all slides present in cells_classified.rds.
## Set FALSE to restore tier2/tier3 QC-filtered behaviour.
USE_ALL_SLIDES_FOR_ANALYSIS <- TRUE

## Include Other cells in grey in spatial maps so every available slide appears.
MAP_INCLUDE_OTHER <- TRUE

## Spatial analysis settings.
MIN_N  <- 20
N_PERM <- 199

## Spatial-bin rho settings.
BIN_UM <- 250
MIN_BINS_FOR_RHO <- 10
MIN_TOTAL_CELLS_FOR_RHO <- 20

## Stromal-cell-centred neighbourhood radius.
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
    warning(
      "figure '", file, "' failed: ",
      conditionMessage(e),
      call. = FALSE
    )
    FALSE
  })

  if (ok) message("  wrote ", file)
  invisible(ok)
}

.figs <- list()

## Box/point plots share one panel shape: square panels of PANEL_IN inches,
## with the PDF sized from the facet layout, so the y axis has the same length
## in every figure.
PANEL_IN <- 2.6
BOX_THEME <- theme(aspect.ratio = 1)

fig_w <- function(ncol) ncol * (PANEL_IN + 0.35) + 1.4
fig_h <- function(nrow, aspect = 1) nrow * (PANEL_IN * aspect + 0.45) + 1.3

keep_fig <- function(p, w = 10, h = 6) {
  .figs[[length(.figs) + 1]] <<- list(plot = p, w = w, h = h)
  invisible(p)
}

###############################################################################
## SLIDE MANIFEST
###############################################################################
## All six slides now share the 8-cycle protocol with secondary-only negative
## controls (AR/DLL3 = Cycle 6). The tiff paths are recorded here for
## provenance only; this script reads cells_classified.rds, not the images.

if (!exists("ROOT")) ROOT <- "/Volumes/T9/Case Studies"

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
    "20261002_161511_1_INNdFn_JJ_full panel_PMO T32/20261002_161511_1_INNdFn_JJ_full panel_PMO T32.tiff",
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
  stop("Run 03_classify.R v11.1+ first: missing cells_classified.rds.", call. = FALSE)
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
    "No cells for: ",
    paste(absent, collapse = ", "),
    " -- excluded because they are absent from cells_classified.rds.",
    call. = FALSE
  )
}

## Freshness guard: every slide should now carry negative controls. A stale
## PM0T32 from the retired 6-cycle acquisition would have has_negctrl = FALSE
## and would make the blank-adjusted sections below meaningless for that slide.
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

###############################################################################
## POPULATION DEFINITIONS
###############################################################################

required_gates <- paste0(
  c("CEBPD", "PDGFRA", "AR", "FOLH1", "panCK", "DLL3"),
  "_pos"
)

missing_gates <- setdiff(required_gates, names(cells))

if (length(missing_gates)) {
  stop(
    "Missing gates required for population logic: ",
    paste(missing_gates, collapse = ", "),
    call. = FALSE
  )
}

## Population definitions (2026-10-13, PDGFRA run).
## PDGFRA+ is the main stromal population; the CEBPD+/- split is exploratory.
TUMOUR_POPS <- c("EPI_CA1", "EPI_CA2", "NE_CA")
SPLIT_POPS  <- c("PDGFRA_CEBPDpos", "PDGFRA_CEBPDneg")
POPS        <- c("PDGFRA", TUMOUR_POPS)

## Expected direction of each population's abundance in NEPC vs ADENO.
POP_DIRECTION <- c(
  PDGFRA = "higher_in_NEPC",
  PDGFRA_CEBPDpos = "higher_in_NEPC",
  PDGFRA_CEBPDneg = "higher_in_NEPC",
  EPI_CA1 = "lower_in_NEPC",
  EPI_CA2 = "lower_in_NEPC",
  NE_CA = "higher_in_NEPC"
)

cells <- cells |>
  mutate(
    PDGFRA          = PDGFRA_pos,
    PDGFRA_CEBPDpos = PDGFRA_pos & CEBPD_pos,
    PDGFRA_CEBPDneg = PDGFRA_pos & !CEBPD_pos,
    EPI_CA1         = AR_pos & FOLH1_pos & panCK_pos,
    EPI_CA2         = AR_pos & !FOLH1_pos & panCK_pos,
    NE_CA           = DLL3_pos
  )

## Blank false-positive population controls: same gates on the matched
## secondary-only channels ("-" gates included).
if (all(paste0(c("PDGFRA", "CEBPD", "panCK", "FOLH1", "AR", "DLL3"),
               "_blank_pos") %in% names(cells))) {
  cells <- cells |>
    mutate(
      PDGFRA_blank          = PDGFRA_blank_pos,
      PDGFRA_CEBPDpos_blank = PDGFRA_blank_pos & CEBPD_blank_pos,
      PDGFRA_CEBPDneg_blank = PDGFRA_blank_pos & !CEBPD_blank_pos,
      EPI_CA1_blank         = AR_blank_pos & FOLH1_blank_pos & panCK_blank_pos,
      EPI_CA2_blank         = AR_blank_pos & !FOLH1_blank_pos & panCK_blank_pos,
      NE_CA_blank           = DLL3_blank_pos
    )
} else {
  warning(
    "Blank false-positive gate columns not all present; blank controls set to NA.",
    call. = FALSE
  )
  for (p in c(POPS, SPLIT_POPS)) cells[[paste0(p, "_blank")]] <- NA
}

## Plotting labels. Overlapping calls (stromal and tumour gates both met) are
## labelled Multi. pop_split is the same label with PDGFRA replaced by its
## CEBPD+/- subset, used only by the exploratory analyses.
cells$n_pop_flags <- rowSums(as.matrix(as.data.frame(cells)[, POPS]), na.rm = TRUE)
cells$multi_pop <- cells$n_pop_flags > 1
cells$pop <- "Other"
for (p in POPS) cells$pop[as.data.frame(cells)[[p]] %in% TRUE] <- p
cells$pop[cells$multi_pop] <- "Multi"
cells$pop_split <- cells$pop
cells$pop_split[cells$pop == "PDGFRA" & cells$CEBPD_pos %in% TRUE] <- "PDGFRA_CEBPDpos"
cells$pop_split[cells$pop == "PDGFRA" & !(cells$CEBPD_pos %in% TRUE)] <- "PDGFRA_CEBPDneg"

###############################################################################
## ANALYSIS INCLUSION
###############################################################################

all_present_slides <- sort(unique(cells$slide_id))

tier2_flagged <- if ("tier2" %in% names(cells)) {
  unique(cells$slide_id[cells$tier2 %in% c(TRUE, "TRUE", 1)])
} else {
  all_present_slides
}

tier3_flagged <- if ("tier3" %in% names(cells)) {
  unique(cells$slide_id[cells$tier3 %in% c(TRUE, "TRUE", 1)])
} else {
  all_present_slides
}

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
    tier2 = if ("tier2" %in% names(cells)) unique(tier2)[1] else TRUE,
    tier3 = if ("tier3" %in% names(cells)) unique(tier3)[1] else TRUE,
    n_PDGFRA = sum(PDGFRA, na.rm = TRUE),
    n_PDGFRA_CEBPDpos = sum(PDGFRA_CEBPDpos, na.rm = TRUE),
    n_PDGFRA_CEBPDneg = sum(PDGFRA_CEBPDneg, na.rm = TRUE),
    n_EPI_CA1 = sum(EPI_CA1, na.rm = TRUE),
    n_EPI_CA2 = sum(EPI_CA2, na.rm = TRUE),
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
  file.path(OUT, "slide_qc_summary_v11_3_all_slides.csv"),
  row.names = FALSE
)

if (any(nslide$n < 2)) {
  warning(
    "A group has <2 slides -- no test is possible; effect sizes only.",
    call. = FALSE
  )
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

  dx <- diff(tx)
  dy <- diff(ty)

  step <- c(dx[dx > 0], dy[dy > 0])
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
## ANALYSES 1-6 ARE RUN TWICE:
##   main:        PDGFRA+ stroma                     (unprefixed files, 04_* figures)
##   exploratory: PDGFRA+ split by CEBPD+ / CEBPD-   (exploratory_* files, 04x_* figures)
###############################################################################

run_set <- function(STROMA_POPS, LABCOL, PFX, FPFX, DEF_TEXT, CAPTION = NULL) {
POPS <- c(STROMA_POPS, TUMOUR_POPS)
if (!is.null(CAPTION)) {
  .save_pdf <- save_pdf
  .keep_fig <- keep_fig
  save_pdf <- function(p, ...) .save_pdf(p + labs(caption = CAPTION), ...)
  keep_fig <- function(p, ...) .keep_fig(p + labs(caption = CAPTION), ...)
}

###############################################################################
## 1. POPULATION ABUNDANCE
###############################################################################

cell_df <- as.data.frame(cells)

per_slide <- function(ids, f) {
  do.call(rbind, lapply(ids, function(sid) {
    d <- cell_df[cell_df$slide_id == sid, , drop = FALSE]
    cbind(
      data.frame(slide_id = sid, group = as.character(d$group[1]),
                 stringsAsFactors = FALSE),
      f(d)
    )
  }))
}

ab_slide <- per_slide(tier2_ok, function(d) {
  out <- data.frame(n_cells = nrow(d))
  tum <- d$EPI_CA1 | d$EPI_CA2 | d$NE_CA
  epi <- d$EPI_CA1 | d$EPI_CA2
  for (p in POPS) out[[paste0("n_", p)]] <- sum(d[[p]], na.rm = TRUE)
  out$n_Multi <- sum(d$multi_pop, na.rm = TRUE)
  for (p in POPS) out[[paste0("frac_", p)]] <- mean(d[[p]], na.rm = TRUE)
  out$frac_Multi <- mean(d$multi_pop, na.rm = TRUE)
  for (s in STROMA_POPS) {
    out[[paste0(s, "_per_1000_cells")]] <- 1000 * mean(d[[s]], na.rm = TRUE)
    out[[paste0(s, "_per_epi_ca")]] <- sum(d[[s]], na.rm = TRUE) / max(sum(epi, na.rm = TRUE), 1)
    out[[paste0(s, "_per_tumour")]] <- sum(d[[s]], na.rm = TRUE) / max(sum(tum, na.rm = TRUE), 1)
  }
  if (length(STROMA_POPS) == 2) {
    s1 <- sum(d[[STROMA_POPS[1]]], na.rm = TRUE)
    s2 <- sum(d[[STROMA_POPS[2]]], na.rm = TRUE)
    out$split_share <- s1 / max(s1 + s2, 1)
  }
  out$EPI_CA1_frac_of_EPI_CA <- sum(d$EPI_CA1, na.rm = TRUE) / max(sum(epi, na.rm = TRUE), 1)
  out
}) |>
  left_join(tile_area_mm2_by_slide, by = "slide_id")

## Cell-bearing tile area, not a true tissue-mask area.
for (p in POPS) {
  ab_slide[[paste0(p, "_per_mm2")]] <-
    ab_slide[[paste0("n_", p)]] / pmax(ab_slide$cellbearing_tile_mm2, 1e-9)
}

cat("\nAbundance by slide:\n")
print(ab_slide)

write.csv(ab_slide, file.path(OUT, paste0(PFX, "population_abundance.csv")), row.names = FALSE)

ab_tests <- do.call(rbind, c(
  lapply(POPS, function(p) rbind(
    cmp(ab_slide, paste0("frac_", p), paste0(p, " fraction of all cells"), POP_DIRECTION[[p]]),
    cmp(ab_slide, paste0(p, "_per_mm2"), paste0(p, " per mm2 cell-bearing tile area"), POP_DIRECTION[[p]])
  )),
  lapply(STROMA_POPS, function(s) rbind(
    cmp(ab_slide, paste0(s, "_per_epi_ca"), paste0(s, " per epithelial cancer cell (EPI_CA1+EPI_CA2)"), "higher_in_NEPC"),
    cmp(ab_slide, paste0(s, "_per_tumour"), paste0(s, " per tumour cell (EPI_CA1+EPI_CA2+NE_CA)"), "higher_in_NEPC")
  )),
  list(
    if (length(STROMA_POPS) == 2)
      cmp(ab_slide, "split_share", paste0(STROMA_POPS[1], " share of PDGFRA+ cells"), "higher_in_NEPC"),
    cmp(ab_slide, "EPI_CA1_frac_of_EPI_CA", "EPI_CA1 share of AR+ panCK+ cells", "lower_in_NEPC")
  )
))

ab_tests$p_adj <- signif(p.adjust(ab_tests$p, "BH"), 3)

cat("\nAbundance: ADENO vs NEPC, slide-level\n")
print(ab_tests, row.names = FALSE)

write.csv(ab_tests, file.path(OUT, paste0(PFX, "population_abundance_tests.csv")), row.names = FALSE)

###############################################################################
## 1B. BLANK FALSE-POSITIVE POPULATION CONTROLS
###############################################################################

blank_slide <- per_slide(tier2_ok, function(d) {
  out <- data.frame(n_cells = nrow(d))
  for (p in POPS) {
    out[[paste0("n_", p)]] <- sum(d[[p]], na.rm = TRUE)
    out[[paste0("n_", p, "_blank")]] <- sum(d[[paste0(p, "_blank")]], na.rm = TRUE)
  }
  for (p in POPS) {
    out[[paste0("frac_", p)]] <- mean(d[[p]], na.rm = TRUE)
    out[[paste0("frac_", p, "_blank")]] <- mean(d[[paste0(p, "_blank")]], na.rm = TRUE)
  }
  out
})

cat("\nBlank false-positive controls by slide:\n")
print(blank_slide)

write.csv(blank_slide, file.path(OUT, paste0(PFX, "population_blank_false_positive_controls.csv")), row.names = FALSE)

###############################################################################
## 1C. BLANK-ADJUSTED MARKER ENRICHMENT
###############################################################################

enrich_slide <- blank_slide
for (p in POPS) {
  f  <- enrich_slide[[paste0("frac_", p)]]
  fb <- enrich_slide[[paste0("frac_", p, "_blank")]]
  enrich_slide[[paste0("frac_", p, "_net")]] <- pmax(f - fb, 0)
  ## Undefined (NA) when the blank call count is zero, rather than f / 1e-9.
  enrich_slide[[paste0(p, "_marker_to_blank")]] <- ifelse(fb > 0, f / fb, NA_real_)
}

cat("\nBlank-adjusted marker enrichment by slide:\n")
print(enrich_slide)

write.csv(enrich_slide, file.path(OUT, paste0(PFX, "population_blank_adjusted_enrichment.csv")), row.names = FALSE)

enrich_tests <- do.call(rbind, lapply(POPS, function(p) rbind(
  cmp(enrich_slide, paste0("frac_", p, "_net"), paste0(p, " net fraction above blank"), POP_DIRECTION[[p]]),
  cmp(enrich_slide, paste0(p, "_marker_to_blank"), paste0(p, " marker-to-blank ratio"), POP_DIRECTION[[p]])
)))

enrich_tests$p_adj <- signif(p.adjust(enrich_tests$p, "BH"), 3)

cat("\nBlank-adjusted enrichment tests, slide-level:\n")
print(enrich_tests, row.names = FALSE)

write.csv(enrich_tests, file.path(OUT, paste0(PFX, "population_blank_adjusted_enrichment_tests.csv")), row.names = FALSE)

p1c <- enrich_slide |>
  select(slide_id, group, ends_with("_marker_to_blank")) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(sub("_marker_to_blank$", "", metric), levels = POPS),
    group = factor(group, levels = GROUPS)
  ) |>
  ggplot(aes(group, value, colour = group)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "Marker-to-blank ratio per population",
    subtitle = "fraction called with marker gates / fraction called with the same gates on blank channels; dashed = 1",
    x = NULL,
    y = "marker / blank"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p1c <- p1c + BOX_THEME
save_pdf(p1c, paste0(FPFX, "blank_adjusted_enrichment.pdf"), fig_w(length(POPS)), fig_h(1))
keep_fig(p1c, fig_w(length(POPS)), fig_h(1))

###############################################################################
## 2. SPATIAL-BIN SPEARMAN RHO ANALYSIS
###############################################################################

make_bin_counts_one_slide <- function(d, sid, bin_um = BIN_UM) {
  if (!nrow(d)) return(NULL)

  bin_px <- bin_um / PX_UM

  d$x_bin <- floor((d$x_px - min(d$x_px, na.rm = TRUE)) / bin_px)
  d$y_bin <- floor((d$y_px - min(d$y_px, na.rm = TRUE)) / bin_px)

  ## Operational stromal/other compartment: no population gate met.
  ## This is not a dedicated fibroblast-marker-positive compartment.
  d$STROMA_OTHER <- d$n_pop_flags == 0

  key <- paste(d$x_bin, d$y_bin)
  agg <- function(v) as.integer(tapply(v, key, sum, na.rm = TRUE))

  out <- data.frame(
    slide_id = sid,
    group = as.character(d$group[1]),
    x_bin = as.integer(tapply(d$x_bin, key, `[`, 1)),
    y_bin = as.integer(tapply(d$y_bin, key, `[`, 1)),
    n_total = as.integer(table(key)[sort(unique(key))]),
    stringsAsFactors = FALSE
  )
  for (p in c(POPS, "STROMA_OTHER")) out[[p]] <- agg(d[[p]])
  out$Multi <- agg(d$multi_pop)
  out
}

spatial_bin_counts <- do.call(rbind, lapply(tier3_ok, function(sid) {
  make_bin_counts_one_slide(cell_df[cell_df$slide_id == sid, , drop = FALSE], sid)
}))

RHO_PAIRS <- rbind(
  expand.grid(x = c(STROMA_POPS, "STROMA_OTHER"), y = TUMOUR_POPS,
              stringsAsFactors = FALSE),
  if (length(STROMA_POPS) == 2)
    data.frame(x = STROMA_POPS[1], y = STROMA_POPS[2], stringsAsFactors = FALSE)
)
RHO_PAIRS$label <- paste(RHO_PAIRS$x, "vs", RHO_PAIRS$y)

if (is.null(spatial_bin_counts) || !nrow(spatial_bin_counts)) {
  warning("No spatial bins available for rho analysis.", call. = FALSE)
  spatial_rho <- NULL
  spatial_rho_tests <- NULL

} else {
  write.csv(spatial_bin_counts, file.path(OUT, paste0(PFX, "spatial_bin_counts.csv")), row.names = FALSE)

  rho_one_pair <- function(b, xcol, ycol, pair_label) {
    x <- b[[xcol]]
    y <- b[[ycol]]

    ok <- is.finite(x) & is.finite(y)

    n_bins <- sum(ok)
    x_total <- sum(x[ok], na.rm = TRUE)
    y_total <- sum(y[ok], na.rm = TRUE)

    res <- data.frame(
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
    )

    if (
      n_bins < MIN_BINS_FOR_RHO ||
      x_total < MIN_TOTAL_CELLS_FOR_RHO ||
      y_total < MIN_TOTAL_CELLS_FOR_RHO ||
      length(unique(x[ok])) < 2 ||
      length(unique(y[ok])) < 2
    ) {
      return(res)
    }

    ct <- suppressWarnings(
      cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
    )

    res$rho <- round(unname(ct$estimate), 3)
    res$p <- signif(ct$p.value, 3)
    res$skipped <- FALSE
    res
  }

  spatial_rho <- do.call(rbind, lapply(
    split(spatial_bin_counts, spatial_bin_counts$slide_id),
    function(b) do.call(rbind, lapply(seq_len(nrow(RHO_PAIRS)), function(i)
      rho_one_pair(b, RHO_PAIRS$x[i], RHO_PAIRS$y[i], RHO_PAIRS$label[i])))
  ))

  spatial_rho$group <- factor(spatial_rho$group, levels = GROUPS)

  cat("\nSpatial-bin Spearman rho by slide:\n")
  print(spatial_rho, row.names = FALSE)

  write.csv(spatial_rho, file.path(OUT, paste0(PFX, "spatial_bin_correlations.csv")), row.names = FALSE)

  spatial_rho_tests <- do.call(rbind, lapply(RHO_PAIRS$label, function(pp) {
    d <- spatial_rho |>
      filter(pair == pp, !skipped)

    if (!nrow(d)) return(NULL)

    expected <- if (startsWith(pp, "STROMA_OTHER")) "lower_in_NEPC" else "higher_in_NEPC"

    cmp(d, "rho", paste0("spatial rho: ", pp), expected)
  }))

  if (!is.null(spatial_rho_tests)) {
    spatial_rho_tests$p_adj <- signif(p.adjust(spatial_rho_tests$p, "BH"), 3)

    cat("\nSpatial-bin rho tests, slide-level:\n")
    print(spatial_rho_tests, row.names = FALSE)

    write.csv(spatial_rho_tests, file.path(OUT, paste0(PFX, "spatial_bin_correlation_tests.csv")), row.names = FALSE)
  }

  p_rho <- spatial_rho |>
    filter(!skipped) |>
    mutate(pair = factor(pair, levels = RHO_PAIRS$label)) |>
    ggplot(aes(group, rho, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.5) +
    facet_wrap(~ pair, ncol = 3) +
    labs(
      title = paste0("Spatial-bin Spearman correlation, ", BIN_UM, " µm bins"),
      subtitle = "rho across occupied bins within each slide; pairs with <20 cells of either population skipped",
      x = NULL,
      y = "Spearman rho"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  p_rho <- p_rho + BOX_THEME
  n_rho <- length(unique(spatial_rho$pair[!spatial_rho$skipped]))
  save_pdf(p_rho, paste0(FPFX, "spatial_bin_rho.pdf"), fig_w(3), fig_h(ceiling(n_rho / 3)))
  keep_fig(p_rho, fig_w(3), fig_h(ceiling(n_rho / 3)))
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

  res <- data.frame(
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
  )

  if (nf < MIN_N || nt < MIN_N) return(res)

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

  res$obs_um <- round(obs, 1)
  res$null_um <- round(mu, 1)
  res$z <- round((obs - mu) / max(sdv, 1e-9), 2)
  res$frac_within_50um <- round(mean(dd <= 50), 3)
  res$frac_within_100um <- round(mean(dd <= 100), 3)
  res$p_perm <- signif((1 + sum(nullv <= obs, na.rm = TRUE)) / (N_PERM + 1), 3)
  res$skipped <- FALSE
  res
}

prox <- do.call(rbind, lapply(tier3_ok, function(sid) {
  d <- cell_df[cell_df$slide_id == sid, , drop = FALSE]
  grp <- as.character(d$group[1])

  do.call(rbind, lapply(STROMA_POPS, function(s)
    do.call(rbind, lapply(TUMOUR_POPS, function(t)
      prox_one(d, d[[s]], d[[t]], s, t, sid, grp)))))
}))

prox$group <- factor(prox$group, levels = GROUPS)

cat("\nProximity by slide:\n")
print(prox, row.names = FALSE)

write.csv(prox, file.path(OUT, paste0(PFX, "stroma_proximity_perslide.csv")), row.names = FALSE)

skip <- prox |>
  filter(skipped)

if (nrow(skip)) {
  warning(
    "Proximity skipped (<", MIN_N, " cells) for ",
    nrow(skip),
    " slide-pair combinations -- see skipped column.",
    call. = FALSE
  )
}

prox_tests <- do.call(rbind, lapply(STROMA_POPS, function(s)
  do.call(rbind, lapply(TUMOUR_POPS, function(tgt) {
    d <- prox |>
      filter(from == s, to == tgt, !skipped)

    if (!nrow(d)) return(NULL)

    lab <- paste0(s, "->", tgt)
    rbind(
      cmp(d, "obs_um", paste0(lab, ": median NN dist (um)"), "lower_in_NEPC"),
      cmp(d, "z", paste0(lab, ": proximity z"), "lower_in_NEPC"),
      cmp(d, "frac_within_50um", paste0(lab, ": frac within 50um"), "higher_in_NEPC"),
      cmp(d, "frac_within_100um", paste0(lab, ": frac within 100um"), "higher_in_NEPC")
    )
  }))))

if (!is.null(prox_tests)) {
  prox_tests$p_adj <- signif(p.adjust(prox_tests$p, "BH"), 3)

  cat("\nProximity: ADENO vs NEPC, slide-level\n")
  print(prox_tests, row.names = FALSE)

  write.csv(prox_tests, file.path(OUT, paste0(PFX, "stroma_proximity_group_tests.csv")), row.names = FALSE)
}

###############################################################################
## 4. NEIGHBOURHOOD COMPOSITION AROUND STROMAL CELLS
###############################################################################

rad_px <- RADIUS_UM / PX_UM
NBR_POPS <- c(POPS, "Multi", "Other")

nbr <- do.call(rbind, lapply(tier3_ok, function(sid) {
  d <- cell_df[cell_df$slide_id == sid, , drop = FALSE]
  xy <- as.matrix(d[, c("x_px", "y_px")])
  k <- min(100L, nrow(d))

  do.call(rbind, lapply(STROMA_POPS, function(s) {
    idx <- which(d[[s]] %in% TRUE)
    if (length(idx) < MIN_N) return(NULL)

    nn <- RANN::nn2(xy, xy[idx, , drop = FALSE], k = k)

    ## Neighbours within radius, excluding the centre cell itself (distance 0).
    within <- nn$nn.dists <= rad_px & nn$nn.dists > 0
    nb_pop <- d[[LABCOL]][nn$nn.idx[within]]
    cnt <- as.numeric(table(factor(nb_pop, levels = NBR_POPS)))

    obs <- cnt / max(sum(cnt), 1)
    exp <- as.numeric(table(factor(d[[LABCOL]], levels = NBR_POPS))) / nrow(d)

    data.frame(
      slide_id = sid,
      group = as.character(d$group[1]),
      centre = s,
      n_centre = length(idx),
      neighbour = NBR_POPS,
      obs_frac = round(obs, 4),
      exp_frac = round(exp, 4),
      log2_enrich = round(log2((obs + 1e-4) / (exp + 1e-4)), 2),
      row.names = NULL
    )
  }))
}))

if (!is.null(nbr) && nrow(nbr)) {
  nbr$group <- factor(nbr$group, levels = GROUPS)

  cat("\nNeighbourhood within ", RADIUS_UM, " um of ", paste(STROMA_POPS, collapse = " / "), " cells:\n", sep = "")
  print(nbr, row.names = FALSE)

  write.csv(nbr, file.path(OUT, paste0(PFX, "stroma_neighbourhood.csv")), row.names = FALSE)
}

###############################################################################
## 5. HYPOTHESIS SUMMARY
###############################################################################

pick <- function(tab, metric_name, prediction) {
  if (is.null(tab)) return(NULL)
  r <- tab[tab$metric == metric_name, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  r$prediction <- prediction
  r
}

hyp <- do.call(rbind, c(
  lapply(STROMA_POPS, function(s) rbind(
    pick(ab_tests, paste0(s, " per mm2 cell-bearing tile area"),
         paste0(s, " more abundant in NEPC")),
    pick(enrich_tests, paste0(s, " net fraction above blank"),
         paste0(s, " enriched above blank in NEPC"))
  )),
  lapply(TUMOUR_POPS, function(t) pick(ab_tests, paste0(t, " fraction of all cells"),
                                       paste0(t, " abundance, expected ", POP_DIRECTION[[t]]))),
  unlist(lapply(STROMA_POPS, function(s) lapply(TUMOUR_POPS, function(t) rbind(
    pick(prox_tests, paste0(s, "->", t, ": proximity z"),
         paste0(s, " closer to ", t, " in NEPC")),
    if (exists("spatial_rho_tests")) pick(spatial_rho_tests, paste0("spatial rho: ", s, " vs ", t),
         paste0(s, " co-localises with ", t, " more in NEPC"))
  ))), recursive = FALSE)
)) |>
  select(prediction, metric, n_ADENO, n_NEPC, med_ADENO, med_NEPC,
         hodges_lehmann, p, supports)

cat("\n================ HYPOTHESIS SUMMARY ================\n")
print(as.data.frame(hyp), row.names = FALSE)
cat(
  "`supports` = direction of effect only. With n<=3 per arm, p cannot\n",
  "reach 0.05; treat every row as descriptive.\n"
)

write.csv(hyp, file.path(OUT, paste0(PFX, "hypothesis_summary.csv")), row.names = FALSE)

###############################################################################
## 6. FIGURES
###############################################################################

message("\nWriting figures (PDF):")

p1 <- ab_slide |>
  select(slide_id, group, all_of(paste0(POPS, "_per_mm2"))) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(sub("_per_mm2$", "", metric), levels = POPS),
    group = factor(group, levels = GROUPS)
  ) |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(
    title = "Population abundance (cells per mm2 cell-bearing tile area)",
    subtitle = paste(
      paste(DEF_TEXT, "EPI_CA1 = AR+ FOLH1+ panCK+;"),
      "EPI_CA2 = AR+ FOLH1- panCK+; NE_CA = DLL3+"
    ),
    x = NULL,
    y = "cells / mm2"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p1 <- p1 + BOX_THEME
save_pdf(p1, paste0(FPFX, "population_abundance.pdf"), max(fig_w(length(POPS)), 11), fig_h(1))
keep_fig(p1, max(fig_w(length(POPS)), 11), fig_h(1))

p1s <- ab_slide |>
  select(slide_id, group, ends_with("_per_tumour"), ends_with("_per_epi_ca"),
         any_of("split_share"), EPI_CA1_frac_of_EPI_CA) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  mutate(group = factor(group, levels = GROUPS)) |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  labs(title = "Stromal ratios", subtitle = "each point = one slide",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

p1s <- p1s + BOX_THEME
n_p1s <- length(unique(p1s$data$metric))
save_pdf(p1s, paste0(FPFX, "stromal_ratios.pdf"), fig_w(ceiling(n_p1s / 2)), fig_h(2))
keep_fig(p1s, fig_w(ceiling(n_p1s / 2)), fig_h(2))

p1b <- blank_slide |>
  select(slide_id, group, starts_with("frac_")) |>
  pivot_longer(-c(slide_id, group), names_to = "metric", values_to = "value") |>
  mutate(
    control = ifelse(grepl("_blank$", metric), "blank", "marker"),
    population = factor(gsub("^frac_|_blank$", "", metric), levels = POPS),
    group = factor(group, levels = GROUPS)
  ) |>
  ggplot(aes(control, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.4) +
  facet_grid(population ~ group, scales = "free_y") +
  labs(
    title = "Blank false-positive controls",
    subtitle = "same population gates applied to the secondary-only channels",
    x = NULL,
    y = "fraction of cells"
  ) +
  theme_bw()

## Many rows: flatter panels so the page stays readable.
p1b <- p1b + theme(aspect.ratio = 0.6)
save_pdf(p1b, paste0(FPFX, "blank_false_positive_controls.pdf"), fig_w(2) + 1, fig_h(length(POPS), 0.6))
keep_fig(p1b, fig_w(2) + 1, fig_h(length(POPS), 0.6))

if (any(!prox$skipped)) {
  pd <- prox |>
    filter(!skipped) |>
    mutate(to = factor(to, levels = TUMOUR_POPS), from = factor(from, levels = STROMA_POPS))

  p2 <- ggplot(pd, aes(group, obs_um, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.6) +
    facet_grid(from ~ to, scales = "free_y") +
    labs(
      title = "Stromal -> tumour: median nearest-neighbour distance",
      subtitle = "slides with <20 cells of either population omitted",
      x = NULL,
      y = "median NN distance (µm)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  p2 <- p2 + BOX_THEME
  save_pdf(p2, paste0(FPFX, "proximity_distance.pdf"), fig_w(length(TUMOUR_POPS)), fig_h(length(STROMA_POPS)))
  keep_fig(p2, fig_w(length(TUMOUR_POPS)), fig_h(length(STROMA_POPS)))

  p3 <- ggplot(pd, aes(group, z, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.6) +
    facet_grid(from ~ to) +
    labs(
      title = "Proximity z vs label-permuted null",
      subtitle = "negative = closer than chance; density-corrected within slide",
      x = NULL,
      y = "z"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  p3 <- p3 + BOX_THEME
  save_pdf(p3, paste0(FPFX, "proximity_z.pdf"), fig_w(length(TUMOUR_POPS)), fig_h(length(STROMA_POPS)))
  keep_fig(p3, fig_w(length(TUMOUR_POPS)), fig_h(length(STROMA_POPS)))

  p4 <- pd |>
    select(slide_id, group, from, to, frac_within_50um, frac_within_100um) |>
    pivot_longer(starts_with("frac_"), names_to = "radius", values_to = "frac") |>
    mutate(radius = sub("frac_within_", "", radius)) |>
    ggplot(aes(group, frac, colour = group, shape = radius)) +
    geom_point(position = position_dodge(width = 0.5), size = 2.4) +
    facet_grid(from ~ to) +
    labs(
      title = "Fraction of stromal cells with a tumour cell within radius",
      x = NULL,
      y = "fraction"
    ) +
    theme_bw() +
    guides(colour = "none")

  p4 <- p4 + BOX_THEME
  save_pdf(p4, paste0(FPFX, "proximity_within_radius.pdf"), fig_w(length(TUMOUR_POPS)) + 0.8, fig_h(length(STROMA_POPS)))
  keep_fig(p4, fig_w(length(TUMOUR_POPS)) + 0.8, fig_h(length(STROMA_POPS)))
}

if (!is.null(nbr) && nrow(nbr)) {
  p5 <- nbr |>
    mutate(neighbour = factor(neighbour, levels = NBR_POPS),
           centre = factor(centre, levels = STROMA_POPS)) |>
    ggplot(aes(neighbour, slide_id, fill = log2_enrich)) +
    geom_tile() +
    geom_text(aes(label = log2_enrich), size = 2.6) +
    facet_grid(group ~ centre, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = paste0("Neighbourhood enrichment within ", RADIUS_UM, " µm of ", paste(STROMA_POPS, collapse = " / "), " cells"),
      subtitle = "log2(observed / slide-wide expected); centres with <20 cells omitted",
      x = NULL,
      y = NULL
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_pdf(p5, paste0(FPFX, "neighbourhood.pdf"), 11, 5)
  keep_fig(p5, 11, 5)
}


invisible(NULL)
}

run_set("PDGFRA", "pop", "", "04_", "PDGFRA = PDGFRA+;")
run_set(SPLIT_POPS, "pop_split", "exploratory_", "04x_exploratory_",
        "PDGFRA_CEBPDpos = PDGFRA+ CEBPD+; PDGFRA_CEBPDneg = PDGFRA+ CEBPD-;",
        CAPTION = "Exploratory: PDGFRA+ cells split by CEBPD")

## 7. SPATIAL MAPS
###############################################################################

map_dat <- cells |>
  filter(slide_id %in% tier3_ok)

if (!MAP_INCLUDE_OTHER) {
  map_dat <- map_dat |>
    filter(pop != "Other")
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
    pop = factor(pop, levels = c(POPS, "Multi", "Other"))
  )

pop_cols <- c(
  PDGFRA = "#FF7F00",
  EPI_CA1 = "#4DAF4A",
  EPI_CA2 = "#A6D854",
  NE_CA = "#377EB8",
  Multi = "#C77CFF",
  Other = "grey85"
)

if (nrow(map_dat)) {
  pmap <- ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
    geom_point(
      data = dplyr::filter(map_dat, pop == "Other"),
      size = 0.08,
      alpha = 0.18
    ) +
    geom_point(
      data = dplyr::filter(map_dat, pop != "Other"),
      size = 0.25,
      alpha = 0.65
    ) +
    facet_wrap(~ panel) +
    coord_equal() +
    scale_y_reverse() +
    scale_colour_manual(values = pop_cols, drop = FALSE) +
    labs(
      title = "Spatial distribution of PDGFRA+ stromal and tumour populations",
      subtitle = "all available slides; each panel translated to its own origin; shared mm scale",
      x = "mm",
      y = "mm"
    ) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  save_pdf(pmap, "04_spatial_maps.pdf", 12, 8)
  keep_fig(pmap, 12, 8)

  PDF_DEV(
    file.path(PLOTDIR, "04_spatial_maps_perslide.pdf"),
    width = 8,
    height = 8,
    onefile = TRUE
  )

  for (sid in unique(map_dat$slide_id)) {
    d <- map_dat |>
      filter(slide_id == sid)

    print(
      ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
        geom_point(
          data = dplyr::filter(d, pop == "Other"),
          size = 0.08,
          alpha = 0.18
        ) +
        geom_point(
          data = dplyr::filter(d, pop != "Other"),
          size = 0.25,
          alpha = 0.65
        ) +
        coord_equal() +
        scale_y_reverse() +
        scale_colour_manual(values = pop_cols, drop = FALSE) +
        labs(
          title = paste0(sid, "  [", unique(d$group), "]"),
          subtitle = paste(
            vapply(c(POPS, "Multi", "Other"), function(p)
              paste(p, sum(d$pop == p, na.rm = TRUE)), character(1)),
            collapse = " | "
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
  warning(
    "No cells available for spatial maps after slide filtering.",
    call. = FALSE
  )
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

  ## Draw each figure at its own size, top-centred, instead of stretching it
  ## to the page.
  for (f in .figs) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      x = 0.5, y = 1, just = c("centre", "top"),
      width = grid::unit(f$w, "in"), height = grid::unit(f$h, "in")
    ))
    print(f$plot, newpage = FALSE)
    grid::popViewport()
  }

  invisible(dev.off())

  message("  wrote 04_figures.pdf (", length(.figs), " pages)")
}

message("\n04_abundance_proximity.R complete -> ", PLOTDIR)
