###############################################################################
# 04_abundance_proximity.R  v11.7
#
# PURPOSE
#   Quantify ASPC-like stromal-cell abundance and spatial relationships to
#   epithelial and neuroendocrine tumour populations using outputs from
#   03_classify.R v11.4+.
#
# v11.7 CHANGE
#   - Removes slide-level marker-usability suppression from 04.
#   - Uses final cell-level calls written by 03_classify.R:
#        ASPC, EPI_CA, NE_CA
#     which are based on pixel-level corrected marker measurements.
#   - Keeps blank-derived calls:
#        ASPC_blank, EPI_CA_blank, NE_CA_blank
#     as audit/QC outputs only.
#   - Does not remove or suppress ASPC based on mean(PDGF_blank_pos).
#   - Keeps ASPC_raw / ASPC_not_callable columns if present for backwards
#     compatibility and audit, but final spatial/abundance analyses use ASPC.
#
# INPUT
#   results_ASPC_NEPC/cells_classified.rds
#   results_ASPC_NEPC/channels.rds
#
# OUTPUTS
#   CSV summaries in OUT
#   PDF figures in OUT/plots
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(RANN)
})

###############################################################################
## NAMESPACE SAFETY
###############################################################################

## Avoid namespace conflicts after sourcing 02_segment_measure.R, which loads MASS.
filter     <- dplyr::filter
select     <- dplyr::select
mutate     <- dplyr::mutate
summarise  <- dplyr::summarise
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
left_join  <- dplyr::left_join
transmute  <- dplyr::transmute
distinct   <- dplyr::distinct
count      <- dplyr::count
case_when  <- dplyr::case_when
any_of     <- dplyr::any_of
first      <- dplyr::first
bind_rows  <- dplyr::bind_rows

###############################################################################
## USER SETTINGS
###############################################################################

if (!exists("OUT")) OUT <- "/Volumes/T9/results_ASPC_NEPC"

PLOTDIR <- file.path(OUT, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

if (!exists("USE_ALL_SLIDES_FOR_ANALYSIS")) USE_ALL_SLIDES_FOR_ANALYSIS <- TRUE
if (!exists("MAP_INCLUDE_OTHER")) MAP_INCLUDE_OTHER <- TRUE

## Keep these display/audit categories if they exist.
if (!exists("MAP_SHOW_NOT_CALLABLE")) MAP_SHOW_NOT_CALLABLE <- TRUE

if (!exists("MIN_N")) MIN_N <- 20L
if (!exists("N_PERM")) N_PERM <- 199L

if (!exists("BIN_UM")) BIN_UM <- 250
if (!exists("MIN_BINS_FOR_RHO")) MIN_BINS_FOR_RHO <- 10L
if (!exists("MIN_TOTAL_CELLS_FOR_RHO")) MIN_TOTAL_CELLS_FOR_RHO <- 20L

if (!exists("RADIUS_UM")) RADIUS_UM <- 50

set.seed(1)

PDF_DEV <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

save_pdf <- function(p, file, width = 10, height = 6) {
  ok <- tryCatch({
    ggplot2::ggsave(
      filename = file.path(PLOTDIR, file),
      plot = p,
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

## Slide table is owned by 00_config.R.
.sf <- readRDS(file.path(OUT, "slides_final.rds"))
MANIFEST <- data.frame(
  slide_id = .sf$slide_id,
  group = as.character(.sf$group),
  tiff = .sf$tiff,
  stringsAsFactors = FALSE
)

GROUPS <- c("ADENO", "NEPC")

###############################################################################
## HELPERS
###############################################################################

.bool <- function(x) {
  y <- x %in% TRUE
  y[is.na(y)] <- FALSE
  y
}

.finite <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

.median_or_na <- function(x) {
  x <- .finite(x)
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}

.qsafe <- function(x, p, default = NA_real_) {
  x <- .finite(x)
  if (!length(x)) return(default)
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
}

.safe_sum <- function(x) sum(.bool(x), na.rm = TRUE)
.safe_mean <- function(x) mean(.bool(x), na.rm = TRUE)

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
    signif(suppressWarnings(stats::wilcox.test(A, N))$p.value, 3)
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

infer_tile_size <- function(d, default = 2048L) {
  if (!all(c("tile_x", "tile_y") %in% names(d))) return(default)

  tx <- sort(unique(d$tile_x))
  ty <- sort(unique(d$tile_y))
  step <- c(diff(tx)[diff(tx) > 0], diff(ty)[diff(ty) > 0])
  step <- step[is.finite(step)]

  if (!length(step)) return(default)
  as.integer(round(median(step)))
}

nn_med <- function(A, B) {
  if (!nrow(A) || !nrow(B)) return(NA_real_)
  median(RANN::nn2(B, A, k = 1)$nn.dists[, 1]) * PX_UM
}

###############################################################################
## LOAD INPUTS
###############################################################################

cells_file <- file.path(OUT, "cells_classified.rds")
channels_file <- file.path(OUT, "channels.rds")

if (!file.exists(cells_file)) {
  stop("Run 03_classify.R first: missing cells_classified.rds.", call. = FALSE)
}

if (!file.exists(channels_file)) {
  stop("Missing channels.rds. Re-run 01_channels.R.", call. = FALSE)
}

cells <- readRDS(cells_file)

if (!is.data.frame(cells) || !nrow(cells)) {
  stop("cells_classified.rds has zero rows or is not a data.frame.", call. = FALSE)
}

z <- readRDS(channels_file)
PX_UM <- if (!is.null(z$px_um)) z$px_um else 12500 / 44643

required_base <- c("slide_id", "x_px", "y_px", "area", "tile_x", "tile_y")
missing_base <- setdiff(required_base, names(cells))

if (length(missing_base)) {
  stop(
    "cells_classified.rds is missing required column(s): ",
    paste(missing_base, collapse = ", "),
    call. = FALSE
  )
}

present <- intersect(MANIFEST$slide_id, unique(cells$slide_id))
absent  <- setdiff(MANIFEST$slide_id, unique(cells$slide_id))
extra   <- setdiff(unique(cells$slide_id), MANIFEST$slide_id)

if (length(absent)) {
  warning(
    "No cells for: ", paste(absent, collapse = ", "),
    " -- excluded because they are absent from cells_classified.rds.",
    call. = FALSE
  )
}

if (length(extra)) {
  warning(
    "cells_classified.rds contains unexpected slide(s): ",
    paste(extra, collapse = ", "),
    " -- excluded from 04 analyses.",
    call. = FALSE
  )
}

if (!length(present)) {
  stop("No manifest slides are present in cells_classified.rds.", call. = FALSE)
}

if ("has_negctrl" %in% names(cells)) {
  .no_nc <- unique(cells$slide_id[!(cells$has_negctrl %in% TRUE)])
  .no_nc <- .no_nc[!is.na(.no_nc)]

  if (length(.no_nc)) {
    warning(
      "Slide(s) without negative controls in cells_classified.rds: ",
      paste(.no_nc, collapse = ", "),
      " -- blank-control audit metrics may be incomplete.",
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
## POPULATION COLUMNS
###############################################################################
## 04 now trusts the final pixel-level corrected calls from 03.
## Blank-derived calls are audit only.

required_pop <- c("ASPC", "EPI_CA", "NE_CA")
missing_pop <- setdiff(required_pop, names(cells))

if (length(missing_pop)) {
  stop(
    "cells_classified.rds is missing final population column(s): ",
    paste(missing_pop, collapse = ", "),
    "\nRe-run updated 03_classify.R first.",
    call. = FALSE
  )
}

cells$ASPC  <- .bool(cells$ASPC)
cells$EPI_CA <- .bool(cells$EPI_CA)
cells$NE_CA  <- .bool(cells$NE_CA)

## Backwards-compatible audit columns.
if (!"ASPC_raw" %in% names(cells)) cells$ASPC_raw <- cells$ASPC
if (!"EPI_CA_raw" %in% names(cells)) cells$EPI_CA_raw <- cells$EPI_CA
if (!"NE_CA_raw" %in% names(cells)) cells$NE_CA_raw <- cells$NE_CA

cells$ASPC_raw <- .bool(cells$ASPC_raw)
cells$EPI_CA_raw <- .bool(cells$EPI_CA_raw)
cells$NE_CA_raw <- .bool(cells$NE_CA_raw)

if (!"ASPC_not_callable" %in% names(cells)) cells$ASPC_not_callable <- FALSE
if (!"EPI_CA_not_callable" %in% names(cells)) cells$EPI_CA_not_callable <- FALSE
if (!"NE_CA_not_callable" %in% names(cells)) cells$NE_CA_not_callable <- FALSE

cells$ASPC_not_callable <- .bool(cells$ASPC_not_callable)
cells$EPI_CA_not_callable <- .bool(cells$EPI_CA_not_callable)
cells$NE_CA_not_callable <- .bool(cells$NE_CA_not_callable)

if (!"ASPC_CEBPg_pos" %in% names(cells)) cells$ASPC_CEBPg_pos <- FALSE
if (!"ASPC_CEBPg_neg" %in% names(cells)) cells$ASPC_CEBPg_neg <- cells$ASPC & !.bool(cells$ASPC_CEBPg_pos)

cells$ASPC_CEBPg_pos <- .bool(cells$ASPC_CEBPg_pos)
cells$ASPC_CEBPg_neg <- .bool(cells$ASPC_CEBPg_neg)

if (!"ASPC_blank" %in% names(cells)) {
  cells$ASPC_blank <- if ("PDGF_blank_pos" %in% names(cells)) .bool(cells$PDGF_blank_pos) else NA
}

if (!"EPI_CA_blank" %in% names(cells)) {
  if (all(c("panCK_blank_pos", "FOLH1_blank_pos", "AR_blank_pos") %in% names(cells))) {
    cells$EPI_CA_blank <- .bool(cells$panCK_blank_pos) &
      (.bool(cells$FOLH1_blank_pos) | .bool(cells$AR_blank_pos))
  } else {
    cells$EPI_CA_blank <- NA
  }
}

if (!"NE_CA_blank" %in% names(cells)) {
  cells$NE_CA_blank <- if ("DLL3_blank_pos" %in% names(cells)) .bool(cells$DLL3_blank_pos) else NA
}

## Rebuild display label from final calls.
cells <- cells |>
  mutate(
    n_pop_flags = rowSums(cbind(ASPC, EPI_CA, NE_CA), na.rm = TRUE),
    multi_pop = n_pop_flags > 1,
    pop = case_when(
      multi_pop ~ "Multi",
      ASPC ~ "ASPC",
      NE_CA ~ "NE_CA",
      EPI_CA ~ "EPI_CA",
      MAP_SHOW_NOT_CALLABLE & ASPC_not_callable ~ "ASPC_not_callable",
      MAP_SHOW_NOT_CALLABLE & EPI_CA_not_callable ~ "EPI_CA_not_callable",
      MAP_SHOW_NOT_CALLABLE & NE_CA_not_callable ~ "NE_CA_not_callable",
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

    n_ASPC = sum(.data$ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(.data$ASPC_not_callable, na.rm = TRUE),

    n_EPI_CA = sum(.data$EPI_CA, na.rm = TRUE),
    n_EPI_CA_raw = sum(.data$EPI_CA_raw, na.rm = TRUE),
    n_EPI_CA_not_callable = sum(.data$EPI_CA_not_callable, na.rm = TRUE),

    n_NE_CA = sum(.data$NE_CA, na.rm = TRUE),
    n_NE_CA_raw = sum(.data$NE_CA_raw, na.rm = TRUE),
    n_NE_CA_not_callable = sum(.data$NE_CA_not_callable, na.rm = TRUE),

    n_Multi = sum(.data$multi_pop, na.rm = TRUE),
    n_Other = sum(.data$pop == "Other", na.rm = TRUE),
    n_non_other = sum(.data$pop != "Other", na.rm = TRUE),

    frac_ASPC = mean(.data$ASPC, na.rm = TRUE),
    frac_EPI_CA = mean(.data$EPI_CA, na.rm = TRUE),
    frac_NE_CA = mean(.data$NE_CA, na.rm = TRUE),

    frac_ASPC_blank = mean(.data$ASPC_blank, na.rm = TRUE),
    frac_EPI_CA_blank = mean(.data$EPI_CA_blank, na.rm = TRUE),
    frac_NE_CA_blank = mean(.data$NE_CA_blank, na.rm = TRUE),

    .groups = "drop"
  )

cat("\nSlides analysed:\n")
print(as.data.frame(nslide))

cat("\nSlide inclusion mode:\n")
cat("  USE_ALL_SLIDES_FOR_ANALYSIS = ", USE_ALL_SLIDES_FOR_ANALYSIS, "\n", sep = "")
cat("  slides used for abundance/tier2 analyses: ", paste(tier2_ok, collapse = ", "), "\n")
cat("  slides used for spatial/tier3 analyses:   ", paste(tier3_ok, collapse = ", "), "\n")

cat("\nSlide QC / call summary:\n")
print(as.data.frame(slide_qc_summary))

write.csv(
  slide_qc_summary,
  file.path(OUT, "slide_qc_summary_v11_7_pixel_calls.csv"),
  row.names = FALSE
)

if (any(nslide$n < 2)) {
  warning("A group has <2 slides -- no test is possible; effect sizes only.", call. = FALSE)
}

###############################################################################
## AREA / TILE METRICS
###############################################################################

tile_area_mm2_by_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id) |>
  summarise(
    inferred_tile_size_px = infer_tile_size(dplyr::pick(dplyr::everything())),
    n_cellbearing_tiles = n_distinct(paste(tile_x, tile_y)),
    cellbearing_tile_mm2 =
      n_cellbearing_tiles * inferred_tile_size_px^2 * PX_UM^2 / 1e6,
    nuclear_area_mm2 = sum(area, na.rm = TRUE) * PX_UM^2 / 1e6,
    .groups = "drop"
  )

## Tissue area from the slide-level tissue mask written by 02 (preferred
## denominator for densities; cell-bearing tile area over-counts tissue at
## section edges).
tissue_area_file <- file.path(OUT, "tissue_area_by_slide.csv")
if (file.exists(tissue_area_file)) {
  tile_area_mm2_by_slide <- left_join(
    tile_area_mm2_by_slide,
    read.csv(tissue_area_file, stringsAsFactors = FALSE)[, c("slide_id", "tissue_mm2")],
    by = "slide_id"
  )
} else {
  tile_area_mm2_by_slide$tissue_mm2 <- tile_area_mm2_by_slide$cellbearing_tile_mm2
}

write.csv(
  tile_area_mm2_by_slide,
  file.path(OUT, "cellbearing_tile_area_by_slide.csv"),
  row.names = FALSE
)

###############################################################################
## 1. ASPC ABUNDANCE
###############################################################################

ab_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),

    n_ASPC = sum(.data$ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
    n_ASPC_not_callable = sum(.data$ASPC_not_callable, na.rm = TRUE),
    n_ASPC_CEBPg_pos = sum(.data$ASPC_CEBPg_pos, na.rm = TRUE),
    n_ASPC_CEBPg_neg = sum(.data$ASPC_CEBPg_neg, na.rm = TRUE),

    n_EPI_CA = sum(.data$EPI_CA, na.rm = TRUE),
    n_EPI_CA_raw = sum(.data$EPI_CA_raw, na.rm = TRUE),
    n_EPI_CA_not_callable = sum(.data$EPI_CA_not_callable, na.rm = TRUE),

    n_NE_CA = sum(.data$NE_CA, na.rm = TRUE),
    n_NE_CA_raw = sum(.data$NE_CA_raw, na.rm = TRUE),
    n_NE_CA_not_callable = sum(.data$NE_CA_not_callable, na.rm = TRUE),

    n_Multi = sum(.data$multi_pop, na.rm = TRUE),

    frac_ASPC = mean(.data$ASPC, na.rm = TRUE),
    frac_ASPC_raw = mean(.data$ASPC_raw, na.rm = TRUE),
    frac_ASPC_not_callable = mean(.data$ASPC_not_callable, na.rm = TRUE),
    frac_ASPC_CEBPg_pos = mean(.data$ASPC_CEBPg_pos, na.rm = TRUE),
    frac_ASPC_CEBPg_neg = mean(.data$ASPC_CEBPg_neg, na.rm = TRUE),

    frac_EPI_CA = mean(.data$EPI_CA, na.rm = TRUE),
    frac_EPI_CA_raw = mean(.data$EPI_CA_raw, na.rm = TRUE),
    frac_NE_CA = mean(.data$NE_CA, na.rm = TRUE),
    frac_NE_CA_raw = mean(.data$NE_CA_raw, na.rm = TRUE),
    frac_Multi = mean(.data$multi_pop, na.rm = TRUE),

    aspc_per_1000_cells = 1000 * mean(.data$ASPC, na.rm = TRUE),
    aspc_raw_per_1000_cells = 1000 * mean(.data$ASPC_raw, na.rm = TRUE),

    aspc_per_epi_ca = sum(.data$ASPC, na.rm = TRUE) /
      pmax(sum(.data$EPI_CA, na.rm = TRUE), 1),

    aspc_per_tumour = sum(.data$ASPC, na.rm = TRUE) /
      pmax(sum(.data$EPI_CA | .data$NE_CA, na.rm = TRUE), 1),

    .groups = "drop"
  ) |>
  left_join(tile_area_mm2_by_slide, by = "slide_id") |>
  mutate(
    aspc_per_mm2 = n_ASPC / pmax(tissue_mm2, 1e-9),
    aspc_raw_per_mm2 = n_ASPC_raw / pmax(tissue_mm2, 1e-9),
    epi_ca_per_mm2 = n_EPI_CA / pmax(tissue_mm2, 1e-9),
    ne_ca_per_mm2 = n_NE_CA / pmax(tissue_mm2, 1e-9)
  )

cat("\nAbundance by slide:\n")
print(as.data.frame(ab_slide))

write.csv(ab_slide, file.path(OUT, "aspc_abundance.csv"), row.names = FALSE)

ab_tests <- rbind(
  cmp(ab_slide, "frac_ASPC", "ASPC fraction of all cells", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_1000_cells", "ASPC per 1000 cells", "higher_in_NEPC"),
  cmp(ab_slide, "aspc_per_mm2", "ASPC per mm2 tissue", "higher_in_NEPC"),
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
## 1B. BLANK FALSE-POSITIVE POPULATION CONTROLS -- AUDIT ONLY
###############################################################################

blank_slide <- cells |>
  filter(slide_id %in% tier2_ok) |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),

    n_ASPC = sum(.data$ASPC, na.rm = TRUE),
    n_ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
    n_ASPC_blank = sum(.data$ASPC_blank, na.rm = TRUE),

    n_EPI_CA = sum(.data$EPI_CA, na.rm = TRUE),
    n_EPI_CA_blank = sum(.data$EPI_CA_blank, na.rm = TRUE),

    n_NE_CA = sum(.data$NE_CA, na.rm = TRUE),
    n_NE_CA_blank = sum(.data$NE_CA_blank, na.rm = TRUE),

    frac_ASPC = mean(.data$ASPC, na.rm = TRUE),
    frac_ASPC_raw = mean(.data$ASPC_raw, na.rm = TRUE),
    frac_ASPC_blank = mean(.data$ASPC_blank, na.rm = TRUE),

    frac_EPI_CA = mean(.data$EPI_CA, na.rm = TRUE),
    frac_EPI_CA_blank = mean(.data$EPI_CA_blank, na.rm = TRUE),

    frac_NE_CA = mean(.data$NE_CA, na.rm = TRUE),
    frac_NE_CA_blank = mean(.data$NE_CA_blank, na.rm = TRUE),

    ASPC_enrichment_over_blank =
      frac_ASPC / pmax(frac_ASPC_blank, 1e-9),
    ASPC_raw_enrichment_over_blank =
      frac_ASPC_raw / pmax(frac_ASPC_blank, 1e-9),
    EPI_CA_enrichment_over_blank =
      frac_EPI_CA / pmax(frac_EPI_CA_blank, 1e-9),
    NE_CA_enrichment_over_blank =
      frac_NE_CA / pmax(frac_NE_CA_blank, 1e-9),

    .groups = "drop"
  )

cat("\nBlank false-positive controls by slide; audit only, not used to suppress calls:\n")
print(as.data.frame(blank_slide))

write.csv(
  blank_slide,
  file.path(OUT, "aspc_blank_false_positive_controls_audit_only.csv"),
  row.names = FALSE
)

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

cat("\nBlank-adjusted enrichment audit by slide:\n")
print(as.data.frame(aspc_enrich_slide))

write.csv(
  aspc_enrich_slide,
  file.path(OUT, "aspc_blank_adjusted_enrichment_audit_only.csv"),
  row.names = FALSE
)

aspc_enrich_tests <- rbind(
  cmp(aspc_enrich_slide, "frac_ASPC_net", "ASPC net fraction above blank", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "ASPC_marker_to_blank", "ASPC marker-to-blank ratio", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "frac_EPI_CA_net", "EPI_CA net fraction above blank", "lower_in_NEPC"),
  cmp(aspc_enrich_slide, "EPI_CA_marker_to_blank", "EPI_CA marker-to-blank ratio", "lower_in_NEPC"),
  cmp(aspc_enrich_slide, "frac_NE_CA_net", "NE_CA net fraction above blank", "higher_in_NEPC"),
  cmp(aspc_enrich_slide, "NE_CA_marker_to_blank", "NE_CA marker-to-blank ratio", "higher_in_NEPC")
)

aspc_enrich_tests$p_adj <- signif(p.adjust(aspc_enrich_tests$p, "BH"), 3)

cat("\nBlank-adjusted enrichment tests; audit only:\n")
print(aspc_enrich_tests, row.names = FALSE)

write.csv(
  aspc_enrich_tests,
  file.path(OUT, "aspc_blank_adjusted_enrichment_tests_audit_only.csv"),
  row.names = FALSE
)

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
      ASPC = sum(.data$ASPC, na.rm = TRUE),
      ASPC_raw = sum(.data$ASPC_raw, na.rm = TRUE),
      ASPC_not_callable = sum(.data$ASPC_not_callable, na.rm = TRUE),
      ASPC_CEBPg_pos = sum(.data$ASPC_CEBPg_pos, na.rm = TRUE),
      ASPC_CEBPg_neg = sum(.data$ASPC_CEBPg_neg, na.rm = TRUE),
      EPI_CA = sum(.data$EPI_CA, na.rm = TRUE),
      NE_CA = sum(.data$NE_CA, na.rm = TRUE),
      Multi = sum(.data$multi_pop, na.rm = TRUE),
      STROMA_OTHER = sum(.data$STROMA_OTHER, na.rm = TRUE),
      .groups = "drop"
    )
}

spatial_bin_counts <- bind_rows(lapply(tier3_ok, function(sid) {
  d <- cells |> filter(slide_id == sid)
  make_bin_counts_one_slide(d, sid, bin_um = BIN_UM)
}))

spatial_rho <- NULL
spatial_rho_tests <- NULL

if (!nrow(spatial_bin_counts)) {
  warning("No spatial bins available for rho analysis.", call. = FALSE)
} else {
  write.csv(
    spatial_bin_counts,
    file.path(OUT, "spatial_bin_counts.csv"),
    row.names = FALSE
  )

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

    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))

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

  spatial_rho <- bind_rows(lapply(split(spatial_bin_counts, spatial_bin_counts$slide_id), function(b) {
    bind_rows(
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

  write.csv(
    spatial_rho,
    file.path(OUT, "spatial_bin_correlations.csv"),
    row.names = FALSE
  )

  spatial_rho_tests <- bind_rows(lapply(unique(spatial_rho$pair), function(pp) {
    d <- spatial_rho |> filter(pair == pp, !skipped)
    if (!nrow(d)) return(NULL)

    expected <- if (pp %in% c("ASPC vs EPI_CA", "ASPC vs NE_CA", "ASPC vs Multi")) {
      "higher_in_NEPC"
    } else {
      "lower_in_NEPC"
    }

    cmp(d, "rho", paste0("spatial rho: ", pp), expected)
  }))

  if (nrow(spatial_rho_tests)) {
    spatial_rho_tests$p_adj <- signif(p.adjust(spatial_rho_tests$p, "BH"), 3)

    cat("\nSpatial-bin rho tests, slide-level:\n")
    print(spatial_rho_tests, row.names = FALSE)

    write.csv(
      spatial_rho_tests,
      file.path(OUT, "spatial_bin_correlation_tests.csv"),
      row.names = FALSE
    )
  }
}

###############################################################################
## 3. SPATIAL PROXIMITY
###############################################################################

prox_one <- function(d, from_lgl, to_lgl, from_lab, to_lab, sid, grp) {
  xy <- as.matrix(d[, c("x_px", "y_px")])

  from_lgl <- .bool(from_lgl)
  to_lgl <- .bool(to_lgl)

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

prox <- bind_rows(lapply(tier3_ok, function(sid) {
  d <- cells |> filter(slide_id == sid)
  if (!nrow(d)) return(NULL)

  grp <- as.character(unique(d$group))[1]

  bind_rows(
    prox_one(d, d$ASPC, d$EPI_CA, "ASPC", "EPI_CA", sid, grp),
    prox_one(d, d$ASPC, d$NE_CA, "ASPC", "NE_CA", sid, grp)
  )
}))

prox_tests <- NULL

if (!nrow(prox)) {
  warning("No proximity results generated.", call. = FALSE)
} else {
  prox$group <- factor(prox$group, levels = GROUPS)

  cat("\nProximity by slide:\n")
  print(prox, row.names = FALSE)

  write.csv(
    prox,
    file.path(OUT, "aspc_proximity_perslide.csv"),
    row.names = FALSE
  )

  skip <- prox |> filter(skipped)
  if (nrow(skip)) {
    warning(
      "Proximity skipped (<", MIN_N, " cells) for ", nrow(skip),
      " slide-pair combinations -- see skipped column.",
      call. = FALSE
    )
  }

  prox_tests <- bind_rows(lapply(c("EPI_CA", "NE_CA"), function(tgt) {
    d <- prox |> filter(to == tgt, !skipped)
    if (!nrow(d)) return(NULL)

    bind_rows(
      cmp(d, "obs_um", paste0("ASPC->", tgt, ": median NN dist (um)"), "lower_in_NEPC"),
      cmp(d, "z", paste0("ASPC->", tgt, ": proximity z"), "lower_in_NEPC"),
      cmp(d, "frac_within_50um", paste0("ASPC->", tgt, ": frac within 50um"), "higher_in_NEPC"),
      cmp(d, "frac_within_100um", paste0("ASPC->", tgt, ": frac within 100um"), "higher_in_NEPC")
    )
  }))

  if (nrow(prox_tests)) {
    prox_tests$p_adj <- signif(p.adjust(prox_tests$p, "BH"), 3)

    cat("\nProximity: ADENO vs NEPC, slide-level\n")
    print(prox_tests, row.names = FALSE)

    write.csv(
      prox_tests,
      file.path(OUT, "aspc_proximity_group_tests.csv"),
      row.names = FALSE
    )
  }
}

###############################################################################
## 4. NEIGHBOURHOOD COMPOSITION AROUND ASPC
###############################################################################

## v11.7 capped the search at the 100 nearest cells (which truncates the
## 50 um disc in dense tumour) and looped over populations in R. Now: an exact
## radius search per slide, which also yields per-cell neighbourhood flags
## used by 05_visualize.R:
##   in_ASPC_nbhd        cell lies within RADIUS_UM of at least one ASPC
##                       (ASPC themselves included)
##   n_ASPC_within_r     number of ASPC (other than itself) within RADIUS_UM
## Composition is reported two ways:
##   pooled_pairs: all (ASPC, neighbour) pairs within the radius, as before;
##   unique_cells: distinct cells inside the union of ASPC discs, compared with
##                 whole-slide composition (the "ASPC neighbourhood" itself).

rad_px <- RADIUS_UM / PX_UM
if (!exists("NBR_K")) NBR_K <- 400L

NBR_POPS <- c("ASPC", "EPI_CA", "NE_CA", "Multi", "Other")

cells$in_ASPC_nbhd <- FALSE
cells$n_ASPC_within_r <- 0L

nbr_list <- list()
nbr_comp_list <- list()

for (sid in tier3_ok) {
  ii <- which(cells$slide_id == sid)
  d <- cells[ii, c("x_px", "y_px", "pop", "ASPC", "group")]
  ia <- which(d$ASPC)

  if (!length(ia)) next

  xy <- as.matrix(d[, c("x_px", "y_px")])
  k <- min(NBR_K, nrow(d))

  nn <- RANN::nn2(xy, xy[ia, , drop = FALSE], k = k,
                  searchtype = "radius", radius = rad_px)
  idx <- nn$nn.idx
  idx[idx == 0L | nn$nn.dists > rad_px] <- NA_integer_
  sat <- mean(!is.na(idx[, k]))
  if (sat > 0.01) {
    warning(sid, ": ", round(100 * sat, 1), "% of ASPC discs hit NBR_K = ", k,
            " neighbours; raise NBR_K.", call. = FALSE)
  }

  ## drop self
  idx[idx == matrix(ia, nrow(idx), ncol(idx))] <- NA_integer_

  nb <- as.vector(idx)
  nb <- nb[!is.na(nb)]

  hits <- tabulate(nb, nbins = nrow(d))
  in_nb <- hits > 0L
  in_nb[ia] <- TRUE

  cells$in_ASPC_nbhd[ii] <- in_nb
  cells$n_ASPC_within_r[ii] <- hits

  grp <- as.character(d$group[1])
  pop_f <- function(v) factor(v, levels = NBR_POPS)

  pair_cnt <- as.numeric(table(pop_f(d$pop[nb])))
  obs <- pair_cnt / max(sum(pair_cnt), 1)
  exp <- as.numeric(table(pop_f(d$pop))) / nrow(d)

  if (length(ia) >= MIN_N) {
    nbr_list[[sid]] <- data.frame(
      slide_id = sid,
      group = grp,
      n_ASPC = length(ia),
      neighbour = NBR_POPS,
      n_pairs = pair_cnt,
      obs_frac = round(obs, 4),
      exp_frac = round(exp, 4),
      log2_enrich = round(log2((obs + 1e-4) / (exp + 1e-4)), 2),
      row.names = NULL
    )
  }

  ## Cells within RADIUS_UM of at least one OTHER ASPC (other ASPC count, so
  ## ASPC-ASPC clustering is visible; a lone ASPC is not its own neighbour).
  uniq <- as.numeric(table(pop_f(d$pop[hits > 0L])))
  nbr_comp_list[[sid]] <- data.frame(
    slide_id = sid,
    group = grp,
    radius_um = RADIUS_UM,
    n_ASPC = length(ia),
    n_cells_slide = nrow(d),
    n_cells_in_nbhd = sum(in_nb),
    frac_slide_cells_in_nbhd = round(mean(in_nb), 4),
    population = NBR_POPS,
    n_in_nbhd = uniq,
    frac_in_nbhd = round(uniq / max(sum(uniq), 1), 4),
    frac_whole_slide = round(exp, 4),
    log2_enrich_vs_slide = round(log2((uniq / max(sum(uniq), 1) + 1e-4) / (exp + 1e-4)), 2),
    row.names = NULL
  )
}

nbr <- bind_rows(nbr_list)
nbr_comp <- bind_rows(nbr_comp_list)

if (nrow(nbr_comp)) {
  write.csv(nbr_comp, file.path(OUT, "aspc_neighbourhood_composition.csv"), row.names = FALSE)
  cat("\nASPC neighbourhood (union of ", RADIUS_UM, " um discs) vs whole slide:\n", sep = "")
  print(nbr_comp, row.names = FALSE)
}

## Slim per-cell table for visualisation / export.
cells_spatial <- cells[, intersect(c(
  "cell_id", "slide_id", "group", "x_px", "y_px", "area_um2", "pop",
  "ASPC", "ASPC_CEBPg_pos", "ASPC_CEBPg_neg", "EPI_CA", "NE_CA", "multi_pop",
  "in_ASPC_nbhd", "n_ASPC_within_r", "tile_x", "tile_y"
), names(cells))]
cells_spatial$px_um <- PX_UM
saveRDS(cells_spatial, file.path(OUT, "cells_spatial.rds"))

if (nrow(nbr)) {
  nbr$group <- factor(nbr$group, levels = GROUPS)

  cat("\nNeighbourhood within ", RADIUS_UM, " um of ASPC:\n", sep = "")
  print(nbr, row.names = FALSE)

  write.csv(
    nbr,
    file.path(OUT, "aspc_neighbourhood.csv"),
    row.names = FALSE
  )
} else {
  warning("No ASPC neighbourhood summaries generated; too few ASPC cells.", call. = FALSE)
}

###############################################################################
## 5. HYPOTHESIS SUMMARY
###############################################################################

hyp_parts <- list(
  ab_tests |>
    filter(metric == "ASPC per mm2 tissue") |>
    mutate(prediction = "(i) ASPC more abundant in NEPC"),

  aspc_enrich_tests |>
    filter(metric == "ASPC net fraction above blank") |>
    mutate(prediction = "(i-b, audit) ASPC above blank")
)

if (!is.null(prox_tests) && nrow(prox_tests)) {
  hyp_parts <- c(hyp_parts, list(
    prox_tests |>
      filter(metric == "ASPC->EPI_CA: proximity z") |>
      mutate(prediction = "(ii) ASPC closer to epithelial cancer in NEPC"),

    prox_tests |>
      filter(metric == "ASPC->NE_CA: proximity z") |>
      mutate(prediction = "(iii) ASPC closer to NE cancer in NEPC")
  ))
}

if (!is.null(spatial_rho_tests) && nrow(spatial_rho_tests)) {
  hyp_parts <- c(hyp_parts, list(
    spatial_rho_tests |>
      filter(metric == "spatial rho: ASPC vs EPI_CA") |>
      mutate(prediction = "(iv) ASPC spatial rho with epithelial cancer"),

    spatial_rho_tests |>
      filter(metric == "spatial rho: ASPC vs NE_CA") |>
      mutate(prediction = "(v) ASPC spatial rho with NE cancer"),

    spatial_rho_tests |>
      filter(metric == "spatial rho: STROMA_OTHER vs EPI_CA") |>
      mutate(prediction = "(vi) non-ASPC stromal spatial rho with epithelial cancer")
  ))
}

hyp <- bind_rows(hyp_parts) |>
  select(
    prediction, metric,
    n_ADENO, n_NEPC,
    med_ADENO, med_NEPC,
    hodges_lehmann, p, supports
  )

cat("\n================ HYPOTHESIS SUMMARY ================\n")
print(as.data.frame(hyp), row.names = FALSE)
cat("`supports` = direction of effect only. With n<=3 per arm, treat every row as descriptive.\n")
cat("Blank-derived metrics are audit/QC context only and do not suppress ASPC calls.\n")

write.csv(
  hyp,
  file.path(OUT, "hypothesis_summary.csv"),
  row.names = FALSE
)

###############################################################################
## 6. FIGURES
###############################################################################

message("\nWriting figures (PDF):")

p1 <- ab_slide |>
  select(
    slide_id, group,
    frac_ASPC, frac_ASPC_raw, frac_ASPC_not_callable,
    frac_ASPC_CEBPg_pos, frac_ASPC_CEBPg_neg,
    aspc_per_1000_cells, aspc_per_mm2,
    aspc_per_epi_ca, aspc_per_tumour
  ) |>
  tidyr::pivot_longer(
    -c(slide_id, group),
    names_to = "metric",
    values_to = "value"
  ) |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.6) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Prediction (i): ASPC abundance",
    subtitle = "ASPC uses pixel-level corrected PDGF-positive cell calls; blank metrics are audit only",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_pdf(p1, "04_aspc_abundance.pdf", 12, 6)
keep_fig(p1, 12, 6)

p1b <- blank_slide |>
  select(
    slide_id, group,
    frac_ASPC, frac_ASPC_raw, frac_ASPC_blank,
    frac_EPI_CA, frac_EPI_CA_blank,
    frac_NE_CA, frac_NE_CA_blank
  ) |>
  tidyr::pivot_longer(
    -c(slide_id, group),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    control = ifelse(grepl("_blank$", metric), "blank_as_marker", "marker"),
    population = gsub("^frac_|_blank$", "", metric)
  ) |>
  ggplot(aes(control, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.4) +
  facet_grid(population ~ group, scales = "free_y") +
  labs(
    title = "Blank false-positive controls",
    subtitle = "negative-control channels processed as pseudo-markers; audit only, not a hard exclusion gate",
    x = NULL,
    y = "fraction of cells"
  ) +
  theme_bw()

save_pdf(p1b, "04_blank_false_positive_controls_audit_only.pdf", 10, 7)
keep_fig(p1b, 10, 7)

p1c <- aspc_enrich_slide |>
  select(
    slide_id, group,
    frac_ASPC_net, ASPC_marker_to_blank,
    frac_EPI_CA_net, EPI_CA_marker_to_blank,
    frac_NE_CA_net, NE_CA_marker_to_blank
  ) |>
  tidyr::pivot_longer(
    -c(slide_id, group),
    names_to = "metric",
    values_to = "value"
  ) |>
  ggplot(aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 2.5) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Blank-adjusted marker enrichment audit",
    subtitle = "audit context only; final cell detection is pixel-level corrected",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

save_pdf(p1c, "04_aspc_blank_adjusted_enrichment_audit_only.pdf", 12, 6)
keep_fig(p1c, 12, 6)

if (!is.null(spatial_rho) && nrow(spatial_rho) && any(!spatial_rho$skipped)) {
  p_rho <- spatial_rho |>
    filter(!skipped) |>
    ggplot(aes(group, rho, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.5) +
    facet_wrap(~ pair, scales = "free_y") +
    labs(
      title = paste0("Spatial-bin Spearman correlation, ", BIN_UM, " um bins"),
      subtitle = "ASPC uses pixel-level corrected calls",
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
      title = paste0("Spatial-bin population co-occurrence, ", BIN_UM, " um bins"),
      subtitle = "each point = one occupied spatial bin; axes are log1p cell counts",
      x = "log1p target population count per bin",
      y = "log1p ASPC or STROMA_OTHER count per bin"
    ) +
    theme_bw()

  save_pdf(p_rho_scatter, "04_spatial_bin_scatter.pdf", 14, 9)
  keep_fig(p_rho_scatter, 14, 9)
}

if (!is.null(prox) && nrow(prox) && any(!prox$skipped)) {
  pd <- prox |> filter(!skipped)

  p2 <- ggplot(pd, aes(group, obs_um, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.6) +
    facet_wrap(~ to, scales = "free_y") +
    labs(
      title = "ASPC -> target: median nearest-neighbour distance",
      subtitle = "ASPC uses pixel-level corrected PDGF-positive calls",
      x = NULL,
      y = "median NN distance (um)"
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
    tidyr::pivot_longer(
      starts_with("frac_"),
      names_to = "radius",
      values_to = "frac"
    ) |>
    mutate(radius = sub("frac_within_", "", radius)) |>
    ggplot(aes(group, frac, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 2.4) +
    facet_grid(radius ~ to) +
    labs(
      title = "Fraction of ASPC with a target cell within radius",
      x = NULL,
      y = "fraction"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  save_pdf(p4, "04_proximity_within_radius.pdf", 9, 6)
  keep_fig(p4, 9, 6)
}

if (nrow(nbr)) {
  p5 <- ggplot(nbr, aes(neighbour, slide_id, fill = log2_enrich)) +
    geom_tile() +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = paste0("Neighbourhood enrichment within ", RADIUS_UM, " um of ASPC"),
      x = NULL,
      y = NULL
    ) +
    theme_bw()

  save_pdf(p5, "04_neighbourhood.pdf", 9, 5)
  keep_fig(p5, 9, 5)
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
    pop = factor(
      pop,
      levels = c(
        "ASPC", "EPI_CA", "NE_CA", "Multi",
        "ASPC_not_callable", "EPI_CA_not_callable", "NE_CA_not_callable",
        "Other"
      )
    )
  )

pop_cols <- c(
  ASPC = "#F8766D",
  EPI_CA = "#7CAE00",
  NE_CA = "#00BFC4",
  Multi = "#C77CFF",
  ASPC_not_callable = "grey35",
  EPI_CA_not_callable = "grey55",
  NE_CA_not_callable = "grey70",
  Other = "grey88"
)

if (nrow(map_dat)) {
  pmap <- ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
    geom_point(
      data = dplyr::filter(map_dat, pop == "Other"),
      size = 0.08,
      alpha = 0.18
    ) +
    geom_point(
      data = dplyr::filter(map_dat, grepl("_not_callable$", as.character(pop))),
      size = 0.10,
      alpha = 0.22
    ) +
    geom_point(
      data = dplyr::filter(map_dat, !pop %in% c(
        "Other", "ASPC_not_callable", "EPI_CA_not_callable", "NE_CA_not_callable"
      )),
      size = 0.25,
      alpha = 0.65
    ) +
    facet_wrap(~ panel) +
    coord_equal() +
    scale_y_reverse() +
    scale_colour_manual(values = pop_cols, drop = FALSE) +
    labs(
      title = "Spatial distribution of ASPC and tumour populations",
      subtitle = "ASPC/EPI_CA/NE_CA use pixel-level corrected calls; blank controls are audit only",
      x = "mm",
      y = "mm"
    ) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3)))

  ## Millions of points: write raster PNGs, not vector PDFs.
  ggplot2::ggsave(file.path(PLOTDIR, "04_spatial_maps.png"), pmap,
                  width = 12, height = 8, dpi = 200)
  message("  wrote 04_spatial_maps.png")

  for (sid in unique(map_dat$slide_id)) {
    d <- map_dat |> filter(slide_id == sid)

    ggplot2::ggsave(
      file.path(PLOTDIR, paste0("04_spatial_map_", sid, ".png")),
      width = 9, height = 9, dpi = 200,
      plot = ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
        geom_point(
          data = dplyr::filter(d, pop == "Other"),
          size = 0.08,
          alpha = 0.18
        ) +
        geom_point(
          data = dplyr::filter(d, grepl("_not_callable$", as.character(pop))),
          size = 0.10,
          alpha = 0.22
        ) +
        geom_point(
          data = dplyr::filter(d, !pop %in% c(
            "Other", "ASPC_not_callable", "EPI_CA_not_callable", "NE_CA_not_callable"
          )),
          size = 0.25,
          alpha = 0.65
        ) +
        coord_equal() +
        scale_y_reverse() +
        scale_colour_manual(values = pop_cols, drop = FALSE) +
        labs(
          title = paste0(sid, "  [", unique(d$group), "]"),
          subtitle = paste0(
            "ASPC ", sum(d$pop == "ASPC", na.rm = TRUE),
            " | EPI_CA ", sum(d$pop == "EPI_CA", na.rm = TRUE),
            " | NE_CA ", sum(d$pop == "NE_CA", na.rm = TRUE),
            " | Multi ", sum(d$pop == "Multi", na.rm = TRUE),
            " | ASPC_not_callable ", sum(d$pop == "ASPC_not_callable", na.rm = TRUE),
            " | Other ", sum(d$pop == "Other", na.rm = TRUE)
          ),
          x = "mm",
          y = "mm"
        ) +
        theme_bw() +
        guides(colour = guide_legend(override.aes = list(size = 3)))
    )
  }

  message("  wrote 04_spatial_map_<slide>.png")
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
