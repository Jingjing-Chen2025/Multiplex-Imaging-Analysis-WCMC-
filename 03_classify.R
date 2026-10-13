###############################################################################
# 03_classify.R  v11.2
#
# Classification using corrected features from 02_segment_measure.R v11.2.
#
# NOTE ON PMOT32 (v11.2):
#   This script reads cells_raw.rds and never touches the TIFFs directly, so
#   the PMOT32 re-acquisition (20260826, 8-cycle protocol with secondary-only
#   negative controls, AR/DLL3 = Cycle 6) requires no path changes here.
#   PMOT32 now carries the same blank-subtracted feature set as every other
#   slide, so the blank false-positive control below is defined for all six
#   slides and no slide-specific handling is needed.
#   Ensure 01_channels.R and 02_segment_measure.R have been re-run against the
#   new acquisition before running this script -- see the freshness guard below.
#
# Uses:
#   - local perinuclear ring features
#   - pixel-level blank subtraction
#   - AF puncta exclusion
#   - radius sensitivity, 4 / 6 / 8 µm
#   - blank false-positive control
#
# Primary classification radius:
#   RING_UM_PRIMARY = 6
#
# Population logic (2026-10-13, PDGFRA run):
#   PDGFRA  = PDGFRA_pos                       (main stromal population)
#   PDGFRA_CEBPDpos = PDGFRA_pos & CEBPD_pos   (exploratory split)
#   PDGFRA_CEBPDneg = PDGFRA_pos & !CEBPD_pos  (exploratory split)
#   EPI_CA1 = AR_pos & FOLH1_pos & panCK_pos
#   EPI_CA2 = AR_pos & !FOLH1_pos & panCK_pos
#   NE_CA   = DLL3_pos
###############################################################################

if (!exists("OUT")) OUT <- "/Volumes/T9/Case Studies/multiplex_analysis_20261013_PDGFRA_v2"

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
## INPUT
###############################################################################

infile <- file.path(OUT, "cells_raw.rds")

if (!file.exists(infile)) {
  stop("Missing cells_raw.rds. Re-run revised 02_segment_measure.R first.",
       call. = FALSE)
}

cells <- readRDS(infile)

if (!nrow(cells)) {
  stop("cells_raw.rds has zero rows.", call. = FALSE)
}

## Expected slide set. All six now share the 8-cycle protocol with NC cycles.
EXPECTED_SLIDES <- c("PM0T32", "PM677_A8", "PM12_B3",
                     "PM1024_B", "PM423", "PM159T5")

.present <- unique(cells$slide_id)
.absent  <- setdiff(EXPECTED_SLIDES, .present)
.extra   <- setdiff(.present, EXPECTED_SLIDES)

if (length(.absent)) {
  warning("cells_raw.rds is missing expected slide(s): ",
          paste(.absent, collapse = ", "),
          " -- was 02_segment_measure.R run on the full set?",
          call. = FALSE)
}

if (length(.extra)) {
  warning("cells_raw.rds contains unexpected slide(s): ",
          paste(.extra, collapse = ", "),
          call. = FALSE)
}

## Freshness guard: PMOT32 must come from the 20260826 re-acquisition, which
## carries negative controls. A stale row from the old 6-cycle run would have
## has_negctrl = FALSE and must not be classified with blank subtraction.
if ("has_negctrl" %in% names(cells)) {
  .no_nc <- unique(cells$slide_id[!cells$has_negctrl])
  if (length(.no_nc)) {
    stop("Slide(s) without negative controls in cells_raw.rds: ",
         paste(.no_nc, collapse = ", "),
         "\nThis indicates a stale segmentation. Re-run 01_channels.R and ",
         "02_segment_measure.R against the current acquisitions.",
         call. = FALSE)
  }
}

PX_UM <- 12500 / 44643
if ("area_um2" %in% names(cells)) {
  PX_UM <- sqrt(median(cells$area_um2 / cells$area, na.rm = TRUE))
}

###############################################################################
## QC FILTER
###############################################################################

MIN_AREA_UM2 <- 8
MAX_AREA_UM2 <- 400
MIN_REG_COR  <- 0.85

cells$area_um2 <- cells$area * PX_UM^2

keep <- is.finite(cells$area_um2) &
        cells$area_um2 >= MIN_AREA_UM2 &
        cells$area_um2 <= MAX_AREA_UM2 &
        (is.na(cells$reg_cor) | cells$reg_cor >= MIN_REG_COR)

cat("QC keep:", sum(keep), "/", nrow(cells), "cells\n")

cells <- cells[keep, , drop = FALSE]

if (!nrow(cells)) {
  stop("No cells remain after QC filtering.", call. = FALSE)
}

###############################################################################
## MARKER ROUTING
###############################################################################

RING_UM_PRIMARY <- 6
RING_UM_SET <- c(4, 6, 8)

MARKER_COMPARTMENT <- c(
  CEBPD = "nuc",
  AR    = "nuc",
  PDGFRA  = "ring",
  panCK = "ring",
  FOLH1 = "ring",
  DLL3  = "ring"
)

## Minimum fraction of the relevant compartment that must be positive.
## These protect against one bright residual AF dot calling a cell positive.
MIN_FRAC <- c(
  CEBPD = 0.08,
  AR    = 0.08,
  PDGFRA  = 0.06,
  panCK = 0.06,
  FOLH1 = 0.06,
  DLL3  = 0.05
)

## Strong morphology override:
## If enough compartment pixels are positive, allow low mean residual signal.
HIGH_FRAC <- c(
  CEBPD = 0.20,
  AR    = 0.20,
  PDGFRA  = 0.15,
  panCK = 0.15,
  FOLH1 = 0.15,
  DLL3  = 0.12
)

## Reject cells where too much of the relevant compartment was AF-masked.
MAX_AF_FRAC <- c(
  CEBPD = 0.05,
  AR    = 0.05,
  PDGFRA  = 0.10,
  panCK = 0.10,
  FOLH1 = 0.10,
  DLL3  = 0.10
)

.ring_tag <- function(r_um) {
  paste0("r", gsub("\\.", "p", as.character(r_um)), "um")
}

.feature_prefix <- function(marker, compartment, r_um = NULL, primary = RING_UM_PRIMARY) {
  if (compartment == "nuc") {
    return(paste0(marker, "_nuc"))
  }

  if (is.null(r_um) || identical(as.numeric(r_um), as.numeric(primary))) {
    return(paste0(marker, "_ring"))
  }

  paste0(marker, "_ring_", .ring_tag(r_um))
}

.required_cols_for <- function(marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  cp <- MARKER_COMPARTMENT[[marker]]
  pre <- .feature_prefix(marker, cp, r_um)

  if (blank) {
    c(
      paste0(pre, "_blank_sub"),
      paste0(pre, "_blank_frac_pos"),
      paste0(pre, "_blank_af_frac")
    )
  } else {
    c(
      paste0(pre, "_sub"),
      paste0(pre, "_frac_pos"),
      paste0(pre, "_af_frac")
    )
  }
}

all_required <- unique(unlist(lapply(names(MARKER_COMPARTMENT), function(m) {
  if (MARKER_COMPARTMENT[[m]] == "nuc") {
    c(
      .required_cols_for(m, blank = FALSE),
      .required_cols_for(m, blank = TRUE)
    )
  } else {
    unlist(lapply(RING_UM_SET, function(r) {
      c(
        .required_cols_for(m, r, blank = FALSE),
        .required_cols_for(m, r, blank = TRUE)
      )
    }))
  }
})))

missing <- setdiff(all_required, names(cells))

if (length(missing)) {
  stop(
    "Missing v11 corrected feature columns from 02:\n  ",
    paste(missing, collapse = "\n  "),
    "\nRe-run 02_segment_measure.R v11.2 before 03.",
    call. = FALSE
  )
}

###############################################################################
## HELPERS
###############################################################################

.qsafe <- function(x, p, default = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

.robust_thr <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 50) return(Inf)

  med <- median(x, na.rm = TRUE)
  sc  <- mad(x, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(sc) || sc <= 0) {
    sc <- sd(x, na.rm = TRUE)
  }

  if (!is.finite(sc) || sc <= 0) {
    return(.qsafe(x, 0.95, default = Inf))
  }

  med + 3 * sc
}

.get_feature_cols <- function(marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  cp <- MARKER_COMPARTMENT[[marker]]
  pre <- .feature_prefix(marker, cp, r_um)

  if (blank) {
    list(
      score = paste0(pre, "_blank_sub"),
      frac  = paste0(pre, "_blank_frac_pos"),
      af    = paste0(pre, "_blank_af_frac")
    )
  } else {
    list(
      score = paste0(pre, "_sub"),
      frac  = paste0(pre, "_frac_pos"),
      af    = paste0(pre, "_af_frac")
    )
  }
}

gate_marker_one_slide <- function(d, marker, r_um = RING_UM_PRIMARY, blank = FALSE) {
  cols <- .get_feature_cols(marker, r_um = r_um, blank = blank)

  score <- d[[cols$score]]
  frac  <- d[[cols$frac]]
  af    <- d[[cols$af]]

  eligible <- is.finite(score) &
              is.finite(frac) &
              is.finite(af) &
              af <= MAX_AF_FRAC[[marker]]

  if (sum(eligible) < 50) {
    pos <- rep(FALSE, nrow(d))

    return(list(
      pos = pos,
      thr = Inf,
      n_eligible = sum(eligible),
      n_pos = 0,
      frac_pos = 0,
      median_score = median(score, na.rm = TRUE),
      p95_score = .qsafe(score, 0.95),
      p99_score = .qsafe(score, 0.99),
      median_frac = median(frac, na.rm = TRUE),
      p95_frac = .qsafe(frac, 0.95),
      p99_frac = .qsafe(frac, 0.99),
      median_af_frac = median(af, na.rm = TRUE),
      p99_af_frac = .qsafe(af, 0.99)
    ))
  }

  thr <- .robust_thr(score[eligible])

  if (!is.finite(thr) || thr <= 0) {
    thr <- .qsafe(score[eligible], 0.95, default = Inf)
  }

  pos <- eligible &
         (
           (score > thr & frac >= MIN_FRAC[[marker]]) |
           (score > 0   & frac >= HIGH_FRAC[[marker]])
         )

  pos[is.na(pos)] <- FALSE

  list(
    pos = pos,
    thr = thr,
    n_eligible = sum(eligible),
    n_pos = sum(pos, na.rm = TRUE),
    frac_pos = mean(pos, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    p95_score = .qsafe(score, 0.95),
    p99_score = .qsafe(score, 0.99),
    median_frac = median(frac, na.rm = TRUE),
    p95_frac = .qsafe(frac, 0.95),
    p99_frac = .qsafe(frac, 0.99),
    median_af_frac = median(af, na.rm = TRUE),
    p99_af_frac = .qsafe(af, 0.99)
  )
}

###############################################################################
## PRIMARY CLASSIFICATION, 6 µm
###############################################################################

for (m in names(MARKER_COMPARTMENT)) {
  cells[[paste0(m, "_pos")]] <- FALSE
  cells[[paste0(m, "_blank_pos")]] <- FALSE
}

gates <- list()
blank_gates <- list()

for (sid in unique(cells$slide_id)) {
  ii <- which(cells$slide_id == sid)
  d  <- cells[ii, , drop = FALSE]

  for (m in names(MARKER_COMPARTMENT)) {
    g <- gate_marker_one_slide(
      d,
      marker = m,
      r_um = RING_UM_PRIMARY,
      blank = FALSE
    )

    cells[ii, paste0(m, "_pos")] <- g$pos

    gates[[length(gates) + 1]] <- data.frame(
      slide_id = sid,
      marker = m,
      radius_um = ifelse(MARKER_COMPARTMENT[[m]] == "nuc", NA, RING_UM_PRIMARY),
      compartment = MARKER_COMPARTMENT[[m]],
      control_type = "marker",
      threshold_sub = g$thr,
      n_eligible = g$n_eligible,
      n_pos = g$n_pos,
      frac_pos = g$frac_pos,
      median_sub = g$median_score,
      p95_sub = g$p95_score,
      p99_sub = g$p99_score,
      median_frac_pixels_pos = g$median_frac,
      p95_frac_pixels_pos = g$p95_frac,
      p99_frac_pixels_pos = g$p99_frac,
      median_af_frac = g$median_af_frac,
      p99_af_frac = g$p99_af_frac,
      min_frac = MIN_FRAC[[m]],
      high_frac = HIGH_FRAC[[m]],
      max_af_frac = MAX_AF_FRAC[[m]],
      stringsAsFactors = FALSE
    )

    bg <- gate_marker_one_slide(
      d,
      marker = m,
      r_um = RING_UM_PRIMARY,
      blank = TRUE
    )

    cells[ii, paste0(m, "_blank_pos")] <- bg$pos

    blank_gates[[length(blank_gates) + 1]] <- data.frame(
      slide_id = sid,
      marker = m,
      radius_um = ifelse(MARKER_COMPARTMENT[[m]] == "nuc", NA, RING_UM_PRIMARY),
      compartment = MARKER_COMPARTMENT[[m]],
      control_type = "blank_as_marker",
      threshold_sub = bg$thr,
      n_eligible = bg$n_eligible,
      n_pos = bg$n_pos,
      frac_pos = bg$frac_pos,
      median_sub = bg$median_score,
      p95_sub = bg$p95_score,
      p99_sub = bg$p99_score,
      median_frac_pixels_pos = bg$median_frac,
      p95_frac_pixels_pos = bg$p95_frac,
      p99_frac_pixels_pos = bg$p99_frac,
      median_af_frac = bg$median_af_frac,
      p99_af_frac = bg$p99_af_frac,
      min_frac = MIN_FRAC[[m]],
      high_frac = HIGH_FRAC[[m]],
      max_af_frac = MAX_AF_FRAC[[m]],
      stringsAsFactors = FALSE
    )
  }
}

qc_gates <- bind_rows(gates)
qc_blank_false_positive <- bind_rows(blank_gates)

write.csv(
  qc_gates,
  file.path(OUT, "qc_gates_primary_6um.csv"),
  row.names = FALSE
)

write.csv(
  qc_blank_false_positive,
  file.path(OUT, "qc_blank_false_positive_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Primary marker gates, 6 µm:\n")
print(qc_gates)

cat("\n=== Blank false-positive control, 6 µm:\n")
print(qc_blank_false_positive)

###############################################################################
## RADIUS SENSITIVITY, 4 / 6 / 8 µm
###############################################################################

sens <- list()

ring_markers <- names(MARKER_COMPARTMENT)[MARKER_COMPARTMENT == "ring"]

for (sid in unique(cells$slide_id)) {
  ii <- which(cells$slide_id == sid)
  d  <- cells[ii, , drop = FALSE]

  for (r_um in RING_UM_SET) {
    for (m in ring_markers) {
      g <- gate_marker_one_slide(
        d,
        marker = m,
        r_um = r_um,
        blank = FALSE
      )

      bg <- gate_marker_one_slide(
        d,
        marker = m,
        r_um = r_um,
        blank = TRUE
      )

      sens[[length(sens) + 1]] <- data.frame(
        slide_id = sid,
        marker = m,
        radius_um = r_um,
        marker_frac_pos = g$frac_pos,
        marker_n_pos = g$n_pos,
        marker_n_eligible = g$n_eligible,
        blank_frac_pos = bg$frac_pos,
        blank_n_pos = bg$n_pos,
        blank_n_eligible = bg$n_eligible,
        enrichment_over_blank = ifelse(bg$frac_pos > 0, g$frac_pos / bg$frac_pos, Inf),
        marker_threshold_sub = g$thr,
        blank_threshold_sub = bg$thr,
        marker_p99_sub = g$p99_score,
        blank_p99_sub = bg$p99_score,
        stringsAsFactors = FALSE
      )
    }
  }
}

qc_radius_sensitivity <- bind_rows(sens)

write.csv(
  qc_radius_sensitivity,
  file.path(OUT, "qc_radius_sensitivity_4_6_8um.csv"),
  row.names = FALSE
)

cat("\n=== Radius sensitivity, ring markers:\n")
print(qc_radius_sensitivity)

###############################################################################
## POPULATION CALLS
###############################################################################

cells$CEBPD_neg <- !cells$CEBPD_pos

## Population definitions (2026-10-13).
cells$PDGFRA          <- cells$PDGFRA_pos
cells$PDGFRA_CEBPDpos <- cells$PDGFRA_pos & cells$CEBPD_pos
cells$PDGFRA_CEBPDneg <- cells$PDGFRA_pos & !cells$CEBPD_pos
cells$EPI_CA1 <- cells$AR_pos & cells$FOLH1_pos & cells$panCK_pos
cells$EPI_CA2 <- cells$AR_pos & !cells$FOLH1_pos & cells$panCK_pos
cells$NE_CA   <- cells$DLL3_pos

## Blank-derived false-positive population controls, same gates on blank channels.
cells$PDGFRA_blank          <- cells$PDGFRA_blank_pos
cells$PDGFRA_CEBPDpos_blank <- cells$PDGFRA_blank_pos & cells$CEBPD_blank_pos
cells$PDGFRA_CEBPDneg_blank <- cells$PDGFRA_blank_pos & !cells$CEBPD_blank_pos
cells$EPI_CA1_blank <- cells$AR_blank_pos & cells$FOLH1_blank_pos & cells$panCK_blank_pos
cells$EPI_CA2_blank <- cells$AR_blank_pos & !cells$FOLH1_blank_pos & cells$panCK_blank_pos
cells$NE_CA_blank   <- cells$DLL3_blank_pos

## POPS drive the display label and overlap (Multi) calls; the CEBPD split is
## a subset of PDGFRA and is counted separately below.
POPS <- c("PDGFRA", "EPI_CA1", "EPI_CA2", "NE_CA")
SPLIT_POPS <- c("PDGFRA_CEBPDpos", "PDGFRA_CEBPDneg")

## Display label.
## Boolean flags above are authoritative; this label is for simple plotting/tables.
cells$pop <- "Other"

for (p in POPS) cells$pop[cells[[p]] %in% TRUE] <- p

cells$n_pop_flags <- rowSums(as.matrix(cells[, POPS]), na.rm = TRUE)

cells$multi_pop <- cells$n_pop_flags > 1
cells$pop[cells$multi_pop] <- "Multi"

###############################################################################
## OUTPUTS
###############################################################################

saveRDS(cells, file.path(OUT, "cells_classified.rds"))

cat("\n=== Cells per slide after QC:\n")
print(table(cells$slide_id))

cat("\n=== Single-marker positive fraction by slide, primary 6 µm:\n")

marker_frac <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_pos")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(marker_frac, 4))

write.csv(
  data.frame(slide_id = rownames(marker_frac), marker_frac, row.names = NULL),
  file.path(OUT, "marker_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Single-marker blank false-positive fraction by slide, primary 6 µm:\n")

blank_frac <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_blank_pos")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(blank_frac, 4))

write.csv(
  data.frame(slide_id = rownames(blank_frac), blank_frac, row.names = NULL),
  file.path(OUT, "marker_blank_false_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

population_counts <- data.frame(
  slide_id = names(table(cells$slide_id)),
  n_cells = as.integer(table(cells$slide_id)),

  sapply(c(POPS, SPLIT_POPS), function(p) as.integer(tapply(cells[[p]], cells$slide_id, sum, na.rm = TRUE))),
  Other = as.integer(tapply(cells$pop == "Other", cells$slide_id, sum, na.rm = TRUE)),
  Multi = as.integer(tapply(cells$multi_pop, cells$slide_id, sum, na.rm = TRUE)),
  setNames(
    as.data.frame(sapply(c(POPS, SPLIT_POPS), function(p)
      as.integer(tapply(cells[[paste0(p, "_blank")]], cells$slide_id, sum, na.rm = TRUE)))),
    paste0(c(POPS, SPLIT_POPS), "_blank")
  ),

  row.names = NULL
)

cat("\n=== Population counts by slide, primary 6 µm:\n")
print(population_counts)

write.csv(
  population_counts,
  file.path(OUT, "population_counts_by_slide_primary_6um.csv"),
  row.names = FALSE
)

population_fractions <- population_counts

for (cl in setdiff(names(population_fractions), c("slide_id", "n_cells"))) {
  population_fractions[[cl]] <- population_fractions[[cl]] / population_fractions$n_cells
}

cat("\n=== Population fractions by slide, primary 6 µm:\n")
print(population_fractions)

write.csv(
  population_fractions,
  file.path(OUT, "population_fractions_by_slide_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Overlapping population flags:\n")
print(table(cells$n_pop_flags))

message("\n03_classify.R complete -> ", file.path(OUT, "cells_classified.rds"))
