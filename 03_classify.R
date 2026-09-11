###############################################################################
# 03_classify.R  v11.3
#
# Classification using corrected features from 02_segment_measure.R v11.3+.
#
# NOTE ON PM0T32 (v11.2+):
#   This script reads cells_raw.rds and never touches the TIFFs directly, so
#   the PM0T32 re-acquisition (20260826, 8-cycle protocol with secondary-only
#   negative controls, AR/DLL3 = Cycle 6) requires no path changes here.
#   PM0T32 now carries the same blank-subtracted feature set as every other
#   slide, so the blank false-positive control below is defined for all six
#   slides and no slide-specific handling is needed.
#   Ensure 01_channels.R and 02_segment_measure.R have been re-run against the
#   new acquisition before running this script -- see the freshness guard below.
#
# Uses:
#   - local perinuclear ring features
#   - pixel-level scaled blank correction from 02_segment_measure.R
#   - AF puncta exclusion
#   - radius sensitivity, 4 / 6 / 8 um
#   - blank false-positive control
#   - marker usability QC before population assignment
#
# Primary classification radius:
#   RING_UM_PRIMARY = 6
#
# Population logic:
#   ASPC_raw = (CEBPg_pos & PDGF_pos) | (CEBPg_neg & PDGF_pos)
#            = PDGF_pos, with CEBPg status retained
#
#   ASPC     = ASPC_raw only on slides where PDGF passes usability QC
#
#   EPI_CA   = (panCK_pos & FOLH1_pos) | (panCK_pos & AR_pos)
#
#   NE_CA    = DLL3_pos
###############################################################################

if (!exists("OUT")) OUT <- "results_ASPC_NEPC"

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

## Freshness guard: every classified slide must carry negative controls. A stale
## row from an old run would have has_negctrl = FALSE or NA and must not be
## classified with blank subtraction.
if ("has_negctrl" %in% names(cells)) {
  .no_nc <- unique(cells$slide_id[!(cells$has_negctrl %in% TRUE)])
  if (length(.no_nc)) {
    stop("Slide(s) without negative controls in cells_raw.rds: ",
         paste(.no_nc, collapse = ", "),
         "\nThis indicates a stale segmentation. Re-run 01_channels.R and ",
         "02_segment_measure.R against the current acquisitions.",
         call. = FALSE)
  }
} else {
  stop("cells_raw.rds is missing has_negctrl; re-run revised 02_segment_measure.R first.",
       call. = FALSE)
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

## v11.3 registration QC treats empty/flat windows as unknown, not failed.
## Keep unknown/non-informative tiles; only drop informative tiles that fail.
if (!"reg_qc_informative" %in% names(cells)) {
  cells$reg_qc_informative <- is.finite(cells$reg_cor)
}

keep <- is.finite(cells$area_um2) &
        cells$area_um2 >= MIN_AREA_UM2 &
        cells$area_um2 <= MAX_AREA_UM2 &
        (
          !(cells$reg_qc_informative %in% TRUE) |
            is.na(cells$reg_cor) |
            cells$reg_cor >= MIN_REG_COR
        )

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
  CEBPg = "nuc",
  AR    = "nuc",
  PDGF  = "ring",
  panCK = "ring",
  FOLH1 = "ring",
  DLL3  = "ring"
)

## Minimum fraction of the relevant compartment that must be positive.
## These protect against one bright residual AF dot calling a cell positive.
MIN_FRAC <- c(
  CEBPg = 0.08,
  AR    = 0.08,
  PDGF  = 0.06,
  panCK = 0.06,
  FOLH1 = 0.06,
  DLL3  = 0.05
)

## Retired v11.2 strong-fraction override.
## The old gate allowed `score > 0 & frac >= HIGH_FRAC`, which is unsafe when
## marker and blank differ by a small acquisition gain. A PM0T32-like positive
## residual pedestal can make many cells positive even without marker-specific
## signal. v11.3 keeps this object only for backwards-compatible QC output; it
## is not used in the positivity gate.
HIGH_FRAC <- c(
  CEBPg = NA_real_,
  AR    = NA_real_,
  PDGF  = NA_real_,
  panCK = NA_real_,
  FOLH1 = NA_real_,
  DLL3  = NA_real_
)

## Reject cells where too much of the relevant compartment was AF-masked.
MAX_AF_FRAC <- c(
  CEBPg = 0.05,
  AR    = 0.05,
  PDGF  = 0.10,
  panCK = 0.10,
  FOLH1 = 0.10,
  DLL3  = 0.10
)

## Marker usability criteria, applied after marker and blank-as-marker gates.
## These prevent a technically failed marker channel from driving populations.
MIN_MARKER_TO_BLANK_ENRICHMENT <- 2.0
MAX_BLANK_FRAC_POS <- 0.10
MAX_GAIN_R2_FOR_WEAK_ENRICHMENT <- 0.90
GAIN_SLOPE_TOL <- 0.05

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
    "\nRe-run 02_segment_measure.R v11.3+ before 03.",
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

.median_col <- function(d, col, default = NA_real_) {
  if (!col %in% names(d)) return(default)
  x <- d[[col]]
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  median(x, na.rm = TRUE)
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

  ## Strict v11.3 gate.
  ## Removed the old rescue clause `(score > 0 & frac >= HIGH_FRAC[[marker]])`
  ## because it converts a weak positive residual pedestal into false positives.
  pos <- eligible &
         score > thr &
         frac >= MIN_FRAC[[marker]]

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
## PRIMARY CLASSIFICATION, 6 um
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

cat("\n=== Primary marker gates, 6 um:\n")
print(qc_gates)

cat("\n=== Blank false-positive control, 6 um:\n")
print(qc_blank_false_positive)

###############################################################################
## MARKER USABILITY QC, 6 um
###############################################################################

marker_usability <- do.call(rbind, lapply(unique(cells$slide_id), function(sid) {
  d <- cells[cells$slide_id == sid, , drop = FALSE]

  do.call(rbind, lapply(names(MARKER_COMPARTMENT), function(m) {
    cp <- MARKER_COMPARTMENT[[m]]
    pre <- .feature_prefix(m, cp, RING_UM_PRIMARY)

    marker_row <- qc_gates[
      qc_gates$slide_id == sid &
        qc_gates$marker == m &
        qc_gates$control_type == "marker",
      , drop = FALSE
    ]

    blank_row <- qc_blank_false_positive[
      qc_blank_false_positive$slide_id == sid &
        qc_blank_false_positive$marker == m &
        qc_blank_false_positive$control_type == "blank_as_marker",
      , drop = FALSE
    ]

    marker_frac <- if (nrow(marker_row)) marker_row$frac_pos[1] else NA_real_
    blank_frac  <- if (nrow(blank_row))  blank_row$frac_pos[1]  else NA_real_

    enrichment <- marker_frac / pmax(blank_frac, 1e-9)

    slope <- .median_col(d, paste0(pre, "_blank_slope"))
    r2    <- .median_col(d, paste0(pre, "_blank_r2"))

    fail_blank <- is.finite(blank_frac) && blank_frac > MAX_BLANK_FRAC_POS
    fail_enrich <- !is.finite(enrichment) ||
                   enrichment < MIN_MARKER_TO_BLANK_ENRICHMENT

    ## Catch channels where the marker is highly explained by the matched blank
    ## and marker-vs-blank enrichment is weak. This is the PM0T32 PDGF/PDGFRa
    ## failure mode, but the rule is data-driven and applies to all slides.
    fail_gain <- is.finite(r2) &&
                 r2 >= MAX_GAIN_R2_FOR_WEAK_ENRICHMENT &&
                 is.finite(slope) &&
                 abs(slope - 1) > GAIN_SLOPE_TOL &&
                 fail_enrich

    marker_usable <- !(fail_blank | fail_enrich | fail_gain)

    fail_reason <- if (marker_usable) {
      "pass"
    } else {
      trimws(paste(
        if (fail_blank) "high_blank_FP" else "",
        if (fail_enrich) "low_marker_to_blank_enrichment" else "",
        if (fail_gain) "high_marker_blank_correlation_gain_mismatch" else ""
      ))
    }

    data.frame(
      slide_id = sid,
      marker = m,
      marker_frac_pos = marker_frac,
      blank_frac_pos = blank_frac,
      enrichment_over_blank = enrichment,
      median_blank_slope = slope,
      median_blank_r2 = r2,
      marker_usable = marker_usable,
      fail_reason = fail_reason,
      stringsAsFactors = FALSE
    )
  }))
}))

cat("\n=== Marker usability, primary 6 um:\n")
print(marker_usability, row.names = FALSE)

write.csv(
  marker_usability,
  file.path(OUT, "qc_marker_usability_primary_6um.csv"),
  row.names = FALSE
)

## Add one usability column per marker onto cells for downstream scripts.
for (m in names(MARKER_COMPARTMENT)) {
  usable_col <- paste0(m, "_usable")
  u <- marker_usability[marker_usability$marker == m,
                        c("slide_id", "marker_usable"),
                        drop = FALSE]
  names(u)[2] <- usable_col

  cells <- dplyr::left_join(cells, u, by = "slide_id")
  cells[[usable_col]][is.na(cells[[usable_col]])] <- FALSE
}

###############################################################################
## RADIUS SENSITIVITY, 4 / 6 / 8 um
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

cells$CEBPg_neg <- !cells$CEBPg_pos

## Raw requested ASPC-like definition.
## This simplifies to PDGF_pos, but is retained for traceability.
cells$ASPC_raw <- (cells$CEBPg_pos & cells$PDGF_pos) |
                  (cells$CEBPg_neg & cells$PDGF_pos)

## PDGF/PDGFRa must be usable on the slide before ASPC can be asserted.
## If not callable, raw calls are retained in ASPC_raw for audit but excluded
## from ASPC and downstream abundance/spatial analyses.
cells$ASPC_callable <- cells$PDGF_usable %in% TRUE
cells$ASPC <- cells$ASPC_callable & cells$ASPC_raw

cells$ASPC_CEBPg_pos <- cells$ASPC & cells$CEBPg_pos
cells$ASPC_CEBPg_neg <- cells$ASPC & !cells$CEBPg_pos
cells$ASPC_raw_CEBPg_pos <- cells$ASPC_raw & cells$CEBPg_pos
cells$ASPC_raw_CEBPg_neg <- cells$ASPC_raw & !cells$CEBPg_pos

if (any(!cells$ASPC_callable, na.rm = TRUE)) {
  bad_aspc <- unique(cells$slide_id[!cells$ASPC_callable])
  warning(
    "ASPC calls are not callable on slide(s): ",
    paste(bad_aspc, collapse = ", "),
    ". Raw PDGF-derived calls are retained in ASPC_raw but excluded from ASPC.",
    call. = FALSE
  )
}

cells$EPI_CA <- (cells$panCK_pos & cells$FOLH1_pos) |
                (cells$panCK_pos & cells$AR_pos)

cells$NE_CA <- cells$DLL3_pos

## Blank-derived false-positive population controls.
cells$ASPC_blank <- cells$PDGF_blank_pos

cells$EPI_CA_blank <- (cells$panCK_blank_pos & cells$FOLH1_blank_pos) |
                      (cells$panCK_blank_pos & cells$AR_blank_pos)

cells$NE_CA_blank <- cells$DLL3_blank_pos

## Display label.
## Boolean flags above are authoritative; this label is for simple plotting/tables.
cells$pop <- "Other"

cells$pop[cells$ASPC]   <- "ASPC"
cells$pop[cells$EPI_CA] <- "EPI_CA"
cells$pop[cells$NE_CA]  <- "NE_CA"

cells$n_pop_flags <- rowSums(cbind(
  ASPC = cells$ASPC,
  EPI_CA = cells$EPI_CA,
  NE_CA = cells$NE_CA
), na.rm = TRUE)

cells$multi_pop <- cells$n_pop_flags > 1

###############################################################################
## OUTPUTS
###############################################################################

saveRDS(cells, file.path(OUT, "cells_classified.rds"))

cat("\n=== Cells per slide after QC:\n")
print(table(cells$slide_id))

cat("\n=== Single-marker positive fraction by slide, primary 6 um:\n")

marker_frac <- sapply(names(MARKER_COMPARTMENT), function(m) {
  tapply(cells[[paste0(m, "_pos")]], cells$slide_id, mean, na.rm = TRUE)
})

print(round(marker_frac, 4))

write.csv(
  data.frame(slide_id = rownames(marker_frac), marker_frac, row.names = NULL),
  file.path(OUT, "marker_positive_fraction_by_slide_primary_6um.csv"),
  row.names = FALSE
)

cat("\n=== Single-marker blank false-positive fraction by slide, primary 6 um:\n")

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

  ASPC = as.integer(tapply(cells$ASPC, cells$slide_id, sum, na.rm = TRUE)),
  ASPC_raw = as.integer(tapply(cells$ASPC_raw, cells$slide_id, sum, na.rm = TRUE)),
  ASPC_not_callable = as.integer(tapply(
    cells$ASPC_raw & !cells$ASPC_callable,
    cells$slide_id,
    sum,
    na.rm = TRUE
  )),
  ASPC_callable = as.logical(tapply(
    cells$ASPC_callable,
    cells$slide_id,
    function(x) all(x %in% TRUE)
  )),
  ASPC_CEBPg_pos = as.integer(tapply(cells$ASPC_CEBPg_pos, cells$slide_id, sum, na.rm = TRUE)),
  ASPC_CEBPg_neg = as.integer(tapply(cells$ASPC_CEBPg_neg, cells$slide_id, sum, na.rm = TRUE)),
  ASPC_raw_CEBPg_pos = as.integer(tapply(cells$ASPC_raw_CEBPg_pos, cells$slide_id, sum, na.rm = TRUE)),
  ASPC_raw_CEBPg_neg = as.integer(tapply(cells$ASPC_raw_CEBPg_neg, cells$slide_id, sum, na.rm = TRUE)),

  EPI_CA = as.integer(tapply(cells$EPI_CA, cells$slide_id, sum, na.rm = TRUE)),
  NE_CA = as.integer(tapply(cells$NE_CA, cells$slide_id, sum, na.rm = TRUE)),
  Other = as.integer(tapply(cells$pop == "Other", cells$slide_id, sum, na.rm = TRUE)),
  Multi = as.integer(tapply(cells$multi_pop, cells$slide_id, sum, na.rm = TRUE)),

  ASPC_blank = as.integer(tapply(cells$ASPC_blank, cells$slide_id, sum, na.rm = TRUE)),
  EPI_CA_blank = as.integer(tapply(cells$EPI_CA_blank, cells$slide_id, sum, na.rm = TRUE)),
  NE_CA_blank = as.integer(tapply(cells$NE_CA_blank, cells$slide_id, sum, na.rm = TRUE)),

  row.names = NULL
)

cat("\n=== Population counts by slide, primary 6 um:\n")
print(population_counts)

write.csv(
  population_counts,
  file.path(OUT, "population_counts_by_slide_primary_6um.csv"),
  row.names = FALSE
)

population_fractions <- population_counts

for (cl in setdiff(names(population_fractions), c("slide_id", "n_cells", "ASPC_callable"))) {
  population_fractions[[cl]] <- population_fractions[[cl]] / population_fractions$n_cells
}

cat("\n=== Population fractions by slide, primary 6 um:\n")
print(population_fractions)

write.csv(
  population_fractions,
  file.path(OUT, "population_fractions_by_slide_primary_6um.csv"),
  row.names = FALSE
)

aspc_callable_summary <- cells |>
  group_by(slide_id) |>
  summarise(
    ASPC_callable = all(ASPC_callable %in% TRUE),
    PDGF_usable = all(PDGF_usable %in% TRUE),
    ASPC = sum(ASPC, na.rm = TRUE),
    ASPC_raw = sum(ASPC_raw, na.rm = TRUE),
    ASPC_not_callable = sum(ASPC_raw & !ASPC_callable, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== ASPC callability by slide:\n")
print(as.data.frame(aspc_callable_summary), row.names = FALSE)

write.csv(
  aspc_callable_summary,
  file.path(OUT, "aspc_callability_by_slide.csv"),
  row.names = FALSE
)

cat("\n=== Overlapping population flags:\n")
print(table(cells$n_pop_flags))

message("\n03_classify.R complete -> ", file.path(OUT, "cells_classified.rds"))
