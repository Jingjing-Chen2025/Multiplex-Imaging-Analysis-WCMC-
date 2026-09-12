###############################################################################
## MARKER USABILITY QC, 6 um -- AUDIT ONLY, NOT USED FOR CELL CALLING
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

    data.frame(
      slide_id = sid,
      marker = m,
      marker_frac_pos = marker_frac,
      blank_frac_pos = blank_frac,
      enrichment_over_blank = enrichment,
      median_blank_slope = slope,
      median_blank_r2 = r2,

      ## Important:
      ## marker_usable is retained only for QC reporting.
      ## It is NOT used to suppress cell calls.
      marker_usable = TRUE,
      fail_reason = "audit_only_not_used_for_calling",

      stringsAsFactors = FALSE
    )
  }))
}))

cat("\n=== Marker usability audit, primary 6 um; NOT used for cell calling:\n")
print(marker_usability, row.names = FALSE)

write.csv(
  marker_usability,
  file.path(OUT, "qc_marker_usability_primary_6um.csv"),
  row.names = FALSE
)

## Add usability columns for compatibility, but keep all markers callable.
for (m in names(MARKER_COMPARTMENT)) {
  usable_col <- paste0(m, "_usable")
  cells[[usable_col]] <- TRUE
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
## POPULATION CALLS -- PIXEL-CORRECTED CELL CALLS, NO SLIDE-LEVEL BLANK QC GATE
###############################################################################

## Raw marker calls are already based on:
##   - pixel-level matched-blank subtraction from 02
##   - AF puncta masking from 02
##   - corrected per-cell sub and frac_pos gates in 03
##
## Therefore do not suppress cells based on slide-level mean blank positivity.

for (m in names(MARKER_COMPARTMENT)) {
  pos_col <- paste0(m, "_pos")
  qc_col <- paste0(m, "_pos_qc")

  ## Compatibility column.
  ## Since slide-level marker usability QC is audit-only now,
  ## the QC-gated call is identical to the corrected marker call.
  cells[[qc_col]] <- cells[[pos_col]] %in% TRUE
}

cells$CEBPg_neg <- !(cells$CEBPg_pos %in% TRUE)
cells$CEBPg_neg_qc <- !(cells$CEBPg_pos_qc %in% TRUE)

## ASPC:
## Requested definition simplifies to PDGF_pos.
## Final ASPC now uses corrected PDGF cell calls directly.
cells$ASPC_raw <- cells$PDGF_pos %in% TRUE
cells$ASPC_callable <- TRUE
cells$ASPC_not_callable <- FALSE
cells$ASPC <- cells$ASPC_raw

cells$ASPC_CEBPg_pos <- cells$ASPC & (cells$CEBPg_pos %in% TRUE)
cells$ASPC_CEBPg_neg <- cells$ASPC & !(cells$CEBPg_pos %in% TRUE)

cells$ASPC_raw_CEBPg_pos <- cells$ASPC_raw & (cells$CEBPg_pos %in% TRUE)
cells$ASPC_raw_CEBPg_neg <- cells$ASPC_raw & !(cells$CEBPg_pos %in% TRUE)

## Compatibility columns for scripts that expect them.
cells$ASPC_CEBPg_known <- cells$ASPC
cells$ASPC_CEBPg_unknown <- FALSE
cells$ASPC_raw_CEBPg_known <- cells$ASPC_raw
cells$ASPC_raw_CEBPg_unknown <- FALSE

## Epithelial cancer.
cells$EPI_CA_raw <- (cells$panCK_pos %in% TRUE) &
                    ((cells$FOLH1_pos %in% TRUE) | (cells$AR_pos %in% TRUE))

cells$EPI_CA_callable <- TRUE
cells$EPI_CA_not_callable <- FALSE
cells$EPI_CA <- cells$EPI_CA_raw

## Neuroendocrine cancer.
cells$NE_CA_raw <- cells$DLL3_pos %in% TRUE
cells$NE_CA_callable <- TRUE
cells$NE_CA_not_callable <- FALSE
cells$NE_CA <- cells$NE_CA_raw

## Blank-derived false-positive population controls.
## These are retained as QC/audit metrics only. They do not suppress marker calls.
cells$ASPC_blank <- cells$PDGF_blank_pos %in% TRUE

cells$EPI_CA_blank <- (cells$panCK_blank_pos %in% TRUE) &
                      ((cells$FOLH1_blank_pos %in% TRUE) |
                       (cells$AR_blank_pos %in% TRUE))

cells$NE_CA_blank <- cells$DLL3_blank_pos %in% TRUE

## Display label.
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
cells$pop[cells$multi_pop] <- "Multi"
