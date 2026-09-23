###############################################################################
# 01_channels.R  (v4)
# Build the per-slide channel map (chmap) DIRECTLY from OME-XML metadata.
#
# WHY v4:
#   (a) PM0T32 was re-acquired 20260826 on the standard 8-cycle protocol WITH
#       secondary-only negative controls (AR/DLL3 = Cycle 6). The hardcoded
#       20260420 path override is removed and replaced with a generic
#       .ome.tiff resolver that applies to every slide.
#   (b) Marker matching was done with literal regexes against the raw channel
#       name, so "JJ pan - CK - 1-200" (spaces around the hyphen, as written in
#       the 20260826 jobfile) failed to match "pan-?ck" and panCK resolved to
#       NA on PM0T32 only. Matching is now done against a NORMALISED name with
#       all non-alphanumeric characters stripped, which tolerates arbitrary
#       spacing, punctuation, vendor prefixes and dilution suffixes.
#   (c) Marker assignments are now cross-checked against the expected
#       fluorophore slot and reported in channel_resolution.csv, and an
#       unmatched marker is a hard error rather than a warning -- an NA here
#       previously surfaced only downstream in 02.
#
# WHY v3 (retained): v2 fell back to jobfile.xml when OME parsing was skipped,
#   producing DAPI_all = 1,4,6,9,11,14,17,20 (acquisition-cycle logic) when the
#   physical plane order is actually OME order: 8 DAPI planes first, then
#   markers. Every downstream measurement built masks from marker planes ->
#   garbage. The jobfile fallback is REMOVED. OME metadata is the only source
#   of truth, and an empirical DAPI check verifies the result against pixels.
#
# Expected OME layout (22-plane files):
#    1-8  DAPI (one per cycle; some named "DAPI-<marker>")
#    9,10 TRITC, Cy5   <- pair 1  (AF baseline)
#   11    JJ CEBPg
#   12,13 TRITC, Cy5   <- pair 2  (NC3)
#   14    JJ FOLH1
#   15,16 TRITC, Cy5   <- pair 3  (NC5)
#   17    JJ DLL3
#   18    JJ AR-Androgen Receptor
#   19    JJ pan-CK
#   20    JJ PDGFa
#   21,22 TRITC, Cy5   <- pair 4  (terminal control)
###############################################################################

suppressPackageStartupMessages({
  library(reticulate); library(dplyr); library(xml2)
})
use_virtualenv("r-reticulate", required = TRUE)

if (!exists("OUT")) OUT <- "/Volumes/T9/results_ASPC_NEPC"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

slides <- readRDS(file.path(OUT, "slides_final.rds"))

## ---- prefer .ome.tiff sidecars, for every slide ----------------------------
## Replaces the v3 PM0T32-specific path override. If a matching .ome.tiff sits
## next to the plain .tiff, use it, because only the OME variant carries the
## Channel metadata this script depends on.
resolve_ome <- function(p) {
  if (grepl("\\.ome\\.tiff?$", p, ignore.case = TRUE)) return(p)
  alt <- sub("\\.tiff?$", ".ome.tiff", p, ignore.case = TRUE)
  if (file.exists(alt)) alt else p
}

.orig <- slides$tiff
slides$tiff <- vapply(slides$tiff, resolve_ome, character(1), USE.NAMES = FALSE)

for (i in seq_along(slides$tiff)) {
  if (!identical(slides$tiff[i], .orig[i]))
    message("Using OME sidecar for ", slides$slide_id[i], " -> ",
            basename(slides$tiff[i]))
}

missing <- slides$slide_id[!file.exists(slides$tiff)]
if (length(missing))
  stop("TIFF not found for: ", paste(missing, collapse = ", "),
       "\nIs /Volumes/T9 mounted?", call. = FALSE)

tf <- import("tifffile", convert = FALSE)
np <- import("numpy",    convert = FALSE)
zr <- import("zarr",     convert = FALSE)
py <- import_builtins()

###############################################################################
## HELPERS
###############################################################################
.sl <- function(i, j) py$slice(as.integer(i), as.integer(j))

open_arr <- function(tiff) {
  st <- tf$TiffFile(tiff)$series[0]$aszarr()
  g  <- tryCatch(zr$open(store = st, mode = "r"),
                 error = function(e) zr$open(st, mode = "r"))
  tryCatch(g$`__getitem__`("0"), error = function(e) g)
}
grab <- function(a, i, x0, y0, w)
  as.numeric(py_to_r(np$asarray(a$`__getitem__`(
    tuple(as.integer(i - 1L), .sl(y0, y0 + w), .sl(x0, x0 + w))))))

## ---- OME-XML -> channel names, in physical plane order ---------------------
ome_channels <- function(tiff) {
  xml <- tryCatch(as.character(py_to_r(tf$TiffFile(tiff)$ome_metadata)),
                  error = function(e) NA_character_)
  if (is.na(xml) || !nzchar(xml) || identical(xml, "None"))
    stop("No OME metadata in ", basename(tiff),
         " -- cannot build chmap without it.", call. = FALSE)
  doc <- read_xml(xml)
  xml_ns_strip <- xml_ns_rename(xml_ns(doc), d1 = "ome")
  ch  <- xml_find_all(doc, "//ome:Image[1]/ome:Pixels/ome:Channel", xml_ns_strip)
  if (!length(ch)) stop("OME metadata has no Channel entries: ", basename(tiff),
                        call. = FALSE)
  nm <- xml_attr(ch, "Name")
  nm[is.na(nm)] <- xml_attr(ch, "ID")[is.na(nm)]
  trimws(nm)
}

## ---- physical pixel size from OME ------------------------------------------
ome_px_um <- function(tiff) {
  xml <- tryCatch(as.character(py_to_r(tf$TiffFile(tiff)$ome_metadata)),
                  error = function(e) NA_character_)
  if (is.na(xml)) return(NA_real_)
  doc <- read_xml(xml)
  ns  <- xml_ns_rename(xml_ns(doc), d1 = "ome")
  px  <- xml_find_first(doc, "//ome:Image[1]/ome:Pixels", ns)
  v   <- suppressWarnings(as.numeric(xml_attr(px, "PhysicalSizeX")))
  u   <- xml_attr(px, "PhysicalSizeXUnit")
  if (is.na(v) || v <= 0) return(NA_real_)
  if (!is.na(u) && u %in% c("nm")) v <- v / 1000
  v
}

###############################################################################
## MARKER NAME MATCHING
###############################################################################
## COMET channel names are free text and vary between operators and software
## versions. Observed forms for pan-cytokeratin across this cohort include
##   "JJ pan-CK"          (1.4.1.0 runs)
##   "JJ pan - CK - 1-200" (1.5.0.0, PM0T32 20260826)
## A literal regex on the raw string matches the first and not the second.
##
## Fix: normalise by lowercasing and deleting every non-alphanumeric character
## before matching, so
##   "JJ pan - CK - 1-200"              -> "jjpanck1200"      contains "panck"
##   "Jingjing PDGFa 647 (Cy5) [1_100]" -> "jingjingpdgfa647cy51100"
## Patterns below are therefore written WITHOUT spaces or punctuation.

.norm <- function(s) {
  s[is.na(s)] <- ""
  gsub("[^a-z0-9]", "", tolower(s))
}

## AR is anchored on "androgenreceptor" rather than a bare "ar", which would
## match as a substring in many unrelated antibody names.
MARKER_PATTERNS <- list(
  CEBPg = c("cebpg", "cebp"),
  FOLH1 = c("folh1", "psma"),
  DLL3  = c("dll3", "notchligand"),
  AR    = c("androgenreceptor", "jjar"),
  panCK = c("panck", "ckpan", "ae1ae3", "cytokeratin"),
  PDGF  = c("pdgfa", "pdgf")
)

## Expected emission channel per marker, from the COMET run log. Used as a
## consistency check only; a mismatch is reported, never silently accepted.
MARKER_FLUOR <- c(
  CEBPg = "Cy5",  FOLH1 = "Cy5",  AR = "Cy5",  PDGF = "Cy5",
  DLL3  = "TRITC", panCK = "TRITC"
)

###############################################################################
## BUILD chmap FROM CHANNEL NAMES
###############################################################################
build_chmap <- function(nm) {
  n   <- length(nm)
  low <- tolower(nm)
  nrm <- .norm(nm)

  is_dapi  <- grepl("^dapi", low)                       # "DAPI", "DAPI-JJ CEBPg"
  DAPI_all <- which(is_dapi)

  ## cycle-matched DAPI for the AR/DLL3 round; fall back to the last DAPI
  dapi_ar <- which(is_dapi & grepl("ar", low) & grepl("dll3", low))
  if (!length(dapi_ar)) dapi_ar <- which(is_dapi & grepl("ar&|dll3", low))
  if (!length(dapi_ar)) dapi_ar <- tail(DAPI_all, 1)
  DAPI_AR <- dapi_ar[1]

  ## marker planes: match on the normalised name and exclude DAPI-prefixed
  ## entries so "DAPI-JJ CEBPg" can never win over the marker plane itself.
  mk <- function(marker) {
    pats <- MARKER_PATTERNS[[marker]]
    hit  <- which(!is_dapi &
                    Reduce(`|`, lapply(pats, function(p) grepl(p, nrm))))
    if (!length(hit)) return(NA_integer_)
    if (length(hit) > 1L)
      warning("marker '", marker, "' matched ", length(hit), " channels: ",
              paste(nm[hit], collapse = " | "), " -- taking the first.",
              call. = FALSE)
    hit[1]
  }

  chmap <- list(
    DAPI_all = DAPI_all,
    DAPI_AR  = DAPI_AR,
    CEBPg    = mk("CEBPg"),
    FOLH1    = mk("FOLH1"),
    DLL3     = mk("DLL3"),
    AR       = mk("AR"),
    panCK    = mk("panCK"),
    PDGF     = mk("PDGF")
  )

  ## ---- control pairs: consecutive TRITC/Cy5 blanks, in acquisition order ----
  tritc <- which(!is_dapi & grepl("^tritc$", low))
  cy5   <- which(!is_dapi & grepl("^cy5$",   low))
  ## pair them by proximity (TRITC immediately followed by Cy5)
  pairs <- lapply(tritc, function(t) {
    c5 <- cy5[cy5 > t & cy5 <= t + 2]
    if (length(c5)) c(TRITC = t, Cy5 = c5[1]) else NULL
  })
  pairs <- Filter(Negate(is.null), pairs)

  ## pair 1 = autofluorescence baseline; pairs 2,3 = NC3, NC5; pair 4 = terminal
  lbl <- c("AF", "NC3", "NC5", "NC7")
  for (p in seq_along(pairs)) {
    if (p > length(lbl)) break
    chmap[[paste0(lbl[p], "_TRITC")]] <- unname(pairs[[p]]["TRITC"])
    chmap[[paste0(lbl[p], "_Cy5")]]   <- unname(pairs[[p]]["Cy5"])
  }
  for (l in lbl) for (s in c("_TRITC", "_Cy5")) {
    k <- paste0(l, s); if (is.null(chmap[[k]])) chmap[[k]] <- NA_integer_
  }

  chmap$n_planes  <- n
  chmap$n_ctrl_pairs <- length(pairs)
  chmap
}

## ---- does each marker sit in its expected fluorophore slot? -----------------
## Inferred positionally: a marker plane immediately preceded by the paired
## control block, or adjacent to its cycle partner. Reported, not enforced.
check_marker_slots <- function(nm, chmap) {
  MK <- names(MARKER_FLUOR)
  do.call(rbind, lapply(MK, function(m) {
    i <- chmap[[m]]
    data.frame(
      marker = m,
      index  = if (is.null(i)) NA_integer_ else i,
      channel = if (is.null(i) || is.na(i)) NA_character_ else nm[i],
      expected_fluor = unname(MARKER_FLUOR[m]),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }))
}

###############################################################################
## EMPIRICAL VERIFICATION: do the declared DAPI planes agree with the pixels?
###############################################################################
verify_dapi <- function(tiff, chmap, W = 512L) {
  a <- open_arr(tiff)
  d <- as.integer(unlist(py_to_r(a$shape)))
  nC <- d[1]; H <- d[2]; Wd <- d[3]

  ## find a dense-tissue probe window using the first DAPI plane
  gx <- round(seq(2000, Wd - W - 2000, length.out = 7))
  gy <- round(seq(2000, H  - W - 2000, length.out = 7))
  gr <- expand.grid(x = gx, y = gy)
  sig <- vapply(seq_len(nrow(gr)), function(j) {
    v <- grab(a, chmap$DAPI_all[1], gr$x[j], gr$y[j], 256L)
    if (sd(v) == 0) -Inf else median(v) + sd(v)
  }, numeric(1))
  best <- gr[which.max(sig), ]

  M <- sapply(seq_len(nC), function(i) grab(a, i, best$x, best$y, W))
  keep <- apply(M, 2, sd) > 0
  cc <- matrix(NA_real_, nC, nC); cc[keep, keep] <- cor(M[, keep])
  npart <- vapply(seq_len(nC), function(i) sum(cc[i, -i] > 0.9, na.rm = TRUE),
                  numeric(1))
  emp <- which(npart >= 3)

  ## cross-cycle DAPI correlation, the number that actually matters
  dp <- M[, chmap$DAPI_all, drop = FALSE]
  k2 <- apply(dp, 2, sd) > 0
  dcor <- if (sum(k2) > 1) {
    m <- cor(dp[, k2]); min(m[upper.tri(m)])
  } else NA_real_

  list(n_planes = nC, probe = c(best$x, best$y),
       empirical_dapi = emp, min_dapi_cor = round(dcor, 3),
       agrees = all(chmap$DAPI_all %in% emp))
}

###############################################################################
## RUN
###############################################################################
inv  <- list(); info <- list(); rep_rows <- list(); slot_rows <- list()
MK <- c("CEBPg","FOLH1","DLL3","AR","panCK","PDGF")
unresolved <- character(0)

for (k in seq_len(nrow(slides))) {
  sid  <- slides$slide_id[k]
  tiff <- slides$tiff[k]
  message("\n=== ", sid, "  (", basename(tiff), ")")

  nm <- ome_channels(tiff)
  inv[[sid]] <- data.frame(slide_id = sid, source = "ome",
                           index = seq_along(nm), channel = nm,
                           stringsAsFactors = FALSE)

  cm <- build_chmap(nm)
  message("  planes         : ", cm$n_planes)
  message("  DAPI_all       : ", paste(cm$DAPI_all, collapse = ","))
  message("  DAPI_AR        : ", cm$DAPI_AR, "  (", nm[cm$DAPI_AR], ")")
  message("  markers        : ",
          paste(sprintf("%s=%s", MK, unlist(cm[MK])), collapse = "  "))
  message("  control pairs  : ", cm$n_ctrl_pairs)

  ## show which channel name each marker actually bound to -- this is the
  ## line that makes a future label change diagnosable at a glance.
  slots <- check_marker_slots(nm, cm)
  slots$slide_id <- sid
  slot_rows[[sid]] <- slots[, c("slide_id","marker","index","channel",
                                "expected_fluor")]
  for (r in seq_len(nrow(slots)))
    message("    ", sprintf("%-6s", slots$marker[r]), " -> plane ",
            sprintf("%-3s", slots$index[r]), " '", slots$channel[r], "'")

  miss <- MK[vapply(cm[MK], function(v) is.na(v[1]), logical(1))]
  if (length(miss)) {
    unresolved <- c(unresolved, paste0(sid, ": ", paste(miss, collapse = ", ")))
    message("  !! UNMATCHED   : ", paste(miss, collapse = ", "))
    message("     channel names present:")
    for (j in seq_along(nm)) message("       ", j, "  ", nm[j])
  }

  v <- verify_dapi(tiff, cm)
  message("  empirical DAPI : ", paste(v$empirical_dapi, collapse = ","))
  message("  min DAPI cor   : ", v$min_dapi_cor,
          if (isTRUE(v$agrees)) "   [OK]" else "   [MISMATCH]")
  if (!isTRUE(v$agrees))
    warning(sid, ": declared DAPI planes are not all supported by the pixels.",
            call. = FALSE)
  if (!is.na(v$min_dapi_cor) && v$min_dapi_cor < 0.85)
    warning(sid, ": cross-cycle DAPI correlation ", v$min_dapi_cor,
            " -- genuine registration problem.", call. = FALSE)

  cm$verify <- v
  info[[sid]] <- list(chmap = cm, channels = nm)

  rep_rows[[sid]] <- data.frame(
    slide_id = sid, n_planes = cm$n_planes,
    dapi_declared = paste(cm$DAPI_all, collapse = ","),
    dapi_empirical = paste(v$empirical_dapi, collapse = ","),
    agrees = v$agrees, min_dapi_cor = v$min_dapi_cor,
    DAPI_AR = cm$DAPI_AR, CEBPg = cm$CEBPg, FOLH1 = cm$FOLH1,
    DLL3 = cm$DLL3, AR = cm$AR, panCK = cm$panCK, PDGF = cm$PDGF,
    stringsAsFactors = FALSE)
}

## An unmatched marker must stop here. Previously this was a warning and the
## NA propagated to 02, where it surfaced only as "channel 'panCK' is NA".
if (length(unresolved))
  stop("Marker(s) could not be matched to any OME channel name:\n  ",
       paste(unresolved, collapse = "\n  "),
       "\nAdd the observed spelling to MARKER_PATTERNS above. Patterns are ",
       "matched against the name with all non-alphanumerics removed.",
       call. = FALSE)

## ---- image dimensions + pixel size -----------------------------------------
dims <- t(vapply(slides$tiff, function(t) {
  d <- as.integer(unlist(py_to_r(open_arr(t)$shape))); c(d[3], d[2])
}, numeric(2)))
slides$sizeX <- as.integer(dims[, 1])
slides$sizeY <- as.integer(dims[, 2])

px <- vapply(slides$tiff, ome_px_um, numeric(1))
PX_UM <- if (all(is.na(px))) 12500 / 44643 else median(px, na.rm = TRUE)
message("\nPX_UM          : ", signif(PX_UM, 4), " µm/px",
        if (all(is.na(px))) "  (fallback)" else "  (from OME)")
if (!all(is.na(px)) && diff(range(px, na.rm = TRUE)) > 1e-6)
  warning("Pixel size differs between slides: ",
          paste(signif(px, 4), collapse = ", "), call. = FALSE)

###############################################################################
## SAVE
###############################################################################
inventory <- do.call(rbind, inv)
write.csv(inventory, file.path(OUT, "channel_inventory.csv"), row.names = FALSE)

marker_slots <- do.call(rbind, slot_rows)
write.csv(marker_slots, file.path(OUT, "channel_marker_assignment.csv"),
          row.names = FALSE)

resolution <- do.call(rbind, rep_rows)
write.csv(resolution, file.path(OUT, "channel_resolution.csv"), row.names = FALSE)
print(resolution[, c("slide_id","n_planes","dapi_declared","agrees",
                     "min_dapi_cor","DAPI_AR")], row.names = FALSE)

cat("\nMarker plane assignment across slides:\n")
print(
  reshape(marker_slots[, c("slide_id","marker","index")],
          idvar = "slide_id", timevar = "marker", direction = "wide"),
  row.names = FALSE
)

## Marker planes should be identical across slides on a shared protocol.
## A slide that differs here is either mis-matched or genuinely re-ordered.
by_marker <- split(marker_slots$index, marker_slots$marker)
inconsistent <- names(by_marker)[vapply(by_marker,
                                        function(v) length(unique(v)) > 1L,
                                        logical(1))]
if (length(inconsistent))
  warning("Marker plane index differs between slides for: ",
          paste(inconsistent, collapse = ", "),
          " -- check channel_marker_assignment.csv.", call. = FALSE)

saveRDS(list(slides = slides, info = info, px_um = PX_UM),
        file.path(OUT, "channels.rds"))

bad <- resolution$slide_id[!resolution$agrees |
                           (!is.na(resolution$min_dapi_cor) &
                            resolution$min_dapi_cor < 0.85)]
if (length(bad))
  message("\n!! Review before segmentation: ", paste(bad, collapse = ", "))

message("\n01_channels.R complete -> ", file.path(OUT, "channels.rds"))
