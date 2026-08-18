###############################################################################
# 01_channels.R  (v3)
# Build the per-slide channel map (chmap) DIRECTLY from OME-XML metadata.
#
# WHY v3: v2 fell back to jobfile.xml when OME parsing was skipped, producing
#   DAPI_all = 1,4,6,9,11,14,17,20 (acquisition-cycle logic) when the physical
#   plane order is actually OME order: 8 DAPI planes first, then markers.
#   Every downstream measurement built masks from marker planes -> garbage.
#   The jobfile fallback is REMOVED. OME metadata is the only source of truth,
#   and an empirical DAPI check verifies the result against the pixels.
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

OUT <- "results_ASPC_NEPC"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

slides <- readRDS(file.path(OUT, "slides_final.rds"))

## ---- PMOT32 path correction: point at the .ome.tiff, not the plain .tiff ----
PMOT32_TIFF <- paste0("/Volumes/T9/Case Studies/",
                      "20260420_JJ_full panel_test_PMOT32/",
                      "20260420_141534_1_0RPELT_JJ_full panel.ome.tiff")
.k <- match("PMOT32", slides$slide_id)
if (!is.na(.k) && file.exists(PMOT32_TIFF)) {
  if (!identical(slides$tiff[.k], PMOT32_TIFF))
    message("PMOT32 path corrected -> ", basename(PMOT32_TIFF))
  slides$tiff[.k] <- PMOT32_TIFF
}

stopifnot(all(file.exists(slides$tiff)))

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
## BUILD chmap FROM CHANNEL NAMES
###############################################################################
build_chmap <- function(nm) {
  n <- length(nm)
  low <- tolower(nm)

  is_dapi  <- grepl("^dapi", low)                       # "DAPI", "DAPI-JJ CEBPg"
  DAPI_all <- which(is_dapi)

  ## cycle-matched DAPI for the AR/DLL3 round; fall back to the last DAPI
  dapi_ar <- which(is_dapi & grepl("ar", low) & grepl("dll3", low))
  if (!length(dapi_ar)) dapi_ar <- which(is_dapi & grepl("ar&|dll3", low))
  if (!length(dapi_ar)) dapi_ar <- tail(DAPI_all, 1)
  DAPI_AR <- dapi_ar[1]

  ## marker planes: exclude DAPI-prefixed names so "DAPI-JJ CEBPg" never wins
  mk <- function(pat) {
    i <- which(!is_dapi & grepl(pat, low))
    if (!length(i)) NA_integer_ else i[1]
  }
  chmap <- list(
    DAPI_all = DAPI_all,
    DAPI_AR  = DAPI_AR,
    CEBPg    = mk("cebp"),
    FOLH1    = mk("folh1|psma"),
    DLL3     = mk("dll3"),
    AR       = mk("\\bar\\b|androgen"),
    panCK    = mk("pan-?ck|pancytokeratin"),
    PDGF     = mk("pdgf")
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
inv  <- list(); info <- list(); rep_rows <- list()

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
          paste(sprintf("%s=%s", c("CEBPg","FOLH1","DLL3","AR","panCK","PDGF"),
                        unlist(cm[c("CEBPg","FOLH1","DLL3","AR","panCK","PDGF")])),
                collapse = "  "))
  message("  control pairs  : ", cm$n_ctrl_pairs)

  miss <- names(which(vapply(cm[c("CEBPg","FOLH1","DLL3","AR","panCK","PDGF")],
                             function(v) is.na(v[1]), logical(1))))
  if (length(miss))
    warning(sid, ": marker(s) not found in OME names: ",
            paste(miss, collapse = ", "), call. = FALSE)

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

resolution <- do.call(rbind, rep_rows)
write.csv(resolution, file.path(OUT, "channel_resolution.csv"), row.names = FALSE)
print(resolution[, c("slide_id","n_planes","dapi_declared","agrees",
                     "min_dapi_cor","DAPI_AR")], row.names = FALSE)

saveRDS(list(slides = slides, info = info, px_um = PX_UM),
        file.path(OUT, "channels.rds"))

bad <- resolution$slide_id[!resolution$agrees |
                           (!is.na(resolution$min_dapi_cor) &
                            resolution$min_dapi_cor < 0.85)]
if (length(bad))
  message("\n!! Review before segmentation: ", paste(bad, collapse = ", "))

message("\n01_channels.R complete -> ", file.path(OUT, "channels.rds"))
