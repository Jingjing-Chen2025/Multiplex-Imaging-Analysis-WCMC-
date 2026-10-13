###############################################################################
# 02_segment_measure.R  v11.2
#
# Tiled segmentation + corrected per-compartment marker feature extraction.
#
# CHANGE IN v11.2:
#   PMOT32 has been re-run on the standard 8-cycle protocol WITH secondary-only
#   negative controls (AR/DLL3 = Cycle 6). It is therefore no longer a protocol
#   variant and is included by default, exactly like every other slide.
#   The INCLUDE_* exclusion switch has been removed.
#
# v11 image-analysis model:
#   1. Nuclear markers:
#        CEBPD, AR measured in eroded nuclear masks.
#
#   2. Membrane / cytoplasmic markers:
#        PDGFRA, panCK, FOLH1, DLL3 measured in local perinuclear rings.
#        Ring radii: 4, 6, 8 µm.
#        Primary ring radius: 6 µm.
#
#   3. Pixel-level matched-blank subtraction:
#        marker - matched secondary-only negative control.
#
#   4. AF puncta exclusion:
#        pixels bright in both marker and matched blank are masked,
#        dilated, and excluded before cell aggregation.
#
#   5. Blank false-positive control:
#        matched secondary-only negative control is treated as a pseudo-marker,
#        corrected against fluorophore-matched AF channel.
#
#   6. Radius sensitivity:
#        Primary 6 µm features are named *_ring_*.
#        4 and 8 µm sensitivity features are named *_ring_r4um_* and *_ring_r8um_*.
#
# Main output columns used by 03_classify.R:
#
#   Nuclear marker example:
#        CEBPD_nuc_sub
#        CEBPD_nuc_frac_pos
#        CEBPD_nuc_af_frac
#        CEBPD_nuc_blank_sub
#        CEBPD_nuc_blank_frac_pos
#        CEBPD_nuc_blank_af_frac
#
#   Primary 6 µm ring marker example:
#        PDGFRA_ring_sub
#        PDGFRA_ring_frac_pos
#        PDGFRA_ring_af_frac
#        PDGFRA_ring_blank_sub
#        PDGFRA_ring_blank_frac_pos
#        PDGFRA_ring_blank_af_frac
#
#   Sensitivity ring examples:
#        PDGFRA_ring_r4um_sub
#        PDGFRA_ring_r4um_frac_pos
#        PDGFRA_ring_r8um_sub
#        PDGFRA_ring_r8um_frac_pos
###############################################################################

###############################################################################
## USER / RUN DEFAULTS
###############################################################################

if (!exists("OUT")) OUT <- "/Volumes/T9/Case Studies/multiplex_analysis_20261013_PDGFRA_v2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (!exists("ROOT")) ROOT <- "/Volumes/T9/Case Studies"

if (!exists("PROCEED_DESPITE_AF")) PROCEED_DESPITE_AF <- TRUE

## These can be overridden before source("02_segment_measure.R").
if (!exists("TILE_SIZE")) TILE_SIZE <- 2048L
if (!exists("MAX_TILES")) MAX_TILES <- Inf
if (!exists("VERBOSE_EVERY")) VERBOSE_EVERY <- 25L

## v11 classifier does not need old raw marker means.
if (!exists("KEEP_RAW_MEANS")) KEEP_RAW_MEANS <- FALSE

###############################################################################
## SLIDE MANIFEST
###############################################################################
## All six slides now share the 8-cycle protocol with secondary-only negative
## controls, so no slide-specific handling is required here.

MANIFEST <- data.frame(
  slide_id = c("PM0T32", "PM677_A8", "PM12_B3",
               "PM1024_B", "PM423", "PM159T5"),
  group    = c("NEPC", "NEPC", "NEPC",
               "ADENO", "ADENO", "ADENO"),
  tiff = file.path(ROOT, c(
    "20261002_161511_1_INNdFn_JJ_full panel_PMO T32/20261002_161511_1_INNdFn_JJ_full panel_PMO T32.tiff",
    "20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8/20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8.tiff",
    "20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3/20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3.tiff",
    "20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B/20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B.tiff",
    "20260625_131559_2_ZxKUqB_JJ_full panel_PM423/20260625_131559_2_ZxKUqB_JJ_full panel_PM423.tiff",
    "20260625_131600_3_VHgSnt_JJ_full panel_PM159T5/20260625_131600_3_VHgSnt_JJ_full panel_PM159T5.tiff")),
  stringsAsFactors = FALSE
)

resolve_tiff <- function(p) {
  if (grepl("\\.ome\\.tiff?$", p, ignore.case = TRUE)) return(p)
  alt <- sub("\\.tiff?$", ".ome.tiff", p, ignore.case = TRUE)
  if (file.exists(alt)) alt else p
}

MANIFEST$tiff_resolved <- vapply(
  MANIFEST$tiff,
  resolve_tiff,
  character(1),
  USE.NAMES = FALSE
)

cat("\nSlide manifest:\n")
for (i in seq_len(nrow(MANIFEST))) {
  ok <- file.exists(MANIFEST$tiff_resolved[i])
  swp <- MANIFEST$tiff_resolved[i] != MANIFEST$tiff[i]

  cat(sprintf(
    "  %-9s %-6s %s%s  %s\n",
    MANIFEST$slide_id[i],
    MANIFEST$group[i],
    if (ok) "[found]" else "[MISSING]",
    if (swp) " [.ome.tiff]" else "",
    basename(MANIFEST$tiff_resolved[i])
  ))
}

if (any(!file.exists(MANIFEST$tiff_resolved))) {
  stop("Missing TIFF(s) -- is /Volumes/T9 mounted?", call. = FALSE)
}

###############################################################################
## LOAD CHANNEL MAP
###############################################################################

channels_file <- file.path(OUT, "channels.rds")

if (!file.exists(channels_file)) {
  stop("Missing channels.rds. Re-run 01_channels.R first.", call. = FALSE)
}

z <- readRDS(channels_file)

slides <- z$slides
info   <- z$info
PX_UM  <- if (!is.null(z$px_um)) z$px_um else 12500 / 44643

unknown <- setdiff(MANIFEST$slide_id, slides$slide_id)
if (length(unknown)) {
  stop(
    "No channel map for: ",
    paste(unknown, collapse = ", "),
    "\nRe-run 01_channels.R (PMOT32 was re-acquired 20260826 and must be re-mapped).",
    call. = FALSE
  )
}

slides <- slides[match(MANIFEST$slide_id, slides$slide_id), , drop = FALSE]
slides$tiff  <- MANIFEST$tiff_resolved
slides$group <- MANIFEST$group
info <- info[MANIFEST$slide_id]

slides$sizeX <- as.integer(unlist(slides$sizeX))
slides$sizeY <- as.integer(unlist(slides$sizeY))

for (cl in c("tier2", "tier3", "has_negctrl")) {
  if (!cl %in% names(slides)) slides[[cl]] <- TRUE
}

for (cl in c("run", "comet_sw")) {
  if (!cl %in% names(slides)) slides[[cl]] <- NA
}

## All slides must now carry negative controls; the correction model requires it.
if (any(!slides$has_negctrl, na.rm = TRUE)) {
  stop(
    "Slide(s) without secondary-only negative controls: ",
    paste(slides$slide_id[!slides$has_negctrl], collapse = ", "),
    "\nv11 matched-blank subtraction requires NC cycles on every slide.",
    call. = FALSE
  )
}

cat("\nGroups:\n")
print(table(slides$group))

###############################################################################
## PARAMETERS
###############################################################################

try(mem.maxVSize(vsize = 64 * 1024), silent = TRUE)

MIN_AREA_UM2 <- 8
MAX_AREA_UM2 <- 400

MIN_AREA_PX <- MIN_AREA_UM2 / PX_UM^2
MAX_AREA_PX <- MAX_AREA_UM2 / PX_UM^2

QC_WIN <- 512L

## Ring radius sensitivity.
if (!exists("RING_UM_PRIMARY")) RING_UM_PRIMARY <- 6
if (!exists("RING_UM_SET")) RING_UM_SET <- c(4, 6, 8)

.odd <- function(x) {
  x <- as.integer(round(x))
  x + as.integer(x %% 2L == 0L)
}

RING_R_PX_SET <- setNames(
  pmax(3L, .odd(RING_UM_SET / PX_UM)),
  paste0("r", gsub("\\.", "p", as.character(RING_UM_SET)), "um")
)

## Pixel sampling for expensive quantiles/backgrounds.
if (!exists("PIX_SAMPLE_N")) PIX_SAMPLE_N <- 200000L

## AF puncta parameters.
if (!exists("AF_Q")) AF_Q <- 0.995
if (!exists("AF_DILATE_R")) AF_DILATE_R <- 5L
if (!exists("MIN_PIX_MAD")) MIN_PIX_MAD <- 25

PY_ENV <- "r-reticulate"

suppressPackageStartupMessages({
  library(EBImage)
  library(dplyr)
  library(purrr)
  library(tibble)
})

###############################################################################
## MARKERS, CONTROLS, COMPARTMENTS
###############################################################################

MARKERS_CORE <- c("CEBPD", "FOLH1", "DLL3", "AR", "panCK", "PDGFRA")

CTRL_CHANS <- c(
  "AF_Cy5", "AF_TRITC",
  "NC3_Cy5", "NC3_TRITC",
  "NC5_Cy5", "NC5_TRITC",
  "NC7_Cy5", "NC7_TRITC"
)

MARKERS <- c(MARKERS_CORE, CTRL_CHANS)

## Corrected from COMET run log:
##   CEBPD  Cy5
##   FOLH1  Cy5
##   AR     Cy5
##   DLL3   TRITC
##   PDGFRA  Cy5
##   panCK  TRITC
CTRL_OF <- c(
  CEBPD = "NC3_Cy5",
  FOLH1 = "NC3_Cy5",
  AR    = "NC5_Cy5",
  DLL3  = "NC5_TRITC",
  PDGFRA  = "NC5_Cy5",
  panCK = "NC5_TRITC"
)

MARKER_COMPARTMENT <- c(
  CEBPD = "nuc",
  AR    = "nuc",
  PDGFRA  = "ring",
  panCK = "ring",
  FOLH1 = "ring",
  DLL3  = "ring"
)

## For blank false-positive control:
## treat the secondary-only negative control as a pseudo-marker and subtract
## fluorophore-matched AF.
BLANK_REF_OF <- c(
  NC3_Cy5   = "AF_Cy5",
  NC5_Cy5   = "AF_Cy5",
  NC7_Cy5   = "AF_Cy5",
  NC3_TRITC = "AF_TRITC",
  NC5_TRITC = "AF_TRITC",
  NC7_TRITC = "AF_TRITC"
)

###############################################################################
## CHANNEL MAP GUARD
###############################################################################

bad <- character(0)

for (sid in names(info)) {
  cm <- info[[sid]]$chmap
  da <- as.integer(cm$DAPI_all)

  if (!identical(da, seq_along(da))) {
    bad <- c(bad, paste0(sid, ": DAPI_all = ", paste(da, collapse = ",")))
  }

  if (!is.null(cm$verify) && isFALSE(cm$verify$agrees)) {
    bad <- c(bad, paste0(sid, ": failed empirical DAPI check in 01"))
  }

  for (m in MARKERS) {
    v <- cm[[m]]

    if (is.null(v) || !length(v) || is.na(v[1])) {
      bad <- c(bad, paste0(sid, ": channel '", m, "' is NA"))
    } else if (v[1] %in% da) {
      bad <- c(bad, paste0(sid, ": channel '", m, "' collides with DAPI plane"))
    }
  }

  for (m in names(CTRL_OF)) {
    ctl <- CTRL_OF[[m]]
    if (is.null(cm[[ctl]]) || is.na(cm[[ctl]][1])) {
      bad <- c(bad, paste0(sid, ": missing matched control for ", m, " -> ", ctl))
    }
  }

  for (ctl in unique(CTRL_OF)) {
    af <- BLANK_REF_OF[[ctl]]
    if (is.null(af) || is.null(cm[[af]]) || is.na(cm[[af]][1])) {
      bad <- c(bad, paste0(sid, ": missing AF reference for ", ctl, " -> ", af))
    }
  }
}

if (length(bad)) {
  stop(
    "channels.rds is stale, unverified, or incomplete:\n  ",
    paste(unique(bad), collapse = "\n  "),
    "\nRe-run 01_channels.R before segmenting.",
    call. = FALSE
  )
}

message("chmap guard: OK")

###############################################################################
## IMAGE READER
###############################################################################

READ_MODE <- NA_character_

.setup_python <- function(auto_install = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)

  library(reticulate)

  if (nzchar(PY_ENV)) {
    try(use_virtualenv(PY_ENV, required = TRUE), silent = TRUE)
  }

  cfg <- tryCatch(py_config(), error = function(e) NULL)
  if (is.null(cfg)) return(NULL)

  message("Python: ", cfg$python)

  have <- function(m) {
    tryCatch({
      import(m)
      TRUE
    }, error = function(e) FALSE)
  }

  need <- c("tifffile", "numpy")[!vapply(c("tifffile", "numpy"), have, logical(1))]

  if (length(need) && auto_install) {
    message("Installing Python modules: ", paste(need, collapse = ", "))
    try(py_install(c(need, "imagecodecs")), silent = TRUE)
    need <- need[!vapply(need, have, logical(1))]
  }

  if (length(need)) return(NULL)

  if (!have("zarr") && auto_install) {
    message("Installing zarr ...")
    try(py_install("zarr>=3"), silent = TRUE)
  }

  list(
    tf = import("tifffile", convert = FALSE),
    np = import("numpy", convert = FALSE),
    py = import_builtins(),
    zarr = if (have("zarr")) import("zarr", convert = FALSE) else NULL
  )
}

.PY <- .setup_python()

.pynum <- function(x) as.numeric(unlist(reticulate::py_to_r(x)))

if (!is.null(.PY)) {
  .tf <- .PY$tf
  .np <- .PY$np
  .py <- .PY$py
  .zr <- .PY$zarr
  .cache <- new.env(parent = emptyenv())

  if (!is.null(.zr)) {
    get_store <- function(tiff) {
      k <- basename(tiff)

      if (!is.null(.cache[[k]])) return(.cache[[k]])

      st <- .tf$TiffFile(tiff)$series[0]$aszarr()

      g <- tryCatch(
        .zr$open(store = st, mode = "r"),
        error = function(e) .zr$open(st, mode = "r")
      )

      a <- tryCatch(g$`__getitem__`("0"), error = function(e) g)

      .cache[[k]] <- a
      a
    }

    img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
      a <- get_store(tiff)

      sl <- function(i, j) .py$slice(as.integer(i), as.integer(j))

      idx <- if (is.null(w)) {
        reticulate::tuple(
          as.integer(c_index - 1L),
          .py$slice(NULL),
          .py$slice(NULL)
        )
      } else {
        reticulate::tuple(
          as.integer(c_index - 1L),
          sl(y0, y0 + h),
          sl(x0, x0 + w)
        )
      }

      m <- reticulate::py_to_r(.np$asarray(a$`__getitem__`(idx)))
      m <- t(m)
      storage.mode(m) <- "double"
      m
    }

    READ_MODE <- "tifffile+zarr"

  } else {
    CACHE_MAX <- 2L

    get_plane <- function(tiff, c_index) {
      k <- paste0(basename(tiff), "#", c_index)

      if (!is.null(.cache[[k]])) return(.cache[[k]])

      pk <- paste0("pages#", basename(tiff))

      if (is.null(.cache[[pk]])) {
        .cache[[pk]] <- .tf$TiffFile(tiff)$series[0]$pages
      }

      pg <- .cache[[pk]]$`__getitem__`(as.integer(c_index - 1L))
      m <- reticulate::py_to_r(.np$asarray(pg$asarray()))

      keys <- setdiff(ls(.cache), grep("^pages#", ls(.cache), value = TRUE))

      if (length(keys) >= CACHE_MAX) {
        rm(list = keys[1], envir = .cache)
        gc(FALSE)
      }

      assign(k, m, envir = .cache)
      m
    }

    img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
      m <- get_plane(tiff, c_index)

      if (!is.null(w)) {
        m <- m[
          (y0 + 1):min(y0 + h, nrow(m)),
          (x0 + 1):min(x0 + w, ncol(m)),
          drop = FALSE
        ]
      }

      m <- t(m)
      storage.mode(m) <- "double"
      m
    }

    READ_MODE <- "tifffile+pagecache"
  }

} else if (requireNamespace("RBioFormats", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(rJava)
    library(RBioFormats)
  })

  img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
    sub <- list(C = as.integer(c_index))

    if (!is.null(w)) {
      sub <- c(
        sub,
        list(
          x = as.integer(x0),
          y = as.integer(y0),
          w = as.integer(w),
          h = as.integer(h)
        )
      )
    }

    a <- as.array(RBioFormats::read.image(
      tiff,
      series = 1,
      subset = sub,
      normalize = FALSE,
      filter.metadata = TRUE
    ))

    if (length(dim(a)) > 2) a <- a[, , 1]
    storage.mode(a) <- "double"
    a
  }

  READ_MODE <- "bioformats"

} else {
  stop(
    "No usable image reader. Install reticulate/tifffile/zarr or RBioFormats.",
    call. = FALSE
  )
}

message("READ MODE: ", READ_MODE)
message("PX_UM: ", signif(PX_UM, 4), " µm/px")
message("TILE_SIZE: ", TILE_SIZE)
message("Area filter: ", round(MIN_AREA_PX), "-", round(MAX_AREA_PX), " px")
message("Ring radii: ", paste(RING_UM_SET, collapse = ", "), " µm")
message("Primary ring: ", RING_UM_PRIMARY, " µm")
message("Ring px: ", paste(names(RING_R_PX_SET), RING_R_PX_SET, sep = "=", collapse = ", "))

###############################################################################
## STORE GEOMETRY CHECK
###############################################################################

if (grepl("^tifffile\\+zarr", READ_MODE)) {
  dimtab <- do.call(rbind, lapply(seq_len(nrow(slides)), function(k) {
    a <- get_store(slides$tiff[k])
    d <- .pynum(a$shape)
    ch <- .pynum(a$chunks)

    data.frame(
      slide_id = slides$slide_id[k],
      group = slides$group[k],
      planes = d[1],
      sizeX = d[3],
      sizeY = d[2],
      chunk = paste(ch[2], ch[3], sep = "x"),
      aligned = TILE_SIZE %% ch[2] == 0 && TILE_SIZE %% ch[3] == 0,
      n_dapi = length(info[[slides$slide_id[k]]]$chmap$DAPI_all),
      row.names = NULL
    )
  }))

  cat("\nStore geometry:\n")
  print(dimtab, row.names = FALSE)

  bad_dim <- which(dimtab$sizeX != slides$sizeX | dimtab$sizeY != slides$sizeY)

  if (length(bad_dim)) {
    warning(
      "Dimensions differ from channels.rds for: ",
      paste(dimtab$slide_id[bad_dim], collapse = ", "),
      " -- using on-disk values.",
      call. = FALSE
    )

    slides$sizeX <- as.integer(dimtab$sizeX)
    slides$sizeY <- as.integer(dimtab$sizeY)
  }
}

###############################################################################
## HELPERS
###############################################################################

.qsafe <- function(x, p, default = NA_real_, n_sample = PIX_SAMPLE_N) {
  x <- x[is.finite(x)]

  if (!length(x)) return(default)

  if (length(x) > n_sample) {
    set.seed(1L)
    x <- sample(x, n_sample)
  }

  as.numeric(quantile(x, p, na.rm = TRUE, names = FALSE))
}

.weighted_sum_by_label <- function(values, labels, nmax) {
  out <- rep(0, nmax)

  if (!length(values)) return(out)

  rs <- rowsum(values, group = labels, reorder = FALSE)
  ii <- suppressWarnings(as.integer(rownames(rs)))
  keep <- is.finite(ii) & ii > 0L & ii <= nmax

  out[ii[keep]] <- as.numeric(rs[keep, 1])
  out
}

.mean_by_label_valid <- function(img, lab, valid, nmax) {
  l <- as.integer(as.vector(lab))
  v <- as.numeric(as.vector(img))

  ok <- as.vector(valid) &
        is.finite(v) &
        l > 0L &
        l <= nmax

  if (!any(ok)) return(rep(NA_real_, nmax))

  ## Important:
  ## base::tabulate() does not support weights=.
  ## Use rowsum() for weighted numerator and base::tabulate() for denominator.
  num <- .weighted_sum_by_label(v[ok], l[ok], nmax)
  den <- base::tabulate(l[ok], nbins = nmax)

  out <- num / pmax(den, 1)
  out[den == 0] <- NA_real_
  out
}

.frac_by_label_valid <- function(mask, lab, valid, nmax) {
  l <- as.integer(as.vector(lab))
  m <- as.logical(as.vector(mask))

  ok <- as.vector(valid) &
        l > 0L &
        l <= nmax

  if (!any(ok)) return(rep(NA_real_, nmax))

  den <- base::tabulate(l[ok], nbins = nmax)

  ok_pos <- ok & m
  if (any(ok_pos)) {
    num <- base::tabulate(l[ok_pos], nbins = nmax)
  } else {
    num <- rep(0, nmax)
  }

  out <- num / pmax(den, 1)
  out[den == 0] <- NA_real_
  out
}

.make_af_mask <- function(marker_img, blank_img, valid_pix) {
  qm <- .qsafe(marker_img[valid_pix], AF_Q)
  qb <- .qsafe(blank_img[valid_pix],  AF_Q)

  af <- is.finite(marker_img) &
        is.finite(blank_img) &
        valid_pix &
        marker_img >= qm &
        blank_img  >= qb

  af <- EBImage::dilate(
    EBImage::Image(af),
    EBImage::makeBrush(AF_DILATE_R, "disc")
  ) > 0

  as.matrix(af)
}

.correct_marker_pixels <- function(marker_img, blank_img, valid_pix) {
  af <- .make_af_mask(marker_img, blank_img, valid_pix)

  resid <- marker_img - blank_img

  bg_pix <- valid_pix &
            !af &
            is.finite(resid)

  vals <- resid[bg_pix]
  vals <- vals[is.finite(vals)]

  if (length(vals) > PIX_SAMPLE_N) {
    set.seed(1L)
    vals <- sample(vals, PIX_SAMPLE_N)
  }

  bg <- median(vals, na.rm = TRUE)
  sc <- mad(vals, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(bg)) bg <- 0
  if (!is.finite(sc) || sc <= 0) sc <- sd(vals, na.rm = TRUE)
  if (!is.finite(sc) || sc <= 0) sc <- MIN_PIX_MAD

  pix_thr <- bg + max(MIN_PIX_MAD, 3 * sc)

  sub <- pmax(resid - bg, 0)
  sub[af] <- NA_real_

  pos_pix <- resid > pix_thr &
             !af &
             valid_pix

  list(
    sub = sub,
    pos_pix = pos_pix,
    af = af,
    bg = bg,
    pix_thr = pix_thr,
    noise = sc
  )
}

.aggregate_corrected_features <- function(cc, lab, valid_lab, nmax) {
  valid_comp <- valid_lab & !cc$af

  data.frame(
    sub      = .mean_by_label_valid(cc$sub,     lab, valid_comp, nmax),
    frac_pos = .frac_by_label_valid(cc$pos_pix, lab, valid_comp, nmax),
    af_frac  = .frac_by_label_valid(cc$af,      lab, valid_lab,  nmax),
    bg       = rep(cc$bg,      nmax),
    pix_thr  = rep(cc$pix_thr, nmax),
    noise    = rep(cc$noise,   nmax),
    stringsAsFactors = FALSE
  )
}

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

.make_ring_labels <- function(nuc, nuc_mask, tis, d_s) {
  max_px <- max(RING_R_PX_SET)

  max_band <- dilate(nuc_mask, makeBrush(max_px, "disc")) &
              !nuc_mask &
              tis

  local_mask <- nuc_mask | max_band

  ## One propagation at largest radius, then subset to 4 / 6 / 8 µm.
  cell_max <- propagate(d_s, seeds = nuc, mask = local_mask, lambda = 1e-4)

  rings <- lapply(seq_along(RING_UM_SET), function(i) {
    r_um <- RING_UM_SET[i]
    r_px <- RING_R_PX_SET[i]

    band <- dilate(nuc_mask, makeBrush(r_px, "disc")) &
            !nuc_mask &
            tis

    ring <- cell_max
    ring[!band] <- 0
    ring
  })

  names(rings) <- vapply(RING_UM_SET, .ring_tag, character(1))

  rm(cell_max, max_band, local_mask)
  rings
}

###############################################################################
## TILE PROCESSING
###############################################################################

USE_CELLPOSE <- FALSE
CP_DIAM <- 36

process_tile <- function(tiff, cm, x0, y0, w, h, sid) {

  ## Registration QC.
  regc <- NA_real_

  if (length(cm$DAPI_all) > 1) {
    qx <- x0 + (w - QC_WIN) %/% 2
    qy <- y0 + (h - QC_WIN) %/% 2

    ref <- as.vector(img_read(tiff, cm$DAPI_all[1], qx, qy, QC_WIN, QC_WIN))

    cors <- vapply(cm$DAPI_all[-1], function(i) {
      v <- as.vector(img_read(tiff, i, qx, qy, QC_WIN, QC_WIN))
      suppressWarnings(cor(ref, v))
    }, numeric(1))

    regc <- suppressWarnings(min(cors, na.rm = TRUE))

    rm(ref, cors)
    gc(FALSE)
  }

  ## DAPI for segmentation.
  dapi_plane <- if (!is.null(cm$DAPI_AR)) cm$DAPI_AR else cm$DAPI_all[1]
  dapi <- img_read(tiff, dapi_plane, x0, y0, w, h)

  mx <- max(dapi, na.rm = TRUE)

  if (!is.finite(mx) || mx <= 0) {
    rm(dapi)
    gc(FALSE)
    return(NULL)
  }

  d <- Image(dapi / mx)
  rm(dapi)

  d_s <- gblur(d, sigma = 1)

  ## Tissue mask.
  tb <- gblur(d, 15)
  tis <- tb > otsu(tb)
  rm(tb)

  tis <- fillHull(opening(tis, makeBrush(15, "disc")))

  if (sum(tis) < 2000) {
    rm(d, d_s, tis)
    gc(FALSE)
    return(NULL)
  }

  ## Nuclear segmentation.
  if (USE_CELLPOSE) {
    cp  <- reticulate::import("cellpose")
    npy <- reticulate::import("numpy")
    mdl <- cp$models$Cellpose(gpu = FALSE, model_type = "nuclei")

    nuc <- Image(reticulate::py_to_r(
      mdl$eval(
        npy$array(as.matrix(d_s)),
        diameter = CP_DIAM,
        channels = c(0L, 0L)
      )[[1]]
    ))

  } else {
    th <- (d_s > otsu(d_s)) & tis
    th <- opening(th, makeBrush(3, "disc"))
    th <- fillHull(th)

    nuc <- watershed(distmap(th), tolerance = 2, ext = 3)
    rm(th)
  }

  sh0 <- computeFeatures.shape(nuc)

  if (is.null(sh0) || !nrow(sh0)) {
    rm(d, d_s, tis, nuc)
    gc(FALSE)
    return(NULL)
  }

  drop <- which(
    sh0[, "s.area"] < MIN_AREA_PX |
      sh0[, "s.area"] > MAX_AREA_PX
  )

  if (length(drop)) {
    nuc <- rmObjects(nuc, drop, reenumerate = TRUE)
  }

  nmax <- max(nuc)

  if (nmax < 10) {
    rm(d, d_s, tis, nuc)
    gc(FALSE)
    return(NULL)
  }

  ## Compartments:
  ##   - eroded nucleus for CEBPD / AR
  ##   - local perinuclear rings for PDGFRA / panCK / FOLH1 / DLL3
  nuc_mask <- nuc > 0
  nuc_ero <- nuc * as.integer(erode(nuc_mask, makeBrush(3, "disc")))

  rings <- .make_ring_labels(nuc, nuc_mask, tis, d_s)

  rm(tis, d_s, d, nuc_mask)
  gc(FALSE)

  ## Corrected marker and blank-as-marker features.
  corr_list <- list()

  for (m in MARKERS_CORE) {
    cp  <- MARKER_COMPARTMENT[[m]]
    ctl <- CTRL_OF[[m]]

    mk  <- img_read(tiff, cm[[m]][1],   x0, y0, w, h)
    blk <- img_read(tiff, cm[[ctl]][1], x0, y0, w, h)

    af_ref_name <- BLANK_REF_OF[[ctl]]
    af_ref <- img_read(tiff, cm[[af_ref_name]][1], x0, y0, w, h)

    if (cp == "nuc") {
      lab <- nuc_ero
      valid_lab <- nuc_ero > 0

      valid_marker <- valid_lab &
                      is.finite(mk) &
                      is.finite(blk)

      cc_marker <- .correct_marker_pixels(mk, blk, valid_marker)

      feat <- .aggregate_corrected_features(
        cc = cc_marker,
        lab = lab,
        valid_lab = valid_lab,
        nmax = nmax
      )

      colnames(feat) <- paste0(.feature_prefix(m, cp), "_", colnames(feat))
      corr_list[[length(corr_list) + 1]] <- feat

      valid_blank <- valid_lab &
                     is.finite(blk) &
                     is.finite(af_ref)

      cc_blank <- .correct_marker_pixels(blk, af_ref, valid_blank)

      blank_feat <- .aggregate_corrected_features(
        cc = cc_blank,
        lab = lab,
        valid_lab = valid_lab,
        nmax = nmax
      )

      colnames(blank_feat) <- paste0(.feature_prefix(m, cp), "_blank_", colnames(blank_feat))
      corr_list[[length(corr_list) + 1]] <- blank_feat

      rm(cc_marker, cc_blank, feat, blank_feat)

    } else {
      union_ring <- Reduce(`|`, lapply(rings, function(x) x > 0))

      valid_marker <- union_ring &
                      is.finite(mk) &
                      is.finite(blk)

      cc_marker <- .correct_marker_pixels(mk, blk, valid_marker)

      valid_blank <- union_ring &
                     is.finite(blk) &
                     is.finite(af_ref)

      cc_blank <- .correct_marker_pixels(blk, af_ref, valid_blank)

      for (r_um in RING_UM_SET) {
        tag <- .ring_tag(r_um)
        lab <- rings[[tag]]
        valid_lab <- lab > 0

        feat <- .aggregate_corrected_features(
          cc = cc_marker,
          lab = lab,
          valid_lab = valid_lab,
          nmax = nmax
        )

        colnames(feat) <- paste0(.feature_prefix(m, cp, r_um), "_", colnames(feat))
        corr_list[[length(corr_list) + 1]] <- feat

        blank_feat <- .aggregate_corrected_features(
          cc = cc_blank,
          lab = lab,
          valid_lab = valid_lab,
          nmax = nmax
        )

        colnames(blank_feat) <- paste0(.feature_prefix(m, cp, r_um), "_blank_", colnames(blank_feat))
        corr_list[[length(corr_list) + 1]] <- blank_feat

        rm(feat, blank_feat)
      }

      rm(union_ring, cc_marker, cc_blank)
    }

    rm(mk, blk, af_ref)
    gc(FALSE)
  }

  corr_df <- do.call(cbind, corr_list)
  rm(corr_list)
  gc(FALSE)

  ## DAPI mean.
  dapi2 <- img_read(tiff, dapi_plane, x0, y0, w, h)
  dv <- .mean_by_label_valid(dapi2, nuc, nuc > 0, nmax)
  rm(dapi2)

  mo <- computeFeatures.moment(nuc)
  sh <- computeFeatures.shape(nuc)

  out <- data.frame(
    slide_id = sid,
    x_px = mo[, "m.cx"] + x0,
    y_px = mo[, "m.cy"] + y0,
    area = sh[, "s.area"],
    perim = sh[, "s.perimeter"],
    dapi_nuc = dv,
    corr_df,
    reg_cor = regc,
    tile_x = x0,
    tile_y = y0,
    stringsAsFactors = FALSE
  )

  rm(nuc, nuc_ero, rings, mo, sh, corr_df)
  gc(FALSE)

  dplyr::filter(out, !is.na(x_px))
}

###############################################################################
## RUN
###############################################################################

if (!exists("RUN_SLIDES")) RUN_SLIDES <- slides$slide_id

stopifnot(all(RUN_SLIDES %in% slides$slide_id))

if (!is.finite(MAX_TILES) && READ_MODE == "tifffile+pagecache") {
  stop(
    "Refusing full run in page-decode mode. Install zarr, restart R.",
    call. = FALSE
  )
}

message(
  "\nRunning slides: ",
  paste(RUN_SLIDES, collapse = ", "),
  " | max tiles = ",
  MAX_TILES
)

.tile_times <- c()

cells_raw <- map_dfr(which(slides$slide_id %in% RUN_SLIDES), function(k) {
  sid <- slides$slide_id[k]
  ## rerun patch: reuse finished per-slide checkpoints when REUSE_CHECKPOINTS
  .ck <- file.path(OUT, paste0("cells_raw_", sid, ".rds"))
  if (isTRUE(get0("REUSE_CHECKPOINTS")) && file.exists(.ck)) {
    message("Reusing checkpoint ", .ck)
    return(readRDS(.ck))
  }
  cm  <- info[[sid]]$chmap

  xs <- seq(0, slides$sizeX[k] - TILE_SIZE, by = TILE_SIZE)
  ys <- seq(0, slides$sizeY[k] - TILE_SIZE, by = TILE_SIZE)

  g <- expand.grid(x = xs, y = ys)
  n_total <- nrow(g)

  if (is.finite(MAX_TILES)) {
    mid <- order(
      (g$x - slides$sizeX[k] / 2)^2 +
        (g$y - slides$sizeY[k] / 2)^2
    )

    g <- g[head(mid, MAX_TILES), ]
  }

  message(
    "=== ",
    sid,
    " [",
    slides$group[k],
    "] : ",
    nrow(g),
    " of ",
    n_total,
    " tiles"
  )

  t0 <- Sys.time()

  out <- map_dfr(seq_len(nrow(g)), function(i) {
    ti <- Sys.time()

    if (i == 1 || i %% VERBOSE_EVERY == 0) {
      message(
        "   starting tile ",
        i,
        "/",
        nrow(g),
        " at x=",
        g$x[i],
        " y=",
        g$y[i]
      )
    }

    r <- tryCatch(
      process_tile(
        slides$tiff[k],
        cm,
        g$x[i],
        g$y[i],
        TILE_SIZE,
        TILE_SIZE,
        sid
      ),
      error = function(e) {
        message("   tile ", i, " FAILED: ", conditionMessage(e))
        NULL
      }
    )

    dt <- as.numeric(difftime(Sys.time(), ti, units = "secs"))
    .tile_times <<- c(.tile_times, dt)

    if (i == 1 || i %% VERBOSE_EVERY == 0) {
      message(
        "   finished tile ",
        i,
        "/",
        nrow(g),
        " [",
        round(dt, 1),
        " s]",
        if (is.null(r)) " -> no cells" else paste0(" -> ", nrow(r), " cells")
      )
    }

    gc(FALSE)
    r
  })

  message(
    "   -> ",
    nrow(out),
    " cells in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1),
    " min"
  )

  saveRDS(out, file.path(OUT, paste0("cells_raw_", sid, ".rds")))

  out
})

if (!nrow(cells_raw)) {
  stop("No cells segmented -- check tissue mask / thresholds.", call. = FALSE)
}

cells_raw <- left_join(
  cells_raw,
  slides[, c(
    "slide_id",
    "group",
    "run",
    "comet_sw",
    "has_negctrl",
    "tier2",
    "tier3"
  )],
  by = "slide_id"
)

out_file <- if (is.finite(MAX_TILES)) "cells_raw_SMOKETEST.rds" else "cells_raw.rds"
saveRDS(cells_raw, file.path(OUT, out_file))

###############################################################################
## QC OUTPUTS
###############################################################################

qc <- cells_raw |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    min_dapi_cor = suppressWarnings(min(reg_cor, na.rm = TRUE)),
    median_area_px = median(area, na.rm = TRUE),
    median_area_um2 = median(area, na.rm = TRUE) * PX_UM^2,
    .groups = "drop"
  )

cat("\nSegmentation QC:\n")
print(as.data.frame(qc))

write.csv(qc, file.path(OUT, "qc_registration.csv"), row.names = FALSE)

## Corrected primary-feature audit.
corr_audit <- do.call(rbind, lapply(unique(cells_raw$slide_id), function(sid) {
  d <- cells_raw[cells_raw$slide_id == sid, ]

  do.call(rbind, lapply(MARKERS_CORE, function(m) {
    cp <- MARKER_COMPARTMENT[[m]]
    pre <- .feature_prefix(m, cp)

    sub_col  <- paste0(pre, "_sub")
    frac_col <- paste0(pre, "_frac_pos")
    af_col   <- paste0(pre, "_af_frac")

    blank_sub_col  <- paste0(pre, "_blank_sub")
    blank_frac_col <- paste0(pre, "_blank_frac_pos")
    blank_af_col   <- paste0(pre, "_blank_af_frac")

    data.frame(
      slide_id = sid,
      marker = m,
      compartment = cp,
      radius_um = ifelse(cp == "ring", RING_UM_PRIMARY, NA),

      median_sub = round(median(d[[sub_col]], na.rm = TRUE), 3),
      p95_sub = round(.qsafe(d[[sub_col]], .95), 3),
      p99_sub = round(.qsafe(d[[sub_col]], .99), 3),
      median_frac_pos = round(median(d[[frac_col]], na.rm = TRUE), 4),
      p95_frac_pos = round(.qsafe(d[[frac_col]], .95), 4),
      p99_frac_pos = round(.qsafe(d[[frac_col]], .99), 4),
      median_af_frac = round(median(d[[af_col]], na.rm = TRUE), 4),
      p99_af_frac = round(.qsafe(d[[af_col]], .99), 4),

      blank_median_sub = round(median(d[[blank_sub_col]], na.rm = TRUE), 3),
      blank_p95_sub = round(.qsafe(d[[blank_sub_col]], .95), 3),
      blank_p99_sub = round(.qsafe(d[[blank_sub_col]], .99), 3),
      blank_median_frac_pos = round(median(d[[blank_frac_col]], na.rm = TRUE), 4),
      blank_p95_frac_pos = round(.qsafe(d[[blank_frac_col]], .95), 4),
      blank_p99_frac_pos = round(.qsafe(d[[blank_frac_col]], .99), 4),
      blank_median_af_frac = round(median(d[[blank_af_col]], na.rm = TRUE), 4),
      blank_p99_af_frac = round(.qsafe(d[[blank_af_col]], .99), 4),

      stringsAsFactors = FALSE
    )
  }))
}))

cat("\nCorrected primary marker-feature audit:\n")
print(corr_audit, row.names = FALSE)

write.csv(
  corr_audit,
  file.path(OUT, "qc_corrected_marker_features_primary_6um.csv"),
  row.names = FALSE
)

## Radius sensitivity feature audit before classification.
radius_feature_audit <- do.call(rbind, lapply(unique(cells_raw$slide_id), function(sid) {
  d <- cells_raw[cells_raw$slide_id == sid, ]

  ring_markers <- names(MARKER_COMPARTMENT)[MARKER_COMPARTMENT == "ring"]

  do.call(rbind, lapply(ring_markers, function(m) {
    do.call(rbind, lapply(RING_UM_SET, function(r_um) {
      pre <- .feature_prefix(m, "ring", r_um)

      sub_col  <- paste0(pre, "_sub")
      frac_col <- paste0(pre, "_frac_pos")
      af_col   <- paste0(pre, "_af_frac")

      blank_sub_col  <- paste0(pre, "_blank_sub")
      blank_frac_col <- paste0(pre, "_blank_frac_pos")
      blank_af_col   <- paste0(pre, "_blank_af_frac")

      data.frame(
        slide_id = sid,
        marker = m,
        radius_um = r_um,

        median_sub = round(median(d[[sub_col]], na.rm = TRUE), 3),
        p95_sub = round(.qsafe(d[[sub_col]], .95), 3),
        p99_sub = round(.qsafe(d[[sub_col]], .99), 3),
        median_frac_pos = round(median(d[[frac_col]], na.rm = TRUE), 4),
        p95_frac_pos = round(.qsafe(d[[frac_col]], .95), 4),
        p99_frac_pos = round(.qsafe(d[[frac_col]], .99), 4),
        median_af_frac = round(median(d[[af_col]], na.rm = TRUE), 4),
        p99_af_frac = round(.qsafe(d[[af_col]], .99), 4),

        blank_median_sub = round(median(d[[blank_sub_col]], na.rm = TRUE), 3),
        blank_p95_sub = round(.qsafe(d[[blank_sub_col]], .95), 3),
        blank_p99_sub = round(.qsafe(d[[blank_sub_col]], .99), 3),
        blank_median_frac_pos = round(median(d[[blank_frac_col]], na.rm = TRUE), 4),
        blank_p95_frac_pos = round(.qsafe(d[[blank_frac_col]], .95), 4),
        blank_p99_frac_pos = round(.qsafe(d[[blank_frac_col]], .99), 4),
        blank_median_af_frac = round(median(d[[blank_af_col]], na.rm = TRUE), 4),
        blank_p99_af_frac = round(.qsafe(d[[blank_af_col]], .99), 4),

        stringsAsFactors = FALSE
      )
    }))
  }))
}))

cat("\nRadius feature audit, ring markers:\n")
print(radius_feature_audit, row.names = FALSE)

write.csv(
  radius_feature_audit,
  file.path(OUT, "qc_radius_feature_summary_4_6_8um.csv"),
  row.names = FALSE
)

if (any(qc$min_dapi_cor < 0.85, na.rm = TRUE)) {
  warning(
    "Low cross-cycle DAPI correlation on: ",
    paste(qc$slide_id[qc$min_dapi_cor < 0.85], collapse = ", "),
    call. = FALSE
  )
}

if (any(qc$median_area_um2 < 15 | qc$median_area_um2 > 150, na.rm = TRUE)) {
  warning(
    "Median nuclear area outside 15-150 um^2 on: ",
    paste(
      qc$slide_id[
        qc$median_area_um2 < 15 |
          qc$median_area_um2 > 150
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(.tile_times)) {
  message("\nMedian tile time: ", round(median(.tile_times), 1), " s")
  message("Total tiles: ", length(.tile_times))
  message("Total tile time: ", round(sum(.tile_times) / 3600, 2), " h")
}

message("\n02_segment_measure.R complete -> ", file.path(OUT, out_file))
