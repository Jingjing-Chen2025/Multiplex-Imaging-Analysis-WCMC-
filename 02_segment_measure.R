###############################################################################
# 02_segment_measure.R  v11.3
#
# Tiled segmentation + corrected per-compartment marker feature extraction.
#
# CHANGE IN v11.3:
#   Registration QC now uses robust informative-window cross-cycle DAPI
#   correlation, treating empty or flat windows as unknown rather than failed.
#
# CHANGE IN v11.2:
#   PM0T32 has been re-run on the standard 8-cycle protocol WITH secondary-only
#   negative controls (AR/DLL3 = Cycle 6). It is therefore no longer a protocol
#   variant and is included by default, exactly like every other slide.
#   The INCLUDE_* exclusion switch has been removed.
#
# v11 image-analysis model:
#   1. Nuclear markers:
#        CEBPg, AR measured in eroded nuclear masks.
#
#   2. Membrane / cytoplasmic markers:
#        PDGF, panCK, FOLH1, DLL3 measured in local perinuclear rings.
#        Ring radii: 4, 6, 8 µm.
#        Primary ring radius: 6 µm.
#
#   3. Pixel-level scaled matched-blank subtraction:
#        marker - (blank intercept + blank slope * matched secondary-only
#        negative control), fit on valid non-AF pixels after upper-tail trimming.
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
#        CEBPg_nuc_sub
#        CEBPg_nuc_frac_pos
#        CEBPg_nuc_af_frac
#        CEBPg_nuc_blank_sub
#        CEBPg_nuc_blank_frac_pos
#        CEBPg_nuc_blank_af_frac
#
#   Primary 6 µm ring marker example:
#        PDGF_ring_sub
#        PDGF_ring_frac_pos
#        PDGF_ring_af_frac
#        PDGF_ring_blank_sub
#        PDGF_ring_blank_frac_pos
#        PDGF_ring_blank_af_frac
#
#   Sensitivity ring examples:
#        PDGF_ring_r4um_sub
#        PDGF_ring_r4um_frac_pos
#        PDGF_ring_r8um_sub
#        PDGF_ring_r8um_frac_pos
###############################################################################

###############################################################################
## USER / RUN DEFAULTS
###############################################################################

if (!exists("OUT")) OUT <- "/Volumes/T9/results_ASPC_NEPC"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

if (!exists("PROCEED_DESPITE_AF")) PROCEED_DESPITE_AF <- TRUE

## These can be overridden before source("02_segment_measure.R").
## TILE_SIZE is the CORE tile edge. Each tile is read with a HALO_PX margin so
## nuclei and rings crossing a tile edge are segmented whole; only cells whose
## centroid falls inside the core are kept, so none are cut or double-counted.
if (!exists("TILE_SIZE")) TILE_SIZE <- 2048L
if (!exists("HALO_PX")) HALO_PX <- 96L
## Re-use per-slide checkpoints (cells_raw_<sid>.rds) unless forced.
if (!exists("FORCE_RERUN")) FORCE_RERUN <- FALSE
## FALSE = only segment RUN_SLIDES and write their checkpoints (used to run
## slides in parallel processes); TRUE = also combine all slides + write QC.
if (!exists("COMBINE")) COMBINE <- TRUE
if (!exists("MAX_TILES")) MAX_TILES <- Inf
if (!exists("VERBOSE_EVERY")) VERBOSE_EVERY <- 25L

## v11 classifier does not need old raw marker means.
if (!exists("KEEP_RAW_MEANS")) KEEP_RAW_MEANS <- FALSE

###############################################################################
## SLIDE MANIFEST
###############################################################################
## All six slides now share the 8-cycle protocol with secondary-only negative
## controls, so no slide-specific handling is required here.

## The slide table is owned by 00_config.R (slides_final.rds); it is not
## re-declared here so slide IDs, groups and paths cannot drift between scripts.
.sf <- readRDS(file.path(OUT, "slides_final.rds"))
MANIFEST <- data.frame(
  slide_id = .sf$slide_id,
  group    = as.character(.sf$group),
  tiff     = .sf$tiff,
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
    "\nRe-run 01_channels.R (PM0T32 was re-acquired 20260826 and must be re-mapped).",
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
if (!exists("DAPI_MIN_SD")) DAPI_MIN_SD <- 10
if (!exists("DAPI_MIN_DYNAMIC_RANGE")) DAPI_MIN_DYNAMIC_RANGE <- 30
if (!exists("REG_LOW_COR_THRESHOLD")) REG_LOW_COR_THRESHOLD <- 0.85

## Ring radius sensitivity.
if (!exists("RING_UM_PRIMARY")) RING_UM_PRIMARY <- 6
if (!exists("RING_UM_SET")) RING_UM_SET <- c(4, 6, 8)

.odd <- function(x) {
  x <- as.integer(round(x))
  x + as.integer(x %% 2L == 0L)
}

## Ring RADII in pixels. Rings are built from a Euclidean distance map, so
## these are true radii (v11.3 passed radius/px to makeBrush(), which takes a
## DIAMETER, so the nominal 6 um ring was only ~3 um wide).
RING_R_PX_SET <- setNames(
  pmax(2, RING_UM_SET / PX_UM),
  paste0("r", gsub("\\.", "p", as.character(RING_UM_SET)), "um")
)

## Pixel sampling for expensive quantiles/backgrounds.
if (!exists("PIX_SAMPLE_N")) PIX_SAMPLE_N <- 200000L

## AF puncta parameters.
if (!exists("AF_Q")) AF_Q <- 0.995
if (!exists("AF_DILATE_R")) AF_DILATE_R <- 5L   # radius in px
if (!exists("MIN_PIX_MAD")) MIN_PIX_MAD <- 25

PY_ENV <- "r-reticulate"

suppressPackageStartupMessages({
  library(EBImage)
  library(MASS)
  library(dplyr)
  library(purrr)
  library(tibble)
})

###############################################################################
## MARKERS, CONTROLS, COMPARTMENTS
###############################################################################

MARKERS_CORE <- c("CEBPg", "FOLH1", "DLL3", "AR", "panCK", "PDGF")

CTRL_CHANS <- c(
  "AF_Cy5", "AF_TRITC",
  "NC3_Cy5", "NC3_TRITC",
  "NC5_Cy5", "NC5_TRITC",
  "NC7_Cy5", "NC7_TRITC"
)

MARKERS <- c(MARKERS_CORE, CTRL_CHANS)

## Matched secondary-only controls, by fluorophore AND secondary species
## (COMET run logs, all 6 slides incl. PM0T32):
##   Cycle 3 NC = anti-Rabbit 647 (Cy5) + anti-Mouse 555 (TRITC)
##   Cycle 5 NC = anti-Mouse 647 (Cy5)  + anti-Rabbit 555 (TRITC)
##   CEBPg  rabbit Cy5    -> NC3_Cy5
##   FOLH1  mouse  Cy5    -> NC5_Cy5
##   AR     mouse  Cy5    -> NC5_Cy5    (cycle 6 uses the cycle-5 mix)
##   DLL3   rabbit TRITC  -> NC5_TRITC  (cycle 6 uses the cycle-5 mix)
##   PDGFa  rabbit Cy5    -> NC3_Cy5    (cycle 7 uses the cycle-3 mix)
##   panCK  mouse  TRITC  -> NC3_TRITC  (cycle 7 uses the cycle-3 mix)
## v11.3 paired FOLH1 with NC3 and PDGF/panCK with NC5, i.e. the
## wrong-species secondary.
CTRL_OF <- c(
  CEBPg = "NC3_Cy5",
  FOLH1 = "NC5_Cy5",
  AR    = "NC5_Cy5",
  DLL3  = "NC5_TRITC",
  PDGF  = "NC3_Cy5",
  panCK = "NC3_TRITC"
)

MARKER_COMPARTMENT <- c(
  CEBPg = "nuc",
  AR    = "nuc",
  PDGF  = "ring",
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

.is_informative_dapi <- function(v,
                                 min_sd = DAPI_MIN_SD,
                                 min_range = DAPI_MIN_DYNAMIC_RANGE) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v) < 1000L) return(FALSE)
  s <- stats::sd(v, na.rm = TRUE)
  r <- as.numeric(stats::quantile(v, 0.99, na.rm = TRUE, names = FALSE) -
                    stats::median(v, na.rm = TRUE))
  is.finite(s) && is.finite(r) && s >= min_sd && r >= min_range
}

.safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 1000L) return(NA_real_)
  if (stats::sd(a[ok], na.rm = TRUE) <= 0) return(NA_real_)
  if (stats::sd(b[ok], na.rm = TRUE) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(a[ok], b[ok]))
}

.dapi_registration_qc <- function(tiff, cm, x0, y0, w, h) {
  out <- list(
    reg_cor = NA_real_,
    reg_cor_min = NA_real_,
    reg_cor_median = NA_real_,
    reg_cor_p05 = NA_real_,
    reg_cor_n = 0L,
    reg_cor_n_low = NA_integer_,
    reg_cor_low_frac = NA_real_,
    reg_qc_informative = FALSE,
    reg_ref_sd = NA_real_,
    reg_ref_range = NA_real_
  )

  if (length(cm$DAPI_all) <= 1L) return(out)

  qx <- x0 + max(0L, (w - QC_WIN) %/% 2L)
  qy <- y0 + max(0L, (h - QC_WIN) %/% 2L)

  ref <- as.vector(img_read(tiff, cm$DAPI_all[1], qx, qy, QC_WIN, QC_WIN))
  ref_ok <- ref[is.finite(ref)]

  if (length(ref_ok) >= 1000L) {
    out$reg_ref_sd <- stats::sd(ref_ok, na.rm = TRUE)
    out$reg_ref_range <- as.numeric(stats::quantile(ref_ok, 0.99, na.rm = TRUE, names = FALSE) -
                                      stats::median(ref_ok, na.rm = TRUE))
  }

  if (!.is_informative_dapi(ref)) return(out)

  cors <- vapply(cm$DAPI_all[-1], function(i) {
    v <- as.vector(img_read(tiff, i, qx, qy, QC_WIN, QC_WIN))
    if (!.is_informative_dapi(v)) return(NA_real_)
    .safe_cor(ref, v)
  }, numeric(1))

  cors <- cors[is.finite(cors)]
  out$reg_cor_n <- length(cors)
  if (!length(cors)) return(out)

  out$reg_qc_informative <- TRUE
  out$reg_cor_min <- suppressWarnings(min(cors, na.rm = TRUE))
  out$reg_cor_median <- suppressWarnings(stats::median(cors, na.rm = TRUE))
  out$reg_cor_p05 <- suppressWarnings(as.numeric(stats::quantile(cors, 0.05, na.rm = TRUE, names = FALSE)))
  out$reg_cor_n_low <- sum(cors < REG_LOW_COR_THRESHOLD, na.rm = TRUE)
  out$reg_cor_low_frac <- mean(cors < REG_LOW_COR_THRESHOLD, na.rm = TRUE)
  out$reg_cor <- out$reg_cor_median
  out
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
    EBImage::makeBrush(2L * as.integer(AF_DILATE_R) + 1L, "disc")
  ) > 0

  as.matrix(af)
}

.fit_blank_affine <- function(marker_img, blank_img, valid_pix, af) {
  out <- list(intercept = 0, slope = 1, r2 = NA_real_)

  ok <- valid_pix &
        !af &
        is.finite(marker_img) &
        is.finite(blank_img)

  if (!any(ok)) return(out)

  marker <- as.numeric(marker_img[ok])
  blank  <- as.numeric(blank_img[ok])

  q_marker <- .qsafe(marker, 0.99, default = Inf)
  q_blank  <- .qsafe(blank,  0.99, default = Inf)

  keep <- is.finite(marker) &
          is.finite(blank) &
          marker <= q_marker &
          blank <= q_blank

  marker <- marker[keep]
  blank  <- blank[keep]

  if (length(marker) > PIX_SAMPLE_N) {
    set.seed(1L)
    idx <- sample.int(length(marker), PIX_SAMPLE_N)
    marker <- marker[idx]
    blank  <- blank[idx]
  }

  if (length(marker) < 3L || stats::sd(blank, na.rm = TRUE) <= 0) {
    return(out)
  }

  fit_df <- data.frame(marker = marker, blank = blank)
  fit <- try(suppressWarnings(MASS::rlm(marker ~ blank, data = fit_df)), silent = TRUE)

  if (inherits(fit, "try-error")) {
    fit <- try(stats::lm(marker ~ blank, data = fit_df), silent = TRUE)
  }

  if (inherits(fit, "try-error")) return(out)

  cf <- try(stats::coef(fit), silent = TRUE)
  if (inherits(cf, "try-error") || length(cf) < 2L || any(!is.finite(cf[1:2]))) {
    return(out)
  }

  fitted <- try(as.numeric(stats::predict(fit, newdata = fit_df)), silent = TRUE)
  if (!inherits(fitted, "try-error") && length(fitted) == nrow(fit_df) && all(is.finite(fitted))) {
    sst <- sum((fit_df$marker - mean(fit_df$marker))^2)
    sse <- sum((fit_df$marker - fitted)^2)
    if (is.finite(sst) && sst > 0 && is.finite(sse)) {
      out$r2 <- max(min(1 - sse / sst, 1), 0)
    }
  }

  out$intercept <- as.numeric(cf[1])
  out$slope <- as.numeric(cf[2])
  out
}

.correct_marker_pixels <- function(marker_img, blank_img, valid_pix) {
  af <- .make_af_mask(marker_img, blank_img, valid_pix)
  blank_fit <- .fit_blank_affine(marker_img, blank_img, valid_pix, af)

  resid <- marker_img - (blank_fit$intercept + blank_fit$slope * blank_img)

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
    noise = sc,
    blank_intercept = blank_fit$intercept,
    blank_slope = blank_fit$slope,
    blank_r2 = blank_fit$r2
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
    blank_intercept = rep(cc$blank_intercept, nmax),
    blank_slope     = rep(cc$blank_slope,     nmax),
    blank_r2        = rep(cc$blank_r2,        nmax),
    stringsAsFactors = FALSE
  )
}

.col_or_na <- function(df, name) {
  if (name %in% names(df)) df[[name]] else rep(NA_real_, nrow(df))
}

.median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(stats::median(x, na.rm = TRUE))
}

.min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}

.quantile_or_na <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE)))
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
  ## Distance (px) of every non-nuclear pixel to the nearest nucleus.
  dm <- distmap(!nuc_mask)

  max_band <- dm > 0 & dm <= max(RING_R_PX_SET) & tis
  local_mask <- nuc_mask | max_band

  ## One propagation at largest radius, then subset to 4 / 6 / 8 µm.
  cell_max <- propagate(d_s, seeds = nuc, mask = local_mask, lambda = 1e-4)

  rings <- lapply(seq_along(RING_UM_SET), function(i) {
    band <- dm > 0 & dm <= RING_R_PX_SET[i] & tis
    ring <- cell_max
    ring[!band] <- 0
    ring
  })

  names(rings) <- vapply(RING_UM_SET, .ring_tag, character(1))

  rm(cell_max, max_band, local_mask, dm)
  rings
}

###############################################################################
## SLIDE-LEVEL TISSUE MASK (low-resolution pyramid level)
###############################################################################
## v11.3 thresholded tissue per tile with Otsu, which splits dense epithelium
## from sparse stroma whenever a tile is entirely tissue -- dropping exactly the
## stromal compartment where ASPC live. The mask is now built once per slide
## on a ~2.2 um/px pyramid level, anchored on the glass/camera background:
## tissue = blurred DAPI > TISSUE_GLASS_MULT x glass level, then closed,
## hole-filled and cleaned of specks. Tiles with no tissue are skipped without
## reading any full-resolution data.

if (!exists("TISSUE_LEVEL_MAX_PX")) TISSUE_LEVEL_MAX_PX <- 6000L
if (!exists("TISSUE_GLASS_MULT")) TISSUE_GLASS_MULT <- 2
if (!exists("TISSUE_MIN_OBJ_MM2")) TISSUE_MIN_OBJ_MM2 <- 0.02

.open_level <- function(tiff, level) {
  st <- .tf$TiffFile(tiff)$series[0]$aszarr()
  g <- tryCatch(
    .zr$open(store = st, mode = "r"),
    error = function(e) .zr$open(st, mode = "r")
  )
  g$`__getitem__`(as.character(level))
}

.pick_mask_level <- function(tiff) {
  lv <- 0L
  repeat {
    a <- tryCatch(.open_level(tiff, lv), error = function(e) NULL)
    if (is.null(a)) return(list(level = lv - 1L, arr = prev))
    d <- .pynum(a$shape)
    if (max(d[2:3]) <= TISSUE_LEVEL_MAX_PX) return(list(level = lv, arr = a))
    prev <- a
    lv <- lv + 1L
  }
}

make_tissue_mask <- function(tiff, dapi_plane, sizeX, sizeY) {
  pl <- .pick_mask_level(tiff)
  m <- reticulate::py_to_r(.np$asarray(pl$arr$`__getitem__`(as.integer(dapi_plane - 1L))))
  m <- t(m)                                   # [x, y], like img_read()
  storage.mode(m) <- "double"

  scale_x <- sizeX / nrow(m)
  scale_y <- sizeY / ncol(m)
  px_l <- PX_UM * scale_x

  acq <- m > 0                                # non-acquired FOVs are exactly 0
  glass <- min(stats::quantile(m[acq], 0.05, names = FALSE), 40)

  b <- gblur(Image(m), sigma = max(1, 9 / px_l))
  tm <- (b > TISSUE_GLASS_MULT * glass) & acq
  tm <- fillHull(closing(tm, makeBrush(2L * round(10 / px_l) + 1L, "disc")))

  lab <- bwlabel(tm)
  ar  <- tabulate(as.integer(lab))
  keep_obj <- which(ar >= TISSUE_MIN_OBJ_MM2 * 1e6 / px_l^2)
  tm <- matrix(as.integer(lab) %in% keep_obj, nrow(m), ncol(m))

  list(
    mask = tm,
    level = pl$level,
    scale_x = scale_x,
    scale_y = scale_y,
    px_um_level = px_l,
    glass = glass,
    tissue_mm2 = sum(tm) * px_l * PX_UM * scale_y / 1e6
  )
}

## Full-resolution tissue mask for a window, by nearest-neighbour upsampling.
.tissue_window <- function(tmk, x0, y0, w, h) {
  ix <- pmin(floor((x0 + seq_len(w) - 1L) / tmk$scale_x) + 1L, nrow(tmk$mask))
  iy <- pmin(floor((y0 + seq_len(h) - 1L) / tmk$scale_y) + 1L, ncol(tmk$mask))
  tmk$mask[ix, iy, drop = FALSE]
}

.tile_has_tissue <- function(tmk, x0, y0, w, h) {
  ix <- unique(pmin(floor((x0 + seq(0, w - 1, by = 4)) / tmk$scale_x) + 1L, nrow(tmk$mask)))
  iy <- unique(pmin(floor((y0 + seq(0, h - 1, by = 4)) / tmk$scale_y) + 1L, ncol(tmk$mask)))
  any(tmk$mask[ix, iy])
}

## Nuclear segmentation parameters (px at 0.28 um/px).
if (!exists("NUC_BG_SIGMA")) NUC_BG_SIGMA <- 25      # local background blur
if (!exists("NUC_OFFSET")) NUC_OFFSET <- 0.03        # contrast above background
if (!exists("NUC_WS_TOL")) NUC_WS_TOL <- 1           # watershed split tolerance

###############################################################################
## TILE PROCESSING
###############################################################################

## Per-tile constants of the pixel correction (background, noise, blank fit)
## are identical for every cell of a tile. They are returned once per tile in
## a separate table (tile_qc_<sid>.csv) instead of being repeated on every
## cell, which kept ~150 redundant columns in memory for millions of cells.
.TILE_CONST <- c("bg", "pix_thr", "noise", "blank_intercept", "blank_slope", "blank_r2")

.feat_cells <- function(cc, lab, valid_lab, nmax, prefix) {
  f <- .aggregate_corrected_features(cc, lab, valid_lab, nmax)
  const <- f[1, .TILE_CONST, drop = FALSE]
  f <- f[, c("sub", "frac_pos", "af_frac"), drop = FALSE]
  colnames(f) <- paste0(prefix, "_", colnames(f))
  colnames(const) <- paste0(prefix, "_", colnames(const))
  list(cells = f, const = const)
}

process_tile <- function(tiff, cm, x0, y0, w, h, sid, tmk, sizeX, sizeY) {

  ## Core [x0, x0+w) x [y0, y0+h); read window adds HALO_PX on every side.
  rx0 <- max(0L, x0 - HALO_PX)
  ry0 <- max(0L, y0 - HALO_PX)
  rx1 <- min(sizeX, x0 + w + HALO_PX)
  ry1 <- min(sizeY, y0 + h + HALO_PX)
  rw  <- rx1 - rx0
  rh  <- ry1 - ry0

  ## Registration QC on the core.
  reg_qc <- .dapi_registration_qc(tiff, cm, x0, y0, w, h)

  ## DAPI for segmentation.
  dapi_plane <- if (!is.null(cm$DAPI_AR)) cm$DAPI_AR else cm$DAPI_all[1]
  dapi <- img_read(tiff, dapi_plane, rx0, ry0, rw, rh)

  mx <- max(dapi, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(NULL)

  tis <- Image(.tissue_window(tmk, rx0, ry0, rw, rh))
  if (sum(tis) < 2000) return(NULL)

  d <- Image(dapi / mx)
  d_s <- gblur(d, sigma = 1)

  ## Nuclear segmentation. v11.3 used a single Otsu threshold per tile, which
  ## kept only the brightest nuclei (~1,500 cells/mm2 in dense tumour) and
  ## missed dim stromal/epithelial nuclei. Now: sqrt-compressed DAPI minus the
  ## glass offset, local-contrast threshold against a broad background blur,
  ## then a smoothed distance-map watershed (restricted to foreground) to
  ## split touching nuclei.
  s <- Image(sqrt(pmax(dapi - tmk$glass, 0)))
  s_hi <- stats::quantile(s[tis > 0], 0.999, names = FALSE)
  if (!is.finite(s_hi) || s_hi <= 0) return(NULL)
  ss <- gblur(s / s_hi, sigma = 1.5)
  fg <- ((ss - gblur(ss, sigma = NUC_BG_SIGMA)) > NUC_OFFSET) & tis
  rm(s, ss)
  fg <- fillHull(opening(fg, makeBrush(5, "disc")))

  nuc <- watershed(gblur(distmap(fg), sigma = 1) * fg, tolerance = NUC_WS_TOL, ext = 3)
  rm(fg)

  sh0 <- computeFeatures.shape(nuc)
  if (is.null(sh0) || !nrow(sh0)) return(NULL)

  drop <- which(sh0[, "s.area"] < MIN_AREA_PX | sh0[, "s.area"] > MAX_AREA_PX)
  if (length(drop)) nuc <- rmObjects(nuc, drop, reenumerate = TRUE)

  nmax <- max(nuc)
  if (nmax < 10) return(NULL)

  ## Compartments:
  ##   - eroded nucleus for CEBPg / AR
  ##   - local perinuclear rings for PDGF / panCK / FOLH1 / DLL3
  nuc_mask <- nuc > 0
  nuc_ero <- nuc * as.integer(erode(nuc_mask, makeBrush(3, "disc")))

  rings <- .make_ring_labels(nuc, nuc_mask, tis, d_s)
  union_ring <- Reduce(`|`, lapply(rings, function(x) x > 0))

  rm(tis, d_s, d, nuc_mask)

  ## Control planes are read once per tile and shared across markers.
  ctrl_needed <- unique(c(unname(CTRL_OF), unname(BLANK_REF_OF[unique(CTRL_OF)])))
  ctrl_img <- lapply(setNames(ctrl_needed, ctrl_needed), function(cn) {
    img_read(tiff, cm[[cn]][1], rx0, ry0, rw, rh)
  })

  corr_list  <- list()
  const_list <- list()

  for (m in MARKERS_CORE) {
    cp  <- MARKER_COMPARTMENT[[m]]
    ctl <- CTRL_OF[[m]]

    mk     <- img_read(tiff, cm[[m]][1], rx0, ry0, rw, rh)
    blk    <- ctrl_img[[ctl]]
    af_ref <- ctrl_img[[BLANK_REF_OF[[ctl]]]]

    if (cp == "nuc") {
      lab <- nuc_ero
      valid_lab <- nuc_ero > 0

      cc_marker <- .correct_marker_pixels(mk, blk, valid_lab & is.finite(mk) & is.finite(blk))
      cc_blank  <- .correct_marker_pixels(blk, af_ref, valid_lab & is.finite(blk) & is.finite(af_ref))

      pre <- .feature_prefix(m, cp)
      f1 <- .feat_cells(cc_marker, lab, valid_lab, nmax, pre)
      f2 <- .feat_cells(cc_blank,  lab, valid_lab, nmax, paste0(pre, "_blank"))

      corr_list  <- c(corr_list,  list(f1$cells, f2$cells))
      const_list <- c(const_list, list(f1$const, f2$const))

    } else {
      cc_marker <- .correct_marker_pixels(mk, blk, union_ring & is.finite(mk) & is.finite(blk))
      cc_blank  <- .correct_marker_pixels(blk, af_ref, union_ring & is.finite(blk) & is.finite(af_ref))

      for (r_um in RING_UM_SET) {
        lab <- rings[[.ring_tag(r_um)]]
        valid_lab <- lab > 0
        pre <- .feature_prefix(m, cp, r_um)

        f1 <- .feat_cells(cc_marker, lab, valid_lab, nmax, pre)
        f2 <- .feat_cells(cc_blank,  lab, valid_lab, nmax, paste0(pre, "_blank"))

        corr_list <- c(corr_list, list(f1$cells, f2$cells))
        ## Correction constants do not depend on ring radius: keep primary only.
        if (r_um == RING_UM_PRIMARY) const_list <- c(const_list, list(f1$const, f2$const))
      }
    }

    rm(mk, cc_marker, cc_blank)
  }

  rm(ctrl_img, union_ring)

  corr_df <- do.call(cbind, corr_list)
  rm(corr_list)

  dv <- .mean_by_label_valid(dapi, nuc, nuc > 0, nmax)
  rm(dapi)

  mo <- computeFeatures.moment(nuc)
  sh <- computeFeatures.shape(nuc)

  out <- data.frame(
    slide_id = sid,
    x_px = mo[, "m.cx"] + rx0,
    y_px = mo[, "m.cy"] + ry0,
    area = sh[, "s.area"],
    area_um2 = sh[, "s.area"] * PX_UM^2,
    perim = sh[, "s.perimeter"],
    dapi_nuc = dv,
    corr_df,
    reg_cor = reg_qc$reg_cor,
    reg_qc_informative = reg_qc$reg_qc_informative,
    tile_x = x0,
    tile_y = y0,
    stringsAsFactors = FALSE
  )

  rm(nuc, nuc_ero, rings, mo, sh, corr_df)

  ## Keep only cells whose centroid lies in the core (halo cells belong to
  ## the neighbouring tile).
  out <- out[is.finite(out$x_px) &
               out$x_px >= x0 & out$x_px < x0 + w &
               out$y_px >= y0 & out$y_px < y0 + h, , drop = FALSE]

  tile_row <- data.frame(
    slide_id = sid,
    tile_x = x0, tile_y = y0, tile_w = w, tile_h = h,
    n_cells = nrow(out),
    as.data.frame(reg_qc),
    do.call(cbind, const_list),
    stringsAsFactors = FALSE
  )

  list(cells = out, tile = tile_row)
}

###############################################################################
## RUN
###############################################################################

if (!exists("RUN_SLIDES")) RUN_SLIDES <- slides$slide_id

stopifnot(all(RUN_SLIDES %in% slides$slide_id))

if (!is.finite(MAX_TILES) && READ_MODE != "tifffile+zarr") {
  stop(
    "Full run needs the tifffile+zarr reader (pyramid levels + chunked reads).",
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

.ckpt_file <- function(sid) {
  file.path(OUT, paste0(if (is.finite(MAX_TILES)) "cells_raw_SMOKETEST_" else "cells_raw_", sid, ".rds"))
}

segment_slide <- function(k) {
  sid <- slides$slide_id[k]
  cm  <- info[[sid]]$chmap
  sizeX <- slides$sizeX[k]
  sizeY <- slides$sizeY[k]

  ckpt <- .ckpt_file(sid)
  if (file.exists(ckpt) && !FORCE_RERUN) {
    message("=== ", sid, ": re-using checkpoint ", basename(ckpt))
    return(invisible(ckpt))
  }

  dapi_plane <- if (!is.null(cm$DAPI_AR)) cm$DAPI_AR else cm$DAPI_all[1]
  tmk <- make_tissue_mask(slides$tiff[k], dapi_plane, sizeX, sizeY)
  saveRDS(tmk, file.path(OUT, paste0("tissue_mask_", sid, ".rds")))
  message(
    "=== ", sid, " tissue mask: level ", tmk$level, " (", round(tmk$px_um_level, 2),
    " um/px), glass ", round(tmk$glass), ", tissue ", round(tmk$tissue_mm2, 1), " mm2"
  )

  ## Full coverage: the last row/column of tiles is a partial tile.
  xs <- seq(0L, sizeX - 1L, by = TILE_SIZE)
  ys <- seq(0L, sizeY - 1L, by = TILE_SIZE)
  g <- expand.grid(x = xs, y = ys)
  g$w <- pmin(TILE_SIZE, sizeX - g$x)
  g$h <- pmin(TILE_SIZE, sizeY - g$y)
  n_total <- nrow(g)

  g <- g[vapply(seq_len(nrow(g)), function(i) {
    .tile_has_tissue(tmk, g$x[i], g$y[i], g$w[i], g$h[i])
  }, logical(1)), , drop = FALSE]

  if (is.finite(MAX_TILES)) {
    mid <- order((g$x - sizeX / 2)^2 + (g$y - sizeY / 2)^2)
    g <- g[head(mid, MAX_TILES), , drop = FALSE]
  }

  message(
    "=== ", sid, " [", slides$group[k], "] : ", nrow(g), " tissue tiles of ",
    n_total, " (", TILE_SIZE, " px core + ", HALO_PX, " px halo)"
  )

  t0 <- Sys.time()
  cells_l <- vector("list", nrow(g))
  tiles_l <- vector("list", nrow(g))

  for (i in seq_len(nrow(g))) {
    ti <- Sys.time()

    r <- tryCatch(
      process_tile(slides$tiff[k], cm, g$x[i], g$y[i], g$w[i], g$h[i],
                   sid, tmk, sizeX, sizeY),
      error = function(e) {
        message("   tile ", i, " (x=", g$x[i], " y=", g$y[i], ") FAILED: ",
                conditionMessage(e))
        NULL
      }
    )

    dt <- as.numeric(difftime(Sys.time(), ti, units = "secs"))
    .tile_times <<- c(.tile_times, dt)

    if (!is.null(r)) {
      cells_l[[i]] <- r$cells
      tiles_l[[i]] <- cbind(r$tile, tile_sec = dt)
    }

    if (i == 1 || i %% VERBOSE_EVERY == 0 || i == nrow(g)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf(
        "   %s tile %d/%d [%.1f s] -> %s cells | elapsed %.1f min, ETA %.1f min",
        sid, i, nrow(g), dt, if (is.null(r)) 0 else nrow(r$cells),
        el, el / i * (nrow(g) - i)
      ))
    }

    if (i %% 20 == 0) gc(FALSE)
  }

  out <- dplyr::bind_rows(cells_l)
  tiles <- dplyr::bind_rows(tiles_l)
  rm(cells_l, tiles_l)

  if (nrow(out)) out$cell_id <- paste0(sid, "_", seq_len(nrow(out)))

  message(
    "   -> ", nrow(out), " cells in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "min")), 1), " min"
  )

  write.csv(tiles, file.path(OUT, paste0("tile_qc_", sid, ".csv")), row.names = FALSE)
  saveRDS(out, ckpt)
  rm(out, tiles)
  gc(FALSE)

  invisible(ckpt)
}

for (k in which(slides$slide_id %in% RUN_SLIDES)) segment_slide(k)

if (!COMBINE) {
  message("\n02_segment_measure.R: segmented ", paste(RUN_SLIDES, collapse = ", "),
          " (COMBINE = FALSE, not combining).")
} else {

.missing_ckpt <- slides$slide_id[!file.exists(vapply(slides$slide_id, .ckpt_file, character(1)))]
if (length(.missing_ckpt)) {
  stop("No checkpoint yet for: ", paste(.missing_ckpt, collapse = ", "),
       " -- run 02 for these slides before combining.", call. = FALSE)
}

cells_raw <- dplyr::bind_rows(lapply(slides$slide_id, function(sid) readRDS(.ckpt_file(sid))))

tile_qc <- dplyr::bind_rows(lapply(slides$slide_id, function(sid) {
  f <- file.path(OUT, paste0("tile_qc_", sid, ".csv"))
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
}))

tissue_area <- do.call(rbind, lapply(slides$slide_id, function(sid) {
  tmk <- readRDS(file.path(OUT, paste0("tissue_mask_", sid, ".rds")))
  data.frame(slide_id = sid, tissue_mm2 = tmk$tissue_mm2, mask_level = tmk$level,
             glass_level = tmk$glass)
}))
write.csv(tissue_area, file.path(OUT, "tissue_area_by_slide.csv"), row.names = FALSE)

if (!nrow(cells_raw)) {
  stop("No cells segmented -- check tissue mask / thresholds.", call. = FALSE)
}

cells_raw <- left_join(
  cells_raw,
  slides[, c("slide_id", "group", "run", "comet_sw", "has_negctrl", "tier2", "tier3")],
  by = "slide_id"
)

out_file <- if (is.finite(MAX_TILES)) "cells_raw_SMOKETEST.rds" else "cells_raw.rds"
saveRDS(cells_raw, file.path(OUT, out_file))
write.csv(tile_qc, file.path(OUT, "tile_qc_all.csv"), row.names = FALSE)

###############################################################################
## QC OUTPUTS
###############################################################################

qc_cells <- cells_raw |>
  group_by(slide_id, group) |>
  summarise(
    n_cells = n(),
    median_area_px = .median_or_na(area),
    median_area_um2 = .median_or_na(area_um2),
    .groups = "drop"
  )

qc_tiles <- tile_qc |>
  group_by(slide_id) |>
  summarise(
    n_tiles = n(),
    n_reg_qc_informative_tiles = sum(reg_qc_informative %in% TRUE, na.rm = TRUE),
    n_informative_dapi_pairs = sum(reg_cor_n, na.rm = TRUE),
    median_dapi_cor = .median_or_na(reg_cor_median),
    p05_dapi_cor = .quantile_or_na(reg_cor_median, 0.05),
    min_dapi_cor = .min_or_na(reg_cor_min),
    frac_low_dapi_cor = suppressWarnings(sum(reg_cor_n_low, na.rm = TRUE) /
                                           pmax(sum(reg_cor_n, na.rm = TRUE), 1)),
    median_ref_sd = .median_or_na(reg_ref_sd),
    median_ref_range = .median_or_na(reg_ref_range),
    total_tile_hours = sum(tile_sec, na.rm = TRUE) / 3600,
    .groups = "drop"
  )

qc <- left_join(qc_cells, qc_tiles, by = "slide_id") |>
  left_join(tissue_area, by = "slide_id") |>
  mutate(cells_per_mm2_tissue = n_cells / pmax(tissue_mm2, 1e-9))

cat("\nSegmentation QC:\n")
print(as.data.frame(qc))

write.csv(qc, file.path(OUT, "qc_registration.csv"), row.names = FALSE)

## Corrected primary-feature audit. Per-cell feature summaries come from
## cells_raw; per-tile correction constants (blank fit, background) from
## tile_qc, summarised as the median over tissue tiles.
.tile_med <- function(sid, col) .median_or_na(.col_or_na(tile_qc[tile_qc$slide_id == sid, , drop = FALSE], col))

corr_audit <- do.call(rbind, lapply(unique(cells_raw$slide_id), function(sid) {
  d <- cells_raw[cells_raw$slide_id == sid, ]

  do.call(rbind, lapply(MARKERS_CORE, function(m) {
    cp <- MARKER_COMPARTMENT[[m]]
    pre <- .feature_prefix(m, cp)
    col <- function(x) .col_or_na(d, paste0(pre, x))

    data.frame(
      slide_id = sid,
      marker = m,
      compartment = cp,
      control = CTRL_OF[[m]],
      radius_um = ifelse(cp == "ring", RING_UM_PRIMARY, NA),

      median_sub = round(.median_or_na(col("_sub")), 3),
      p95_sub = round(.qsafe(col("_sub"), .95), 3),
      p99_sub = round(.qsafe(col("_sub"), .99), 3),
      median_frac_pos = round(.median_or_na(col("_frac_pos")), 4),
      p95_frac_pos = round(.qsafe(col("_frac_pos"), .95), 4),
      p99_frac_pos = round(.qsafe(col("_frac_pos"), .99), 4),
      median_af_frac = round(.median_or_na(col("_af_frac")), 4),
      p99_af_frac = round(.qsafe(col("_af_frac"), .99), 4),
      median_tile_bg = round(.tile_med(sid, paste0(pre, "_bg")), 3),
      median_tile_noise = round(.tile_med(sid, paste0(pre, "_noise")), 3),
      median_blank_intercept = round(.tile_med(sid, paste0(pre, "_blank_intercept")), 4),
      median_blank_slope = round(.tile_med(sid, paste0(pre, "_blank_slope")), 4),
      median_blank_r2 = round(.tile_med(sid, paste0(pre, "_blank_r2")), 4),

      blank_median_sub = round(.median_or_na(col("_blank_sub")), 3),
      blank_p95_sub = round(.qsafe(col("_blank_sub"), .95), 3),
      blank_p99_sub = round(.qsafe(col("_blank_sub"), .99), 3),
      blank_median_frac_pos = round(.median_or_na(col("_blank_frac_pos")), 4),
      blank_p95_frac_pos = round(.qsafe(col("_blank_frac_pos"), .95), 4),
      blank_p99_frac_pos = round(.qsafe(col("_blank_frac_pos"), .99), 4),
      blank_median_af_frac = round(.median_or_na(col("_blank_af_frac")), 4),
      blank_p99_af_frac = round(.qsafe(col("_blank_af_frac"), .99), 4),

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
      col <- function(x) .col_or_na(d, paste0(pre, x))

      data.frame(
        slide_id = sid,
        marker = m,
        radius_um = r_um,

        median_sub = round(.median_or_na(col("_sub")), 3),
        p95_sub = round(.qsafe(col("_sub"), .95), 3),
        p99_sub = round(.qsafe(col("_sub"), .99), 3),
        median_frac_pos = round(.median_or_na(col("_frac_pos")), 4),
        p95_frac_pos = round(.qsafe(col("_frac_pos"), .95), 4),
        p99_frac_pos = round(.qsafe(col("_frac_pos"), .99), 4),
        median_af_frac = round(.median_or_na(col("_af_frac")), 4),
        p99_af_frac = round(.qsafe(col("_af_frac"), .99), 4),

        blank_median_sub = round(.median_or_na(col("_blank_sub")), 3),
        blank_p95_sub = round(.qsafe(col("_blank_sub"), .95), 3),
        blank_p99_sub = round(.qsafe(col("_blank_sub"), .99), 3),
        blank_median_frac_pos = round(.median_or_na(col("_blank_frac_pos")), 4),
        blank_p99_frac_pos = round(.qsafe(col("_blank_frac_pos"), .99), 4),

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

bad_reg <- qc$slide_id[
  (is.finite(qc$median_dapi_cor) & qc$median_dapi_cor < REG_LOW_COR_THRESHOLD) |
    (is.finite(qc$frac_low_dapi_cor) & qc$frac_low_dapi_cor > 0.25)
]
if (length(bad_reg)) {
  warning(
    "Low informative-window cross-cycle DAPI registration QC on: ",
    paste(bad_reg, collapse = ", "),
    call. = FALSE
  )
}

no_reg_info <- qc$slide_id[qc$n_reg_qc_informative_tiles %in% 0]
if (length(no_reg_info)) {
  warning(
    "No DAPI-informative registration QC windows on: ",
    paste(no_reg_info, collapse = ", "),
    " -- registration QC is unknown, not failed.",
    call. = FALSE
  )
}

if (any(qc$median_area_um2 < 15 | qc$median_area_um2 > 150, na.rm = TRUE)) {
  warning(
    "Median nuclear area outside 15-150 um^2 on: ",
    paste(qc$slide_id[qc$median_area_um2 < 15 | qc$median_area_um2 > 150], collapse = ", "),
    call. = FALSE
  )
}

if (length(.tile_times)) {
  message("\nMedian tile time (this session): ", round(median(.tile_times), 1), " s")
  message("Total tiles (this session): ", length(.tile_times))
}

message("\n02_segment_measure.R complete -> ", file.path(OUT, out_file))

}  # end if (COMBINE)
