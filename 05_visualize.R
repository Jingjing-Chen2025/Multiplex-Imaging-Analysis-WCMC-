###############################################################################
# 05_visualize.R  v1.0
#
# Per-biopsy visualisation of ASPC, epithelial cancer (EPI_CA), neuroendocrine
# cancer (NE_CA) cells and the ASPC neighbourhood, on the image itself.
#
# INPUT (all in OUT)
#   cells_spatial.rds          per-cell calls + neighbourhood flags (04)
#   channels.rds               channel map + pixel size (01)
#   tissue_mask_<sid>.rds      slide tissue mask (02)
#   aspc_neighbourhood_composition.csv (04)
#   population_counts_by_slide_primary_6um.csv (03)
#   original COMET OME-TIFFs   read through tifffile/zarr pyramid levels
#
# OUTPUT  OUT/plots/per_slide/<sid>/
#   <sid>_A_marker_overview.png    control-subtracted marker composite + panels
#   <sid>_B_celltype_map.png       ASPC / EPI_CA / NE_CA / Multi on the tissue
#   <sid>_C_aspc_neighbourhood.png ASPC, 50 um neighbourhood, neighbours by type
#   <sid>_D_hotspot_<k>.png        full-resolution crops of top ASPC hotspots
#   <sid>_E_nbhd_composition.png   neighbourhood vs whole-slide composition
#   cells_<sid>.csv                per-cell x/y (um), call, neighbourhood flags
#                                  (for QuPath / napari)
#         OUT/plots/
#   05_summary_*.png, 05_all_slides.pdf
#
# Display images are control-subtracted: marker - (a + b * matched secondary-
# only control), fit on tissue pixels, clipped at 0 -- the same matched
# controls as 02. Display only; calls come from 02/03.
###############################################################################

suppressPackageStartupMessages({
  library(reticulate)
  library(EBImage)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

if (!exists("OUT")) OUT <- "/Volumes/T9/results_ASPC_NEPC"
if (!exists("RADIUS_UM")) RADIUS_UM <- 50
if (!exists("DISPLAY_PX")) DISPLAY_PX <- 2000L      # long edge of overview images
if (!exists("HOTSPOT_UM")) HOTSPOT_UM <- 500
if (!exists("N_HOTSPOTS")) N_HOTSPOTS <- 3L

use_virtualenv("r-reticulate", required = TRUE)
.tf <- import("tifffile", convert = FALSE)
.np <- import("numpy", convert = FALSE)
.zr <- import("zarr", convert = FALSE)
.py <- import_builtins()

PLOTDIR <- file.path(OUT, "plots")
VISDIR  <- file.path(PLOTDIR, "per_slide")
dir.create(VISDIR, showWarnings = FALSE, recursive = TRUE)

z <- readRDS(file.path(OUT, "channels.rds"))
PX_UM <- if (!is.null(z$px_um)) z$px_um else 12500 / 44643
slides <- readRDS(file.path(OUT, "slides_final.rds"))
cells <- readRDS(file.path(OUT, "cells_spatial.rds"))

nbr_comp_file <- file.path(OUT, "aspc_neighbourhood_composition.csv")
nbr_comp <- if (file.exists(nbr_comp_file)) read.csv(nbr_comp_file, stringsAsFactors = FALSE) else NULL

## Matched secondary-only controls (same as 02).
CTRL_OF <- c(CEBPg = "NC3_Cy5", FOLH1 = "NC5_Cy5", AR = "NC5_Cy5",
             DLL3 = "NC5_TRITC", PDGF = "NC3_Cy5", panCK = "NC3_TRITC")
MARKER_LABEL <- c(PDGF = "PDGFa", panCK = "pan-CK", DLL3 = "DLL3",
                  AR = "AR", FOLH1 = "FOLH1 (PSMA)", CEBPg = "CEBPg")

## Population colours. Maps are scatter plots, so only three categorical
## slots are used (validated all-pairs, CVD-safe); Multi and Other are neutral.
POP_LEVELS <- c("ASPC", "EPI_CA", "NE_CA", "Multi", "Other")
POP_COLS <- c(
  EPI_CA = "#2a78d6",   # slot 1 blue
  ASPC   = "#eb6834",   # slot 2 orange
  NE_CA  = "#1baf7a",   # slot 3 aqua
  Multi  = "#52514e",   # neutral dark
  Other  = "#c9c8c2"    # neutral light
)
POP_LABEL <- c(ASPC = "ASPC (PDGFa+)", EPI_CA = "Epithelial ca. (panCK+ FOLH1/AR+)",
               NE_CA = "NE ca. (DLL3+)", Multi = "Multi (>1 call)", Other = "Other")

## Fluorescence composite colours (dark background, additive):
## PDGFa magenta, pan-CK green, DLL3 yellow, DAPI dim blue.
COMP <- list(PDGF = c(1, 0, 1), panCK = c(0, 1, 0), DLL3 = c(1, 0.85, 0))

###############################################################################
## IMAGE ACCESS
###############################################################################

.group_cache <- new.env(parent = emptyenv())

.level <- function(tiff, lev) {
  k <- paste0(tiff, "#", lev)
  if (!is.null(.group_cache[[k]])) return(.group_cache[[k]])
  st <- .tf$TiffFile(tiff)$series[0]$aszarr()
  g <- .zr$open(store = st, mode = "r")
  a <- g$`__getitem__`(as.character(lev))
  assign(k, a, envir = .group_cache)
  a
}

## Read plane (1-based) at a level, optional window in that level's pixels.
## Returns [x, y] double matrix.
read_plane <- function(tiff, plane, lev = 0L, x0 = NULL, y0 = NULL, w = NULL, h = NULL) {
  a <- .level(tiff, lev)
  idx <- if (is.null(w)) {
    tuple(as.integer(plane - 1L), .py$slice(NULL), .py$slice(NULL))
  } else {
    tuple(as.integer(plane - 1L),
          .py$slice(as.integer(y0), as.integer(y0 + h)),
          .py$slice(as.integer(x0), as.integer(x0 + w)))
  }
  m <- t(py_to_r(.np$asarray(a$`__getitem__`(idx))))
  storage.mode(m) <- "double"
  m
}

## Marker minus scaled matched control, fit on tissue pixels.
ctrl_subtract <- function(mk, ctl, tis) {
  ok <- which(tis & is.finite(mk) & is.finite(ctl))
  if (length(ok) > 200000L) { set.seed(1); ok <- sample(ok, 200000L) }
  if (length(ok) < 100L) return(pmax(mk - ctl, 0))
  cf <- tryCatch(coef(MASS::rlm(mk[ok] ~ ctl[ok], maxit = 50)),
                 error = function(e) coef(lm(mk[ok] ~ ctl[ok])))
  pmax(mk - (cf[1] + cf[2] * ctl), 0)
}

stretch <- function(m, tis = NULL, hi_q = 0.998, lo = 0) {
  v <- if (is.null(tis)) m else m[tis]
  hi <- stats::quantile(v[is.finite(v)], hi_q, names = FALSE)
  if (!is.finite(hi) || hi <= lo) hi <- lo + 1
  x <- (m - lo) / (hi - lo)
  x[!is.finite(x)] <- 0
  pmin(pmax(x, 0), 1)
}

## [x, y] matrices in 0..1 -> raster (rows = y) for grid / ggplot.
to_raster <- function(r, g, b) {
  as.raster(array(c(t(r), t(g), t(b)), dim = c(ncol(r), nrow(r), 3)))
}

downsize <- function(m, w) {
  if (nrow(m) <= w) return(m)
  as.matrix(resize(Image(m), w = w, h = round(ncol(m) * w / nrow(m))))
}

###############################################################################
## PLOT HELPERS
###############################################################################

theme_map <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(colour = "#52514e", size = 8),
      axis.title = element_text(colour = "#52514e", size = 9),
      plot.title = element_text(face = "bold", colour = "#0b0b0b"),
      plot.subtitle = element_text(colour = "#52514e", size = 9),
      legend.text = element_text(colour = "#0b0b0b", size = 9),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

save_png <- function(p, file, w = 9, h = 9, dpi = 200) {
  ggsave(file, p, width = w, height = h, dpi = dpi, bg = "white")
  message("  wrote ", basename(file))
}

###############################################################################
## PER-SLIDE FIGURES
###############################################################################

hotspot_windows <- function(d, sizeX, sizeY, win_um = HOTSPOT_UM, n = N_HOTSPOTS) {
  a <- d[d$ASPC, , drop = FALSE]
  if (nrow(a) < 5) return(NULL)
  bin_px <- (win_um / 5) / PX_UM
  nx <- ceiling(sizeX / bin_px); ny <- ceiling(sizeY / bin_px)
  bx <- pmin(floor(a$x_px / bin_px) + 1, nx)
  by <- pmin(floor(a$y_px / bin_px) + 1, ny)
  H <- matrix(tabulate(bx + (by - 1) * nx, nbins = nx * ny), nx, ny)
  ## 5x5 box sum = one window
  S <- as.matrix(filter2(Image(H), matrix(1, 5, 5), boundary = 0))
  out <- list()
  for (k in seq_len(n)) {
    if (max(S) < 3) break
    ij <- which(S == max(S), arr.ind = TRUE)[1, ]
    cx <- (ij[1] - 0.5) * bin_px; cy <- (ij[2] - 0.5) * bin_px
    w <- round(win_um / PX_UM)
    x0 <- max(0, min(sizeX - w, round(cx - w / 2)))
    y0 <- max(0, min(sizeY - w, round(cy - w / 2)))
    out[[k]] <- data.frame(k = k, x0 = x0, y0 = y0, w = w, n_ASPC_bins = max(S))
    S[max(1, ij[1] - 7):min(nx, ij[1] + 7), max(1, ij[2] - 7):min(ny, ij[2] + 7)] <- 0
  }
  bind_rows(out)
}

visualize_slide <- function(k) {
  sid  <- slides$slide_id[k]
  grp  <- as.character(slides$group[k])
  tiff <- z$slides$tiff[match(sid, z$slides$slide_id)]
  cm   <- z$info[[sid]]$chmap
  sizeX <- z$slides$sizeX[match(sid, z$slides$slide_id)]
  sizeY <- z$slides$sizeY[match(sid, z$slides$slide_id)]
  d <- cells[cells$slide_id == sid, , drop = FALSE]
  if (!nrow(d)) { message("No cells for ", sid); return(invisible(NULL)) }

  sdir <- file.path(VISDIR, sid)
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  message("\n=== ", sid, " [", grp, "]  ", nrow(d), " cells")

  d$x_mm <- d$x_px * PX_UM / 1000
  d$y_mm <- d$y_px * PX_UM / 1000
  d$pop <- factor(d$pop, levels = POP_LEVELS)
  counts <- table(d$pop)

  ## ---- export for QuPath / napari -------------------------------------------
  write.csv(
    data.frame(cell_id = d$cell_id, x_um = round(d$x_px * PX_UM, 2),
               y_um = round(d$y_px * PX_UM, 2), area_um2 = round(d$area_um2, 1),
               population = d$pop, ASPC_CEBPg_pos = d$ASPC_CEBPg_pos,
               in_ASPC_nbhd = d$in_ASPC_nbhd, n_ASPC_within_r = d$n_ASPC_within_r),
    file.path(sdir, paste0("cells_", sid, ".csv")), row.names = FALSE
  )

  ## ---- low-resolution planes ------------------------------------------------
  tmk <- readRDS(file.path(OUT, paste0("tissue_mask_", sid, ".rds")))
  lev <- tmk$level
  dapi_l <- read_plane(tiff, cm$DAPI_AR, lev)
  tis_l  <- tmk$mask > 0
  sc <- tmk$scale_x                         # full-res px per level px

  disp_w <- min(DISPLAY_PX, nrow(dapi_l))
  f <- nrow(dapi_l) / disp_w                # level px per display px
  ext_mm <- c(0, sizeX * PX_UM / 1000, 0, sizeY * PX_UM / 1000)

  mk_l <- list()
  for (m in names(CTRL_OF)) {
    mk  <- read_plane(tiff, cm[[m]][1], lev)
    ctl <- read_plane(tiff, cm[[CTRL_OF[[m]]]][1], lev)
    mk_l[[m]] <- downsize(stretch(ctrl_subtract(mk, ctl, tis_l), tis_l, 0.999), disp_w)
  }
  dapi_d <- downsize(stretch(dapi_l, tis_l, 0.995), disp_w)
  tis_d  <- downsize(tis_l * 1, disp_w) > 0.5
  rm(dapi_l)

  ## ---- A. marker overview ---------------------------------------------------
  comp <- list(r = dapi_d * 0, g = dapi_d * 0, b = dapi_d * 0.35)
  for (m in names(COMP)) {
    comp$r <- comp$r + COMP[[m]][1] * mk_l[[m]]
    comp$g <- comp$g + COMP[[m]][2] * mk_l[[m]]
    comp$b <- comp$b + COMP[[m]][3] * mk_l[[m]]
  }
  comp <- lapply(comp, function(x) pmin(x, 1))
  comp_r <- to_raster(comp$r, comp$g, comp$b)

  fA <- file.path(sdir, paste0(sid, "_A_marker_overview.png"))
  asp <- ncol(dapi_d) / nrow(dapi_d)
  png(fA, width = 3000, height = round(3000 * (0.62 * asp + 0.08)), res = 200, bg = "black")
  grid.newpage()
  lay <- grid.layout(3, 4, heights = unit(c(0.07, 0.62, 0.31), "null"))
  pushViewport(viewport(layout = lay))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:4))
  grid.text(paste0(sid, " [", grp, "]  control-subtracted markers"),
            x = 0.01, just = "left", gp = gpar(col = "white", fontsize = 14, fontface = "bold"))
  grid.text("composite: PDGFa magenta · pan-CK green · DLL3 yellow · DAPI blue",
            x = 0.99, just = "right", gp = gpar(col = "#c3c2b7", fontsize = 10))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:4))
  grid.raster(comp_r, interpolate = TRUE)
  popViewport()
  pan <- c("PDGF", "panCK", "DLL3", "AR")
  for (i in seq_along(pan)) {
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = i))
    g <- mk_l[[pan[i]]]
    grid.raster(to_raster(g, g, g), interpolate = TRUE)
    grid.text(MARKER_LABEL[[pan[i]]], x = 0.04, y = 0.94, just = c("left", "top"),
              gp = gpar(col = "white", fontsize = 11, fontface = "bold"))
    popViewport()
  }
  popViewport()
  dev.off()
  message("  wrote ", basename(fA))

  ## Light tissue background for maps: white glass, pale grey tissue with a
  ## little DAPI texture.
  bg_v <- ifelse(tis_d, 0.93 - 0.18 * dapi_d, 1)
  bg_r <- to_raster(bg_v, bg_v, bg_v)

  ## ---- B. cell-type map -----------------------------------------------------
  called <- d[d$pop != "Other", , drop = FALSE]
  called <- called[order(called$pop == "ASPC"), ]        # ASPC drawn last
  pB <- ggplot(called, aes(x_mm, y_mm, colour = pop)) +
    annotation_raster(bg_r, ext_mm[1], ext_mm[2], ext_mm[3], ext_mm[4]) +
    geom_point(size = 0.12, alpha = 0.75, stroke = 0, shape = 16) +
    scale_colour_manual(
      values = POP_COLS, drop = FALSE,
      labels = function(b) paste0(POP_LABEL[b], "  n=", format(as.integer(counts[b]), big.mark = ","))
    ) +
    scale_y_reverse(limits = c(ext_mm[4], 0), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, ext_mm[2]), expand = c(0, 0)) +
    coord_fixed() +
    labs(title = paste0(sid, " [", grp, "]  cell types"),
         subtitle = paste0(format(nrow(d), big.mark = ","), " segmented cells; ",
                           "Other (", format(as.integer(counts["Other"]), big.mark = ","),
                           ") shown as tissue background"),
         x = "mm", y = "mm") +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
    theme_map()
  save_png(pB, file.path(sdir, paste0(sid, "_B_celltype_map.png")), 9, 9.6)

  ## ---- C. ASPC neighbourhood ------------------------------------------------
  a <- d[d$ASPC, , drop = FALSE]
  if (nrow(a)) {
    ## union of RADIUS_UM discs, rasterised at display resolution
    disp_um <- PX_UM * sc * f
    pts <- matrix(0, nrow(dapi_d), ncol(dapi_d))
    pi_ <- cbind(pmin(floor(a$x_px / (sc * f)) + 1, nrow(pts)),
                 pmin(floor(a$y_px / (sc * f)) + 1, ncol(pts)))
    pts[pi_] <- 1
    nb_d <- as.matrix(distmap(Image(pts == 0))) <= RADIUS_UM / disp_um
    nb_v_r <- ifelse(nb_d, 0.99, bg_v)
    nb_v_g <- ifelse(nb_d, 0.86, bg_v)
    nb_v_b <- ifelse(nb_d, 0.76, bg_v)
    nb_r <- to_raster(nb_v_r, nb_v_g, nb_v_b)

    nbd <- d[d$in_ASPC_nbhd & !d$ASPC, , drop = FALSE]
    nbd <- nbd[order(nbd$pop != "Other"), ]
    ncnt <- table(nbd$pop)
    pC <- ggplot(mapping = aes(x_mm, y_mm, colour = pop)) +
      annotation_raster(nb_r, ext_mm[1], ext_mm[2], ext_mm[3], ext_mm[4]) +
      geom_point(data = nbd, size = 0.15, alpha = 0.8, stroke = 0, shape = 16) +
      geom_point(data = a, size = 0.35, alpha = 1, stroke = 0, shape = 16) +
      scale_colour_manual(
        values = POP_COLS, drop = FALSE,
        labels = function(b) ifelse(
          b == "ASPC",
          paste0("ASPC  n=", format(nrow(a), big.mark = ",")),
          paste0(POP_LABEL[b], " in nbhd  n=", format(as.integer(ncnt[b]), big.mark = ","))
        )
      ) +
      scale_y_reverse(limits = c(ext_mm[4], 0), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, ext_mm[2]), expand = c(0, 0)) +
      coord_fixed() +
      labs(title = paste0(sid, " [", grp, "]  ASPC neighbourhood (", RADIUS_UM, " µm)"),
           subtitle = paste0("peach = within ", RADIUS_UM, " µm of an ASPC; ",
                             format(sum(d$in_ASPC_nbhd), big.mark = ","), " cells (",
                             round(100 * mean(d$in_ASPC_nbhd), 1), "% of slide) in the neighbourhood"),
           x = "mm", y = "mm") +
      guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1))) +
      theme_map()
    save_png(pC, file.path(sdir, paste0(sid, "_C_aspc_neighbourhood.png")), 9, 9.6)
  }

  ## ---- D. hotspot zooms at full resolution -----------------------------------
  hs <- hotspot_windows(d, sizeX, sizeY)
  if (!is.null(hs) && nrow(hs)) {
    write.csv(transform(hs, x0_um = x0 * PX_UM, y0_um = y0 * PX_UM, w_um = w * PX_UM),
              file.path(sdir, paste0(sid, "_D_hotspots.csv")), row.names = FALSE)
    for (i in seq_len(nrow(hs))) {
      x0 <- hs$x0[i]; y0 <- hs$y0[i]; w <- hs$w[i]
      dp <- read_plane(tiff, cm$DAPI_AR, 0L, x0, y0, w, w)
      tis0 <- dp > 2 * tmk$glass
      ch <- lapply(names(COMP), function(m) {
        mk <- read_plane(tiff, cm[[m]][1], 0L, x0, y0, w, w)
        ctl <- read_plane(tiff, cm[[CTRL_OF[[m]]]][1], 0L, x0, y0, w, w)
        stretch(gblur(Image(ctrl_subtract(mk, ctl, tis0)), 0.8), tis0, 0.998)
      })
      names(ch) <- names(COMP)
      dpn <- stretch(dp, tis0, 0.995)
      rr <- 0 * dpn; gg <- 0 * dpn; bb <- 0.45 * dpn
      for (m in names(COMP)) {
        rr <- rr + COMP[[m]][1] * ch[[m]]
        gg <- gg + COMP[[m]][2] * ch[[m]]
        bb <- bb + COMP[[m]][3] * ch[[m]]
      }
      ds <- function(x) downsize(pmin(as.matrix(x), 1), 1000L)
      img <- to_raster(ds(rr), ds(gg), ds(bb))

      dz <- d[d$x_px >= x0 & d$x_px < x0 + w & d$y_px >= y0 & d$y_px < y0 + w &
                d$pop != "Other", , drop = FALSE]
      dz$xu <- (dz$x_px - x0) * PX_UM
      dz$yu <- (dz$y_px - y0) * PX_UM
      wu <- w * PX_UM
      zc <- table(dz$pop)

      base <- ggplot() +
        annotation_raster(img, 0, wu, 0, wu) +
        scale_y_reverse(limits = c(wu, 0), expand = c(0, 0)) +
        scale_x_continuous(limits = c(0, wu), expand = c(0, 0)) +
        coord_fixed() + labs(x = "µm", y = "µm") + theme_map() +
        theme(plot.title = element_text(size = 10))
      p1 <- base + labs(title = "PDGFa magenta · pan-CK green · DLL3 yellow · DAPI blue")
      p2 <- base +
        geom_point(data = dz, aes(xu, yu, colour = pop), shape = 21, size = 2.2, stroke = 0.6) +
        scale_colour_manual(values = POP_COLS, drop = TRUE,
                            labels = function(b) paste0(b, " n=", as.integer(zc[b]))) +
        labs(title = "cell calls (rings at nuclear centroids)")

      fD <- file.path(sdir, paste0(sid, "_D_hotspot_", i, ".png"))
      png(fD, width = 3600, height = 2000, res = 200, bg = "white")
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.06, 0.94), "null"))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
      grid.text(sprintf("%s [%s]  ASPC hotspot %d  (x %.2f–%.2f mm, y %.2f–%.2f mm)",
                        sid, grp, i, x0 * PX_UM / 1000, (x0 + w) * PX_UM / 1000,
                        y0 * PX_UM / 1000, (y0 + w) * PX_UM / 1000),
                x = 0.01, just = "left", gp = gpar(fontsize = 13, fontface = "bold"))
      popViewport()
      print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
      print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
      popViewport()
      dev.off()
      message("  wrote ", basename(fD))
    }
  }

  ## ---- E. neighbourhood composition ----------------------------------------
  if (!is.null(nbr_comp) && sid %in% nbr_comp$slide_id) {
    nc <- nbr_comp[nbr_comp$slide_id == sid, ]
    ed <- rbind(
      data.frame(where = paste0("Within ", RADIUS_UM, " µm of ASPC"),
                 population = nc$population, frac = nc$frac_in_nbhd),
      data.frame(where = "Whole slide", population = nc$population, frac = nc$frac_whole_slide)
    )
    ed$population <- factor(ed$population, levels = rev(POP_LEVELS))
    ed$where <- factor(ed$where, levels = c("Whole slide", paste0("Within ", RADIUS_UM, " µm of ASPC")))
    pE <- ggplot(ed, aes(frac, where, fill = population)) +
      geom_col(width = 0.6, colour = "white", linewidth = 0.5) +
      geom_text(aes(label = ifelse(frac >= 0.04, paste0(round(100 * frac), "%"), "")),
                position = position_stack(vjust = 0.5), size = 3,
                colour = ifelse(ed$population %in% c("Other"), "#0b0b0b", "white")) +
      scale_fill_manual(values = POP_COLS, breaks = POP_LEVELS, labels = POP_LABEL[POP_LEVELS]) +
      scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
      labs(title = paste0(sid, " [", grp, "]  who surrounds ASPC?"),
           subtitle = paste0("Neighbourhood = cells within ", RADIUS_UM, " µm of another ASPC; ",
                             nc$n_ASPC[1], " ASPC"),
           x = "share of cells", y = NULL) +
      theme_map() +
      theme(panel.grid.major.x = element_line(colour = "#e6e5e0"),
            legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 2))
    save_png(pE, file.path(sdir, paste0(sid, "_E_nbhd_composition.png")), 9, 3.6)
  }

  invisible(sdir)
}

for (k in seq_len(nrow(slides))) {
  tryCatch(visualize_slide(k), error = function(e)
    warning(slides$slide_id[k], ": visualisation failed: ", conditionMessage(e), call. = FALSE))
  gc(FALSE)
}

###############################################################################
## CROSS-SLIDE SUMMARIES
###############################################################################

slide_levels <- slides |> arrange(group, slide_id) |> pull(slide_id)
grp_of <- setNames(as.character(slides$group), slides$slide_id)

comp_df <- cells |>
  count(slide_id, pop) |>
  group_by(slide_id) |>
  mutate(frac = n / sum(n)) |>
  ungroup() |>
  mutate(group = grp_of[slide_id],
         slide_id = factor(slide_id, levels = rev(slide_levels)),
         pop = factor(pop, levels = rev(POP_LEVELS)))

pS1 <- ggplot(comp_df |> filter(pop != "Other"), aes(frac, slide_id, fill = pop)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = POP_COLS, breaks = POP_LEVELS[1:4], labels = POP_LABEL[POP_LEVELS[1:4]]) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = "Called populations as a share of all segmented cells",
       subtitle = "remainder of each bar is 'Other' (unclassified) cells", x = NULL, y = NULL) +
  theme_map() + theme(panel.grid.major.x = element_line(colour = "#e6e5e0")) +
  guides(fill = guide_legend(nrow = 2))
save_png(pS1, file.path(PLOTDIR, "05_summary_population_share.png"), 9, 5)

ta_file <- file.path(OUT, "tissue_area_by_slide.csv")
if (file.exists(ta_file)) {
  ta <- read.csv(ta_file, stringsAsFactors = FALSE)
  dens <- cells |>
    group_by(slide_id) |>
    summarise(ASPC = sum(ASPC), EPI_CA = sum(EPI_CA), NE_CA = sum(NE_CA), .groups = "drop") |>
    left_join(ta, by = "slide_id") |>
    tidyr::pivot_longer(c(ASPC, EPI_CA, NE_CA), names_to = "pop", values_to = "n") |>
    mutate(per_mm2 = n / tissue_mm2, group = grp_of[slide_id],
           pop = factor(pop, levels = c("ASPC", "EPI_CA", "NE_CA")))
  write.csv(dens, file.path(OUT, "population_density_per_mm2_tissue.csv"), row.names = FALSE)
  pS2 <- ggplot(dens, aes(group, per_mm2)) +
    geom_point(aes(colour = pop), size = 3) +
    geom_text(aes(label = slide_id), hjust = -0.2, size = 2.8, colour = "#52514e") +
    facet_wrap(~ pop, scales = "free_y") +
    scale_colour_manual(values = POP_COLS, guide = "none") +
    labs(title = "Cell density per mm² of tissue, by slide",
         subtitle = "n = 3 slides per group: descriptive only", x = NULL, y = "cells / mm²") +
    theme_map() + theme(panel.grid.major.y = element_line(colour = "#e6e5e0"))
  save_png(pS2, file.path(PLOTDIR, "05_summary_density.png"), 10, 4.5)
}

if (!is.null(nbr_comp) && nrow(nbr_comp)) {
  ne <- nbr_comp |>
    mutate(group = grp_of[slide_id],
           slide_id = factor(slide_id, levels = slide_levels),
           population = factor(population, levels = POP_LEVELS))
  pS3 <- ggplot(ne, aes(slide_id, log2_enrich_vs_slide, colour = population)) +
    geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    scale_colour_manual(values = POP_COLS, labels = POP_LABEL) +
    labs(title = paste0("Composition of the ASPC neighbourhood (", RADIUS_UM, " µm) vs whole slide"),
         subtitle = "log2(share in neighbourhood / share on slide); > 0 = enriched near ASPC",
         x = NULL, y = "log2 enrichment") +
    theme_map() + theme(panel.grid.major.y = element_line(colour = "#e6e5e0")) +
    guides(colour = guide_legend(nrow = 2))
  save_png(pS3, file.path(PLOTDIR, "05_summary_nbhd_enrichment.png"), 10, 5)
}

## One multi-page PDF of every per-slide PNG + summaries, for sharing.
pngs <- c(
  file.path(PLOTDIR, c("05_summary_population_share.png", "05_summary_density.png",
                       "05_summary_nbhd_enrichment.png")),
  unlist(lapply(slide_levels, function(sid) {
    sort(list.files(file.path(VISDIR, sid), pattern = "\\.png$", full.names = TRUE))
  }))
)
pngs <- pngs[file.exists(pngs)]
if (length(pngs)) {
  pdf(file.path(PLOTDIR, "05_all_slides.pdf"), width = 11, height = 11, onefile = TRUE)
  for (f in pngs) {
    im <- readImage(f)
    grid.newpage()
    ar <- dim(im)[2] / dim(im)[1]
    grid.raster(as.raster(aperm(as.array(im)[, , 1:3], c(2, 1, 3))),
                width = if (ar <= 1) unit(1, "npc") else unit(1 / ar, "npc"),
                height = if (ar <= 1) unit(ar, "npc") else unit(1, "npc"))
  }
  dev.off()
  message("\nwrote 05_all_slides.pdf (", length(pngs), " pages)")
}

message("\n05_visualize.R complete -> ", VISDIR)
