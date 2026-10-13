#' Read one channel, optionally a tile. Returns a numeric matrix [w x h].
#' Coordinates are 0-based; x = column, y = row.
#'
#' Path-agnostic: the `tiff` argument is supplied by the caller from
#' `slides$tiff` (see 00_config.R). No slide-specific logic lives here --
#' PMOT32 now uses the same 8-cycle protocol as every other slide, so
#' channel indices are resolved upstream via `cyc_AR_DLL3` / `cyc_CK_PDGFRA`.
img_read <- function(tiff, c_index, x0 = 0, y0 = 0, w = NULL, h = NULL) {
  if (READER_ACTIVE == "bioformats") {
    sub <- list(C = c_index)
    if (!is.null(w) && !is.null(h))
      sub <- c(sub, list(x = as.integer(x0), y = as.integer(y0),
                         w = as.integer(w),  h = as.integer(h)))
    a <- tryCatch(
      as.array(RBioFormats::read.image(tiff, series = 1, subset = sub,
                                       normalize = FALSE, filter.metadata = TRUE)),
      error = function(e) {
        # fallback: whole plane then crop (older RBioFormats without x/y/w/h)
        z <- as.array(RBioFormats::read.image(tiff, series = 1,
                                              subset = list(C = c_index),
                                              normalize = FALSE,
                                              filter.metadata = TRUE))
        if (length(dim(z)) > 2) z <- z[, , 1]
        if (!is.null(w) && !is.null(h))
          z <- z[x0 + seq_len(min(w, nrow(z) - x0)),
                 y0 + seq_len(min(h, ncol(z) - y0)), drop = FALSE]
        z
      })
    if (length(dim(a)) > 2) a <- a[, , 1]
  } else {
    arr <- .tf$imread(tiff, series = 0L, key = as.integer(c_index - 1L))
    a <- t(reticulate::py_to_r(.np$asarray(arr)))
    if (!is.null(w) && !is.null(h))
      a <- a[x0 + seq_len(min(w, nrow(a) - x0)),
             y0 + seq_len(min(h, ncol(a) - y0)), drop = FALSE]
  }
  storage.mode(a) <- "double"
  a
}
