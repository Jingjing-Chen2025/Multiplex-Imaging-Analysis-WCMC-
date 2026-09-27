###############################################################################
# CONFIG + BATCH TABLE -- ASPC / NEPC multiplex study
# Run this first. It sets global parameters, installs/checks dependencies,
# builds the slide/batch table, and reports confounding structure.
###############################################################################

## --- Java heap for BioFormats: MUST be set before rJava/RBioFormats loads ----
options(java.parameters = "-Xmx24g")

## ---------------------------------------------------------------- PATHS -----
## Both can be overridden before source(); OUT is absolute so every script
## writes to the same place regardless of the R working directory.
if (!exists("ROOT")) ROOT <- "/Volumes/T9/Case Studies"
if (!exists("OUT"))  OUT  <- "/Volumes/T9/results_ASPC_NEPC"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------ PARAMETERS ----
## Pixel size, tile size and analysis radii are NOT set here: pixel size is
## read from OME metadata in 01 (0.28 um/px on this cohort), tile size is owned
## by 02, and spatial radii by 04. Setting them here previously leaked stale
## values (PX_UM 0.23, TILE_SIZE 4096) into later scripts sourced in-session.

## Antibody caveat: the protocol stains PDGFa (LIGAND, CST 5241T),
## NOT PDGFRA (receptor). ASPC definition uses the channel actually stained.
PDGF_IS_LIGAND <- TRUE

## Reader backend:
##   "auto"       -> try RBioFormats (Java); fall back to tifffile (Python)
##   "bioformats" -> force RBioFormats
##   "tifffile"   -> force Python tifffile (no Java needed)
READER <- "tifffile"   # rJava/RBioFormats are not installed on this machine
PY_ENV <- "r-reticulate"

## Install anything missing without prompting?
AUTO_INSTALL <- TRUE

## -------------------------------------------------- DEPENDENCY BOOTSTRAP ----
.ensure <- function(pkgs, bioc = FALSE) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(need)) return(invisible(TRUE))
  if (!AUTO_INSTALL) {
    message("Missing: ", paste(need, collapse = ", ")); return(invisible(FALSE))
  }
  message("Installing: ", paste(need, collapse = ", "))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    try(BiocManager::install(need, ask = FALSE, update = FALSE), silent = TRUE)
  } else {
    try(install.packages(need, repos = "https://cloud.r-project.org"), silent = TRUE)
  }
  invisible(all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)))
}

## core (always required)
.ensure(c("xml2","dplyr","tidyr","purrr","tibble","ggplot2","MASS","RANN"))
.ensure("EBImage", bioc = TRUE)

.core_ok <- all(vapply(c("xml2","dplyr","tidyr","purrr","tibble",
                         "ggplot2","MASS","RANN","EBImage"),
                       requireNamespace, logical(1), quietly = TRUE))
if (!.core_ok)
  stop("Core packages failed to install. Install manually:\n",
       '  install.packages(c("xml2","dplyr","tidyr","purrr","tibble","ggplot2","MASS","RANN"))\n',
       '  BiocManager::install("EBImage")', call. = FALSE)

## ------------------------------------------------------ READER SELECTION ----
.try_bioformats <- function() {
  .ensure("rJava")
  if (!requireNamespace("rJava", quietly = TRUE)) return(FALSE)
  ok <- tryCatch({ rJava::.jinit(); TRUE }, error = function(e) FALSE)
  if (!ok) return(FALSE)
  .ensure("RBioFormats", bioc = TRUE)
  requireNamespace("RBioFormats", quietly = TRUE)
}

.try_tifffile <- function() {
  .ensure("reticulate")
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  try(reticulate::use_virtualenv(PY_ENV, required = TRUE), silent = TRUE)
  ok <- tryCatch({ reticulate::import("tifffile"); TRUE }, error = function(e) FALSE)
  if (!ok) ok <- tryCatch({
    reticulate::py_install(c("tifffile","numpy","imagecodecs"))
    reticulate::import("tifffile"); TRUE
  }, error = function(e) FALSE)
  ok
}

READER_ACTIVE <- switch(
  READER,
  "bioformats" = if (.try_bioformats()) "bioformats"
                 else stop("RBioFormats/Java unavailable. Set READER <- \"tifffile\".",
                           call. = FALSE),
  "tifffile"   = if (.try_tifffile()) "tifffile"
                 else stop("Python tifffile unavailable.", call. = FALSE),
  "auto"       = if (.try_bioformats()) "bioformats"
                 else if (.try_tifffile()) { 
                   message("rJava/RBioFormats unavailable -- using Python tifffile."); "tifffile"
                 } else stop(
                   "No usable reader.\n",
                   "Java route : brew install --cask temurin && sudo R CMD javareconf\n",
                   "             then install.packages('rJava'); BiocManager::install('RBioFormats')\n",
                   "Python route: install.packages('reticulate'); ",
                   "reticulate::py_install(c('tifffile','numpy','imagecodecs'))",
                   call. = FALSE),
  stop("READER must be 'auto', 'bioformats', or 'tifffile'.", call. = FALSE)
)

suppressPackageStartupMessages({
  library(EBImage); library(xml2)
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2)
  if (READER_ACTIVE == "bioformats") { library(rJava); library(RBioFormats) }
  else                                library(reticulate)
})
message("Reader backend active: ", READER_ACTIVE)

## ------------------------------------------------------------ BATCH TABLE ---
## Derived from the COMET protocol logs.
##   All slides: 8-cycle protocol, negative controls (AR/DLL3 = Cycle 6)
##   Negative-control pairing (species/fluorophore matched):
##     Cycle 5 (anti-Mouse647 / anti-Rabbit555) -> controls Cycle 6 (AR, DLL3)
##     Cycle 3 (anti-Mouse555 / anti-Rabbit647) -> controls Cycle 7 (panCK, PDGFa)
slides <- tribble(
  ~slide_id,  ~group,           ~stainer, ~run,       ~comet_sw, ~n_cycles, ~has_negctrl, ~runtime_min, ~dir, ~file,

  "PM0T32",   "NEPC",           1L, "20260826", "1.5.0.0", 8L, TRUE, NA_integer_,
    "20260826_155725_1_ybucNp_JJ_full panel_PM0 T32",
    "20260826_155725_1_ybucNp_JJ_full panel_PM0 T32.tiff",

  "PM677_A8", "NEPC",           1L, "20260608", "1.4.1.0", 8L, TRUE, 760L,
    "20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8",
    "20260608_164617_1_ahXvmt_JJ_full panel_PM677 A8.tiff",

  "PM12_B3",  "NEPC",           2L, "20260608", "1.4.1.0", 8L, TRUE, 775L,
    "20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3",
    "20260608_164617_2_9ZdZ8V_JJ_full panel_PM12 B3.tiff",

  "PM1024_B", "ADENO",          4L, "20260608", "1.4.1.0", 8L, TRUE, 799L,
    "20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B",
    "20260608_164617_4_CA2z2F_JJ_full panel_PM1024 B.tiff",

  "PM423",    "ADENO",          2L, "20260625", "1.5.0.0", 8L, TRUE, 500L,
    "20260625_131559_2_ZxKUqB_JJ_full panel_PM423",
    "20260625_131559_2_ZxKUqB_JJ_full panel_PM423.tiff",

  "PM159T5",  "ADENO",          3L, "20260625", "1.5.0.0", 8L, TRUE, 529L,
    "20260625_131600_3_VHgSnt_JJ_full panel_PM159T5",
    "20260625_131600_3_VHgSnt_JJ_full panel_PM159T5.tiff"
) |>
  mutate(
    tiff  = file.path(ROOT, dir, file),
    xml   = file.path(ROOT, dir, "jobfile.xml"),
    group = factor(group, levels = c("ADENO", "NEPC")),

    ## protocol-dependent cycle indices
    cyc_AR_DLL3  = if_else(n_cycles == 6L, 4L, 6L),
    cyc_CK_PDGF  = if_else(n_cycles == 6L, 5L, 7L),
    ctrl_AR_DLL3 = if_else(has_negctrl, 5L, NA_integer_),
    ctrl_CK_PDGF = if_else(has_negctrl, 3L, NA_integer_),

    ## analysis tiers
    tier1 = TRUE,                 # all 6 slides, AF-only correction  (PRIMARY)
    tier2 = has_negctrl,          # 6 slides, negative-control correction
    tier3 = run == "20260608"     # within-batch: 2 NEPC vs 1 adeno, same sw/lots
  )

## ------------------------------------------------------- FILE EXISTENCE -----
.bad <- slides |>
  mutate(tiff_ok = file.exists(tiff), xml_ok = file.exists(xml)) |>
  filter(!tiff_ok | !xml_ok)
if (nrow(.bad)) {
  print(as.data.frame(.bad[, c("slide_id","tiff_ok","xml_ok","tiff")]))
  stop("Missing file(s) above -- check that the T9 volume is mounted.", call. = FALSE)
}

## -------------------------------------------------- CONFOUNDING QC ---------
cat("\n--- group x run ---\n");             print(with(slides, table(group, run)))
cat("\n--- group x COMET software ---\n");  print(with(slides, table(group, comet_sw)))
cat("\n--- group x negative control ---\n");print(with(slides, table(group, has_negctrl)))
cat("\n--- group x stainer ---\n");         print(with(slides, table(group, stainer, useNA = "ifany")))
cat("\n--- group x n_cycles ---\n");        print(with(slides, table(group, n_cycles)))

.sw_nepc  <- unique(slides$comet_sw[slides$group == "NEPC"])
.sw_adeno <- unique(slides$comet_sw[slides$group == "ADENO"])
if (!setequal(.sw_nepc, .sw_adeno))
  warning("COMET software version differs between groups. Rely on Tier 3 ",
          "(20260608 within-batch) as the technical control.", call. = FALSE)

if (any(!slides$has_negctrl))
  warning("One or more slides lack negative-control cycles. Tier 1 uses AF-only ",
          "correction on ALL slides for comparability; Tier 2 drops them.",
          call. = FALSE)

if (PDGF_IS_LIGAND)
  message("\nNOTE: panel stains PDGFa (ligand, CST 5241T), not PDGFRA (receptor). ",
          "ASPC is called as PDGFa+ (03_classify.R); CEBPg splits ASPC into ",
          "CEBPg+ / CEBPg- subsets and panCK is not used to exclude ASPC.")

## ------------------------------------------------------------ SAVE ---------
attr(slides, "reader") <- READER_ACTIVE
saveRDS(slides, file.path(OUT, "slides_final.rds"))

cat("\n--- slide table ---\n")
print(as.data.frame(slides[, c("slide_id","group","stainer","run","comet_sw",
                               "n_cycles","has_negctrl","cyc_AR_DLL3",
                               "ctrl_AR_DLL3","ctrl_CK_PDGF","tier2","tier3")]))
cat("\nTier 1 (primary, AF-only): ", sum(slides$tier1), " slides\n",
    "Tier 2 (neg-control corr): ", sum(slides$tier2), " slides\n",
    "Tier 3 (within-batch)    : ", sum(slides$tier3), " slides\n", sep = "")
message("\n00_config.R complete -> ", file.path(OUT, "slides_final.rds"))
