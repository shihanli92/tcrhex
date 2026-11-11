#' tcrhex: paired TCRdist hex rows from single-cell TCR data
#'
#' Draws "paired TCRdist hex rows" with ggplot2: one row per query TCR, one
#' hexagonal frame per paired-TCRdist bin, one small hexagon per clone, coloured
#' by any categorical column. See `vignette("tcrhex")`.
#'
#' The V-gene loop table in `inst/extdata/tcrdist_vloops_human.tsv` is derived
#' from tcrdist3's `alphabeta_gammadelta_db.tsv` (MIT License, Copyright (c) 2020
#' Koshlan Mayer-Blackwell); the distance kernel ports pwseqdist's
#' `nb_vector_tcrdist`. The logo glyph outlines in
#' `inst/extdata/glyphs_dejavu_sans_mono_bold.tsv` are traced from DejaVu Sans
#' Mono Bold (Bitstream Vera / DejaVu license). See `inst/COPYRIGHTS`.
#'
#' @keywords internal
#' @useDynLib tcrhex, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data .env := enquo quo_is_null eval_tidy as_label
#' @importFrom stats setNames
"_PACKAGE"

utils::globalVariables(c(".clone_id", ".row", ".bin", ".level", ".pattern", "donor"))
