#' Align CDR3s for a logo
#'
#' Trims `trim` residues at each end, then aligns every sequence to the modal
#' length `L` (ties: the shorter), TCRdist-style: longer sequences lose
#' residues from the middle (keep the first `ceiling(L / 2)` and the last
#' `floor(L / 2)`), shorter ones get gaps (`"-"`) inserted in the middle,
#' after their first `ceiling(len / 2)` residues.
#'
#' @param cdr3 Character vector of CDR3s.
#' @param trim Residues trimmed at each end.
#' @return Character vector of aligned sequences, all of the modal length.
#' @export
#' @examples
#' logo_align(c("CASSLGQAYEQYF", "CASSLGQAYEQYF", "CASSLGAYEQYF", "CASSLGQGAYEQYF"))
logo_align <- function(cdr3, trim = 1) {
  s <- cdr3[!is.na(cdr3)]
  if (!length(s)) return(character())
  n <- nchar(s)
  s <- substr(s, trim + 1, n - trim)
  len <- nchar(s)
  tab <- table(len)
  L <- as.integer(names(tab)[tab == max(tab)][1])
  if (L <= 0) return(rep("", length(s)))
  vapply(seq_along(s), function(i) {
    x <- s[i]
    l <- len[i]
    if (l == L) return(x)
    if (l > L) {
      h <- ceiling(L / 2)
      return(paste0(substr(x, 1, h), substr(x, l - floor(L / 2) + 1, l)))
    }
    h <- ceiling(l / 2)
    paste0(substr(x, 1, h), strrep("-", L - l), substr(x, h + 1, l))
  }, "")
}

#' Column frequencies of aligned sequences
#'
#' @param aligned Output of [logo_align()].
#' @return A tibble `pos`, `letter`, `freq` (count / all sequences; gaps are
#'   not letters, so gappy columns are shorter), sorted within a column from
#'   least to most frequent (ties alphabetical), i.e. bottom-to-top stacking.
#' @export
#' @examples
#' logo_freqs(logo_align(c("CASSLGQAYEQYF", "CASSLGAYEQYF")))
logo_freqs <- function(aligned) {
  if (!length(aligned) || all(nchar(aligned) == 0)) {
    return(tibble::tibble(pos = integer(), letter = character(), freq = numeric()))
  }
  L <- nchar(aligned[1])
  m <- do.call(rbind, strsplit(aligned, ""))
  dim(m) <- c(length(aligned), L)
  out <- lapply(seq_len(L), function(j) {
    t <- table(m[, j])
    t <- t[names(t) != "-"]
    if (!length(t)) return(NULL)
    f <- as.numeric(t) / length(aligned)
    o <- order(f, -rank(names(t)))  # least frequent first; ties: Z..A so A ends on top
    tibble::tibble(pos = j, letter = names(t)[o], freq = f[o])
  })
  dplyr::bind_rows(out)
}

# Glyph outlines: built-in table (DejaVu Sans Mono Bold) or, if requested and
# available, systemfonts::glyph_outline() for another font.
.glyphs <- function(font = NULL) {
  key <- paste0("glyphs_", if (is.null(font)) "builtin" else font)
  if (!is.null(.tcrhex_cache[[key]])) return(.tcrhex_cache[[key]])
  g <- NULL
  if (!is.null(font) && .has_glyph_outline()) {
    g <- tryCatch(.glyphs_systemfonts(font), error = function(e) NULL)
    if (is.null(g)) {
      warning("Could not outline font '", font, "'; using the built-in DejaVu Sans ",
              "Mono Bold outlines.", call. = FALSE)
    }
  }
  if (is.null(g)) {
    f <- system.file("extdata", "glyphs_dejavu_sans_mono_bold.tsv", package = "tcrhex")
    g <- readr::read_tsv(f, col_types = "ciidd", progress = FALSE, na = character())
  }
  # normalise every letter to its ink box [0, 1] x [0, 1]
  g <- dplyr::mutate(dplyr::group_by(g, .data$letter),
                     x = (.data$x - min(.data$x)) / (max(.data$x) - min(.data$x)),
                     y = (.data$y - min(.data$y)) / (max(.data$y) - min(.data$y)))
  g <- dplyr::ungroup(g)
  .tcrhex_cache[[key]] <- g
  g
}

.has_glyph_outline <- function() {
  requireNamespace("systemfonts", quietly = TRUE) &&
    "glyph_outline" %in% getNamespaceExports("systemfonts")
}

.glyphs_systemfonts <- function(font) {
  sf <- asNamespace("systemfonts")
  letters20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  info <- sf$match_fonts(font, weight = "bold")
  gi <- if ("weight" %in% names(formals(sf$glyph_info))) {
    sf$glyph_info(letters20, family = font, weight = "bold")
  } else {
    sf$glyph_info(letters20, family = font, bold = TRUE)
  }
  out <- sf$glyph_outline(gi$index, info$path, info$index, size = 100)
  # glyph_outline returns glyph, contour, x, y
  tibble::tibble(letter = letters20[out$glyph], contour = as.integer(out$contour),
                 order = stats::ave(out$x, out$glyph, out$contour, FUN = seq_along),
                 x = out$x, y = out$y)
}

# Logo polygons for one bin and one chain.
# freqs: logo_freqs(); edge from p0 to p1 (data units); baseline `offset` off
# the edge (outwards = left normal of p0->p1), full height `height`.
.logo_polys <- function(freqs, p0, p1, offset = 0.35, height = 2.6, pad = 0.04,
                        glyphs = .glyphs()) {
  if (!nrow(freqs)) return(NULL)
  L <- max(freqs$pos)
  dx <- p1[1] - p0[1]
  dy <- p1[2] - p0[2]
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  nx <- -uy
  ny <- ux
  w <- len / L
  freqs <- dplyr::mutate(dplyr::group_by(freqs, .data$pos),
                         top = cumsum(.data$freq) * height,
                         bottom = .data$top - .data$freq * height)
  freqs <- dplyr::ungroup(freqs)
  freqs$glyph_id <- seq_len(nrow(freqs))
  d <- dplyr::inner_join(freqs, glyphs, by = "letter", relationship = "many-to-many")
  u <- (d$pos - 1 + pad + d$x * (1 - 2 * pad)) * w
  v <- offset + d$bottom + d$y * (d$top - d$bottom)
  tibble::tibble(glyph = d$glyph_id, contour = d$contour, order = d$order,
                 letter = d$letter, x = p0[1] + u * ux + v * nx,
                 y = p0[2] + u * uy + v * ny)
}
