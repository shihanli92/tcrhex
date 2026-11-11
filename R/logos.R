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
    g <- readr::read_tsv(f, col_types = "ciiddd", progress = FALSE, na = character())
    g$advance <- NULL
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

# Glyph outlines in font units for words: x / advance width (so a character
# fills its monospace cell with its own side bearings), y / cap height (so
# "-" stays small and mid-height); "|" is scaled to the cap height.
.word_glyphs <- function() {
  if (!is.null(.tcrhex_cache$word_glyphs)) return(.tcrhex_cache$word_glyphs)
  f <- system.file("extdata", "glyphs_dejavu_sans_mono_bold.tsv", package = "tcrhex")
  g <- readr::read_tsv(f, col_types = "ciiddd", progress = FALSE, na = character())
  cap <- max(g$y[g$letter == "H"])
  g$x <- g$x / g$advance
  g$y <- g$y / cap
  bar <- g$letter == "|"
  g$y[bar] <- (g$y[bar] - min(g$y[bar])) / (max(g$y[bar]) - min(g$y[bar]))
  g$advance <- NULL
  .tcrhex_cache$word_glyphs <- g
  g
}

#' Gene-usage words for a logo
#'
#' Turns the V (or J) genes of a bin's clones into the words of a gene logo
#' (`"TRAV12-2"` -> `"AV12-2"`) and returns the most frequent ones. V and J
#' are counted independently, so each gene appears once.
#'
#' @param genes V or J genes of the clones in a bin.
#' @param min_freq Keep genes found in at least this fraction of the clones.
#' @param max_words Keep at most this many genes (most frequent first, ties
#'   alphabetical).
#' @return A tibble `word`, `freq` (fraction of all clones, so missing genes
#'   make the stack shorter), most frequent first.
#' @export
#' @examples
#' gene_words(c("TRAV12-2", "TRAV12-2", "TRAV21"))
#' gene_words(c("TRAJ42", "TRAJ42", NA))
gene_words <- function(genes, min_freq = 0.1, max_words = 4) {
  n <- length(genes)
  w <- sub("^TR", "", genes[!is.na(genes)])
  if (!length(w)) return(tibble::tibble(word = character(), freq = numeric()))
  tab <- table(w)
  out <- tibble::tibble(word = names(tab), freq = as.numeric(tab) / n)
  out <- out[order(-out$freq, out$word), ]
  out <- out[out$freq >= min_freq, ]
  utils::head(out, max_words)
}

# Word logo along an edge p0 -> p1, hanging on the outer side (right of the
# direction), letters upright (their up = the edge's left normal, i.e. into
# the hexagon for the lower edges). Most frequent word nearest the edge; each
# word spans the part [from, to] of the edge (fractions), height = freq * height.
.word_logo_polys <- function(words, p0, p1, offset = 0.35, height = 2.6, pad = 0.04,
                             glyphs = .word_glyphs(), from = 0, to = 1) {
  if (is.null(words) || !nrow(words)) return(NULL)
  dx <- p1[1] - p0[1]
  dy <- p1[2] - p0[2]
  len <- sqrt(dx^2 + dy^2)
  ux <- dx / len
  uy <- dy / len
  nx <- -uy
  ny <- ux
  out <- vector("list", nrow(words))
  vt <- -offset
  gid <- 0
  for (i in seq_len(nrow(words))) {
    h <- words$freq[i] * height
    vb <- vt - h
    ch <- strsplit(words$word[i], "")[[1]]
    span <- len * (to - from)
    cw <- span * (1 - 2 * pad) / length(ch)
    parts <- vector("list", length(ch))
    for (k in seq_along(ch)) {
      gg <- glyphs[glyphs$letter == ch[k], ]
      if (!nrow(gg)) next
      gid <- gid + 1
      u <- len * from + span * pad + (k - 1) * cw + gg$x * cw
      v <- vb + gg$y * h
      parts[[k]] <- tibble::tibble(glyph = gid, contour = gg$contour, order = gg$order,
                                   letter = ch[k], x = p0[1] + u * ux + v * nx,
                                   y = p0[2] + u * uy + v * ny)
    }
    out[[i]] <- dplyr::bind_rows(parts)
    vt <- vb
  }
  dplyr::bind_rows(out)
}
