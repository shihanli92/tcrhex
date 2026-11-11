#' Honeycomb positions for k small hexagons
#'
#' All centres of a pointy-top hex lattice (`x = sqrt(3) * (q + r / 2)`,
#' `y = -1.5 * r`, circumradius 1) are ranked by distance to the origin, ties
#' broken by angle (`distance + 1e-6 * atan2(y, x)`); the `k` closest are kept
#' and ordered row by row, top to bottom, then left to right. Filling them in
#' this order with clones sorted by fill level gives horizontal bands.
#'
#' @param k Number of hexagons.
#' @return A tibble with `slot` (1..k), `x`, `y`.
#' @export
#' @examples
#' hex_layout(7)
hex_layout <- function(k) {
  if (k <= 0) return(tibble::tibble(slot = integer(), x = numeric(), y = numeric()))
  R <- ceiling(sqrt(k)) + 2
  g <- expand.grid(q = -R:R, r = -R:R)
  x <- sqrt(3) * (g$q + g$r / 2)
  y <- -1.5 * g$r
  key <- sqrt(x^2 + y^2) + 1e-6 * atan2(y, x)
  o <- order(key)[seq_len(k)]
  x <- x[o]
  y <- y[o]
  yr <- round(y, 6)
  o2 <- order(-yr, round(x, 6))
  tibble::tibble(slot = seq_len(k), x = x[o2], y = y[o2])
}

# Radius of the patch: largest centre distance from the origin.
.patch_extent <- function(k) {
  if (k <= 0) return(0)
  l <- hex_layout(k)
  max(sqrt(l$x^2 + l$y^2))
}

#' Vertices of pointy-top hexagons
#'
#' @param x,y Centres.
#' @param r Circumradius (scalar or one per hexagon).
#' @return A tibble with `id` (index of the hexagon), `x`, `y` (6 rows per
#'   hexagon, starting at the apex, clockwise).
#' @export
#' @examples
#' hex_vertices(0, 0, 1)
hex_vertices <- function(x, y, r = 1) {
  n <- length(x)
  r <- rep_len(r, n)
  ang <- pi / 2 - (0:5) * pi / 3
  tibble::tibble(id = rep(seq_len(n), each = 6),
                 x = rep(x, each = 6) + rep(r, each = 6) * rep(cos(ang), n),
                 y = rep(y, each = 6) + rep(r, each = 6) * rep(sin(ang), n))
}

#' Split k slots over levels by the largest-remainder method
#'
#' Each level gets `floor(k * n / sum(n))` slots; the slots left over go to
#' the levels with the largest remainders (ties: larger count first, then
#' level order). If `sum(n) <= k`, every level keeps all its members.
#'
#' @param n Integer counts per level (in level order).
#' @param k Number of slots.
#' @return Integer vector of slots per level (sums to `min(k, sum(n))`).
#' @export
#' @examples
#' allot(c(a = 30, b = 10, c = 2), 19)
allot <- function(n, k) {
  tot <- sum(n)
  if (tot <= k) return(stats::setNames(as.integer(n), names(n)))
  quota <- k * n / tot
  base <- floor(quota + 1e-9)
  rem <- quota - base
  left <- k - sum(base)
  if (left > 0) {
    o <- order(-round(rem, 9), -n, seq_along(n))
    base[o[seq_len(left)]] <- base[o[seq_len(left)]] + 1
  }
  stats::setNames(as.integer(base), names(n))
}

#' Paired TCRdist bins
#'
#' Bin 0 holds `d == 0` (the query itself); bin `i >= 1` holds
#' `(i - 1) * width < d <= i * width`, the last bin capped at the radius.
#'
#' @param d Distances.
#' @param width Bin width.
#' @param radius Search radius.
#' @return For `tcrdist_bin()`, integer bins (`NA` beyond the radius); for
#'   `bin_labels()`, a tibble of `bin`, `lo`, `hi`, `label`.
#' @export
#' @examples
#' tcrdist_bin(c(0, 1, 12, 13, 108), 12, 108)
#' bin_labels(12, 108)
tcrdist_bin <- function(d, width = 12, radius = 108) {
  b <- ifelse(d <= 0, 0L, as.integer(ceiling(d / width)))
  b[d > radius] <- NA_integer_
  b
}

#' @rdname tcrdist_bin
#' @export
bin_labels <- function(width = 12, radius = 108) {
  nb <- ceiling(radius / width)
  i <- seq_len(nb)
  lo <- (i - 1) * width + 1
  hi <- pmin(i * width, radius)
  lab <- ifelse(lo == hi, as.character(lo), paste0(lo, "\u2013", hi))
  tibble::tibble(bin = c(0L, i), lo = c(0, lo), hi = c(0, hi), label = c("0", lab))
}
