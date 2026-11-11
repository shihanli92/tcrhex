#' Paired TCRdist hex rows
#'
#' One row per query TCR; along the row, a hexagonal frame per paired-TCRdist
#' bin (0, 1-12, 13-24, ... up to `radius`); inside each frame one small
#' hexagon per clone (donor x TRAV + CDR3a x TRBV + CDR3b), up to `max_hex`,
#' coloured by any categorical column. Returns a ggplot object, so
#' `+ scale_fill_manual()`, `+ labs()`, `+ theme()` work as usual. Save it
#' with [ggsave_hexrows()], which sizes the page.
#'
#' Fill: any categorical expression, the ggplot way: `aes(fill = disease)`,
#' `aes(fill = paste(tissue, timepoint))`, `aes(fill = tissue_class(tissue))`.
#' Values are computed per cell and summarised per clone with
#' `summarise_fill` (one value stays, several give `"multiple"`). Hexagons
#' in a bin are ordered by the fill's level order (factor levels, else most
#' common first), giving horizontal bands. Without a fill scale, the default
#' palette is used: up to 8 fixed hues, `"Other"` (#b4b3ac), `"unknown"` / NA
#' (#d6d5ce), with clone counts in the legend.
#'
#' When a bin holds more than `max_hex` clones, the slots are split over the
#' fill levels by the largest-remainder method (so the shown hexagons keep the
#' bin's composition); within a level the closest clones (smallest d, then
#' most cells) are shown. The number below-right of a frame counts all
#' clones of the bin (darker when some are hidden).
#'
#' @inheritParams hexrows_data
#' @param title,subtitle Plot title and subtitle (`NULL` for the defaults;
#'   change them later with `+ labs()`).
#' @return A ggplot object, with the [hexrows_data()]
#'   result attached as `attr(p, "hexrows")`.
#' @export
#' @examples
#' cells <- tcrhex_example_cells()   # synthetic data
#' hexrows(cells, ggplot2::aes(fill = disease), queries = top_public(cells, 3))
#' hexrows(cells, ggplot2::aes(fill = tissue_class(tissue)),
#'         queries = "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF",
#'         summarise_fill = summarise_tissue) + scale_fill_tissue()
hexrows <- function(cells, mapping = ggplot2::aes(), queries = NULL, radius = 108,
                    bin_width = 12, max_hex = 19, logos = TRUE, logo_trim = 1,
                    logo_font = NULL, summarise_fill = summarise_unique,
                    fold_levels = 8, dominant = 0.10, row_tag = NULL,
                    clones = NULL, hits = NULL, title = NULL, subtitle = NULL,
                    unit_in = 0.072, organism = NULL) {
  hd <- hexrows_data(cells, mapping, queries = queries, radius = radius,
                     bin_width = bin_width, max_hex = max_hex, logos = logos,
                     logo_trim = logo_trim, logo_font = logo_font,
                     summarise_fill = summarise_fill, fold_levels = fold_levels,
                     dominant = dominant, row_tag = row_tag, clones = clones,
                     hits = hits, unit_in = unit_in, organism = organism)
  hexrows_plot(hd, title = title, subtitle = subtitle)
}

#' Draw a hexrows_data object
#'
#' @param hd Result of [hexrows_data()].
#' @inheritParams hexrows
#' @return A ggplot object.
#' @export
hexrows_plot <- function(hd, title = NULL, subtitle = NULL) {
  m <- hd$meta
  env <- .default_scale_env(hd)
  hexes <- hd$hexes
  names(hexes)[names(hexes) == "level"] <- ".level"
  names(hexes)[names(hexes) == "pattern"] <- ".pattern"
  has_fill <- !is.null(m$fill_label)
  has_pat <- !is.null(m$pattern_label)
  hex_aes <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$group)
  if (has_fill) hex_aes$fill <- rlang::quo(.data$.level)
  if (has_pat) hex_aes$pattern <- rlang::quo(.data$.pattern)
  p <- ggplot2::ggplot(environment = env) +
    ggplot2::geom_segment(data = hd$lines, ggplot2::aes(x = .data$x, y = .data$y,
                          xend = .data$xend, yend = .data$yend),
                          colour = .hx$row_line, linewidth = .pt2lw(0.6)) +
    ggplot2::geom_polygon(data = hd$frames, ggplot2::aes(x = .data$x, y = .data$y,
                          group = .data$group), fill = .hx$frame_fill,
                          colour = .hx$frame_line, linewidth = .pt2lw(.hx$frame_lwd_pt))
  hex_layer <- if (has_fill) {
    geom_hextile(data = hexes, mapping = hex_aes, colour = .hx$frame_fill,
                 linewidth = .pt2lw(.hx$hex_lwd_pt))
  } else {
    geom_hextile(data = hexes, mapping = hex_aes, colour = .hx$frame_fill,
                 fill = .hx$hues[1], linewidth = .pt2lw(.hx$hex_lwd_pt))
  }
  p <- p + hex_layer
  if (!is.null(hd$logos) && nrow(hd$logos)) {
    p <- p + ggplot2::geom_polygon(data = hd$logos, ggplot2::aes(x = .data$x, y = .data$y,
                                   group = .data$group, subgroup = .data$subgroup),
                                   fill = .hx$ink, colour = NA)
  }
  lab <- hd$labels
  p <- p +
    ggplot2::geom_text(data = hd$counts, ggplot2::aes(x = .data$x, y = .data$y,
                       label = .data$label, colour = I(.data$colour)),
                       size = .hx$size_count / ggplot2::.pt, hjust = 0, vjust = 1,
                       family = "mono") +
    ggplot2::geom_text(data = lab[lab$part %in% c("line1", "line2", "tag"), ],
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                    colour = I(.data$colour), size = I(.data$size / ggplot2::.pt)),
                       hjust = 0, vjust = 0.5, family = "mono") +
    ggplot2::geom_text(data = lab[lab$part == "right", ],
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                    colour = I(.data$colour), size = I(.data$size / ggplot2::.pt)),
                       hjust = 0, vjust = 0.5) +
    ggplot2::geom_text(data = lab[lab$part %in% c("bin", "axis"), ],
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                    colour = I(.data$colour), size = I(.data$size / ggplot2::.pt)),
                       hjust = 0.5, vjust = 0) +
    ggplot2::coord_fixed(xlim = m$xlim, ylim = m$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(title = title %||% .default_title(hd),
                  subtitle = subtitle %||% .default_subtitle(hd),
                  fill = m$fill_label, pattern = m$pattern_label) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, byrow = TRUE),
                    pattern = ggplot2::guide_legend(ncol = 3, byrow = TRUE)) +
    theme_hexrows()
  attr(p, "hexrows") <- hd
  p
}

#' Theme of the hexrows figure
#' @param base_size Base font size (pt).
#' @return A ggplot2 theme.
#' @export
theme_hexrows <- function(base_size = 8) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = .hx$page, colour = NA),
      panel.background = ggplot2::element_rect(fill = .hx$page, colour = NA),
      plot.title = ggplot2::element_text(colour = .hx$ink, face = "bold",
                                         size = base_size + 2, hjust = 0,
                                         margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(colour = .hx$ink2, size = base_size - 1.5,
                                            hjust = 0, lineheight = 1.15,
                                            margin = ggplot2::margin(b = 6)),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = ggplot2::element_text(colour = .hx$ink2, size = base_size - 1),
      legend.text = ggplot2::element_text(colour = .hx$ink, size = base_size - 1.5),
      legend.key.size = grid::unit(9, "pt"),
      legend.margin = ggplot2::margin(t = 2),
      plot.margin = ggplot2::margin(10, 12, 8, 12))
}

.default_title <- function(hd) {
  n <- hd$meta$n_rows
  if ("donors" %in% names(hd$rows)) {
    sprintf("The %d most public TCRs and their paired-TCRdist neighbourhoods", n)
  } else {
    sprintf("Paired-TCRdist neighbourhoods of %d query TCR%s", n, if (n > 1) "s" else "")
  }
}

.default_subtitle <- function(hd) {
  m <- hd$meta
  s1 <- sprintf(paste0(
    "One hexagon = one clone (donor \u00d7 V+CDR3\u03b1 \u00d7 V+CDR3\u03b2), at most %d per bin; ",
    "when a bin has more, the slots are split over the fill levels by largest remainder ",
    "and the closest clones are shown."), m$max_hex)
  s2 <- "The number below each frame counts all clones in the bin (darker: some not shown)."
  s3 <- if (isTRUE(m$logos)) sprintf(paste0(
    "Logos: CDR3\u03b1 (upper-left edge) and CDR3\u03b2 (upper-right edge) of the bin's clones, ",
    "%d residue%s trimmed at each end, aligned to the modal length."),
    m$logo_trim, if (m$logo_trim == 1) "" else "s") else NULL
  paste(c(s1, paste(c(s2, s3), collapse = " ")), collapse = "\n")
}

# Environment holding the default discrete fill / pattern scales for this
# plot; ggplot2 looks up scale_<aes>_discrete in the plot environment first,
# so any scale the user adds replaces these silently.
.default_scale_env <- function(hd) {
  env <- new.env(parent = globalenv())
  lv <- hd$levels
  counts <- stats::setNames(lv$clones, as.character(lv$level))
  env$scale_fill_discrete <- function(...) {
    scale_fill_hexrows(counts = counts)
  }
  pl <- hd$pattern_levels
  if (!is.null(pl)) {
    pc <- stats::setNames(pl$clones, as.character(pl$level))
    env$scale_pattern_discrete <- function(...) {
      vals <- rep_len(c("none", "stripe", "crosshatch", "circle"), length(pc))
      ggplot2::discrete_scale("pattern", palette = function(n) vals[seq_len(n)],
                              labels = .count_labels(pc), name = ggplot2::waiver())
    }
  }
  env
}

.count_labels <- function(counts) {
  function(breaks) {
    n <- counts[breaks]
    ifelse(is.na(n), breaks, sprintf("%s (%s)", breaks, .fmt(n)))
  }
}

#' Default hexrows fill scale
#'
#' Fixed hues in level order (#2a78d6 #eb6834 #1baf7a #eda100 #e87ba4 #008300
#' #6250d6 #e34948, then ggplot hues), `"Other"` #b4b3ac, `"unknown"` / NA
#' #d6d5ce. Used automatically by [hexrows()] when no fill scale is added.
#'
#' @param counts Optional named clone counts per level, in level order; the
#'   names fix the colours (so `"Other"` and `"unknown"` get their greys) and
#'   the legend reads `"level (n)"`.
#' @param ... Passed to [ggplot2::scale_fill_manual()] /
#'   [ggplot2::discrete_scale()].
#' @return A ggplot2 scale.
#' @export
scale_fill_hexrows <- function(counts = NULL, ...) {
  if (is.null(counts)) {
    return(ggplot2::discrete_scale("fill", palette = .hexrows_pal_n,
                                   na.value = .hx$unknown, ...))
  }
  lv <- names(counts)
  ggplot2::scale_fill_manual(values = stats::setNames(.hexrows_pal(lv), lv),
                             breaks = lv, labels = .count_labels(counts),
                             na.value = .hx$unknown, ...)
}

.hexrows_pal <- function(lv) {
  cols <- character(length(lv))
  named <- lv %in% c("Other", "unknown")
  k <- sum(!named)
  hues <- if (k <= length(.hx$hues)) .hx$hues[seq_len(k)] else
    c(.hx$hues, scales::hue_pal()(k - length(.hx$hues)))
  cols[!named] <- hues
  cols[lv == "Other"] <- .hx$other
  cols[lv == "unknown"] <- .hx$unknown
  cols
}

.hexrows_pal_n <- function(n) {
  if (n <= length(.hx$hues)) .hx$hues[seq_len(n)] else
    c(.hx$hues, scales::hue_pal()(n - length(.hx$hues)))
}

#' Save a hexrows figure with a page sized to its rows and bins
#'
#' The page is `unit_in` inches per hexagon unit (0.072 as in the Python
#' figure), plus room for the title, subtitle and legend.
#'
#' @param filename Output file (png, pdf, svg, ...).
#' @param plot A [hexrows()] plot.
#' @param unit_in Inches per hexagon unit (must match the one used to build
#'   the plot for the labels to fit).
#' @param dpi Resolution for raster output.
#' @param ... Passed to [ggplot2::ggsave()].
#' @return The file name (invisibly).
#' @export
ggsave_hexrows <- function(filename, plot = ggplot2::last_plot(), unit_in = NULL,
                           dpi = 300, ...) {
  hd <- attr(plot, "hexrows")
  if (is.null(hd)) stop("`plot` is not a hexrows() plot.", call. = FALSE)
  u <- unit_in %||% hd$meta$unit_in
  w <- diff(hd$meta$xlim) * u + 24 / 72
  b <- ggplot2::ggplot_build(plot)
  nleg <- tryCatch({
    sc <- b$plot$scales$get_scales("fill")
    length(sc$get_breaks())
  }, error = function(e) 0)
  leg_rows <- if (nleg > 0) ceiling(nleg / 3) else 0
  sub <- plot$labels$subtitle %||% ""
  sub_lines <- if (is.null(sub) || identical(sub, "")) 0 else length(strsplit(sub, "\n")[[1]])
  # long subtitle lines wrap in neither ggplot nor here: count characters
  head_in <- 0.35 + 0.14 * sub_lines
  leg_in <- if (leg_rows > 0) 0.3 + 0.15 * leg_rows else 0
  h <- diff(hd$meta$ylim) * u + head_in + leg_in + 18 / 72
  ggplot2::ggsave(filename, plot, width = w, height = h, units = "in", dpi = dpi,
                  bg = .hx$page, limitsize = FALSE, ...)
  invisible(filename)
}
