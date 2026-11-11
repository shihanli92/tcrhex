#' Paired TCRdist radar for one query
#'
#' A compact, single-query version of [hexrows()]: the query sits in a box in
#' the centre; each paired-TCRdist bin is a hexagonal frame placed on its
#' dashed distance ring (0, 12, 24, ... up to `radius`), spread evenly around
#' the query and joined to it by a spoke. Inside each frame, one small hexagon
#' per clone (up to `max_hex`), coloured by any categorical column, banded by
#' level, exactly as in the rows (largest-remainder split, closest clones
#' first). Each frame is labelled with its bin, its number of clones and
#' donors, and how many are shown when some are hidden.
#'
#' @inheritParams hexrows_data
#' @param query One query TCR: `"TRAV:CDR3a:TRBV:CDR3b"` (optionally named) or
#'   a one-row query table, e.g. `top_public(cells, 1)` or
#'   `top_expanded(cells, 1)`.
#' @param radius Search radius; the default 48 keeps the plot compact.
#' @param max_hex Maximum hexagons per frame.
#' @param logos,gene_logos CDR3 logos on the upper edges / V and J gene logos
#'   on the lower edges of each frame (off by default; frames are small).
#' @param donor_noun Word used in the labels, e.g. `"people"` or `"mice"`.
#' @param title,caption Plot title and caption (none by default). The legend
#'   title reads "<fill>: clones within TCRdist <radius>"; change it with
#'   `+ labs(fill = ...)`.
#' @return A ggplot object, with the data attached as `attr(p, "hexrows")`
#'   (so [ggsave_hexrows()] sizes the page).
#' @export
#' @examples
#' cells <- tcrhex_example_cells(n_donors = 20)   # synthetic
#' hexradar(cells, ggplot2::aes(fill = hla_carrier(hla, "A*02")),
#'          query = "TRAV27:CAGAGSQGNLIF:TRBV6-5:CASSYSGGSYEQYF")
hexradar <- function(cells, mapping = ggplot2::aes(), query, radius = 48, bin_width = 12,
                     max_hex = 50, cumulative = FALSE, logos = FALSE, gene_logos = FALSE,
                     logo_trim = 1, summarise_fill = summarise_unique, fold_levels = 8,
                     clones = NULL, hits = NULL, organism = NULL, donor_noun = "donors",
                     title = NULL, caption = NULL, unit_in = 0.072) {
  rd <- hexradar_data(cells, mapping, query = query, radius = radius, bin_width = bin_width,
                      max_hex = max_hex, cumulative = cumulative, logos = logos,
                      gene_logos = gene_logos, logo_trim = logo_trim,
                      summarise_fill = summarise_fill, fold_levels = fold_levels,
                      clones = clones, hits = hits, organism = organism,
                      donor_noun = donor_noun, unit_in = unit_in)
  .hexradar_plot(rd, title = title, caption = caption)
}

#' @rdname hexradar
#' @return `hexradar_data()`: a list of tibbles (`frames`, `hexes`, `logos`,
#'   `labels`, `rings`, `ring_labels`, `spokes`, `box`, `box_text`, `bins`,
#'   `hits`, `levels`) and `meta`.
#' @export
hexradar_data <- function(cells, mapping = ggplot2::aes(), query, radius = 48, bin_width = 12,
                          max_hex = 50, cumulative = FALSE, logos = FALSE, gene_logos = FALSE,
                          logo_trim = 1, summarise_fill = summarise_unique, fold_levels = 8,
                          clones = NULL, hits = NULL, organism = NULL,
                          donor_noun = "donors", unit_in = 0.072) {
  organism <- organism %||% attr(query, "organism") %||% attr(hits, "organism") %||% "human"
  q <- as_queries(query, organism = organism)
  if (nrow(q) != 1) {
    stop("hexradar() draws one query; got ", nrow(q), ". Use hexrows() for several.",
         call. = FALSE)
  }
  hd <- hexrows_data(cells, mapping, queries = query, radius = radius, bin_width = bin_width,
                     max_hex = max_hex, logos = FALSE, gene_logos = FALSE,
                     summarise_fill = summarise_fill, fold_levels = fold_levels,
                     clones = clones, hits = hits, unit_in = unit_in, organism = organism,
                     cumulative = cumulative)
  h <- hd$hits
  bl <- bin_labels(bin_width, radius)
  nb <- nrow(bl)
  c30 <- cos(pi / 6)
  logo_out <- .hx$logo_offset + .hx$logo_height

  # ---- per-bin frames: size from the shown clones (as in hexrows(logos = FALSE)) ----
  fr <- tibble::tibble(bin = bl$bin, edge = bl$hi, label = bl$label)
  fr$n_all <- vapply(fr$bin, function(b) sum(h$bin == b), 1L)
  fr$n_shown <- vapply(fr$bin, function(b) sum(h$bin == b & h$shown), 1L)
  fr$donors <- vapply(fr$bin, function(b) dplyr::n_distinct(h$donor[h$bin == b]), 1L)
  fr$H <- (vapply(pmax(fr$n_shown, 1L), .patch_extent, 1) + 1 + .hx$pad) / c30
  fr$theta <- pi / 2 - 2 * pi * (fr$bin) / nb          # clockwise from 12 o'clock

  # ---- query box (text in mono) ----
  cw <- function(size) 0.6 * size / 72 / unit_in
  lh <- function(size) 1.3 * size / 72 / unit_in
  vw <- max(nchar(q$va), nchar(q$vb))
  box_lines <- c(if (!identical(q$name, paste(q$va, q$ca, q$vb, q$cb, sep = ":"))) q$name else "query",
                 paste0(formatC(q$va, width = -vw), "  ", q$ca),
                 paste0(formatC(q$vb, width = -vw), "  ", q$cb))
  box_w <- max(nchar(box_lines)) * cw(.hx$size_label + 0.5) + 2.4
  box_h <- length(box_lines) * lh(.hx$size_label + 0.5) + 1.6
  box_r <- function(th) pmin(box_w / 2 / pmax(abs(cos(th)), 1e-9), box_h / 2 / pmax(abs(sin(th)), 1e-9))

  # ---- place frames on their rings: r = R0 + edge * s, no overlaps ----
  ext_up <- if (logos) logo_out else 0
  ext_dn <- if (gene_logos) logo_out else 0
  clear <- fr$H + pmax(ext_up, ext_dn)
  R0 <- box_r(fr$theta[1]) + clear[1] + 1.2
  s <- 2 * max(clear) / bin_width / 2
  pos <- function(R0, s) list(x = (R0 + fr$edge * s) * cos(fr$theta),
                              y = (R0 + fr$edge * s) * sin(fr$theta))
  for (it in 1:400) {
    p <- pos(R0, s)
    r <- sqrt(p$x^2 + p$y^2)
    ok_box <- all(r >= box_r(fr$theta) + clear + 1.2)
    nxt <- c(seq_len(nb)[-1], 1)
    dd <- sqrt((p$x - p$x[nxt])^2 + (p$y - p$y[nxt])^2)
    ok_pair <- all(dd >= clear + clear[nxt] + 1.6)
    if (ok_box && ok_pair) break
    if (!ok_box) R0 <- R0 * 1.03
    if (!ok_pair) s <- s * 1.04
  }
  p <- pos(R0, s)
  fr$x <- p$x
  fr$y <- p$y
  ring_r <- R0 + c(0, bl$hi[-1]) * s

  # ---- frames, hexagons (moved from hexrows_data's row layout), logos ----
  fv <- hex_vertices(fr$x, fr$y, fr$H)
  frames <- tibble::tibble(bin = fr$bin[fv$id], group = fv$id, x = fv$x, y = fv$y)
  hexes <- hd$hexes
  if (nrow(hexes)) {
    k <- match(hexes$bin, fr$bin)
    hexes$x <- hexes$x - hexes$bin * hd$meta$px + fr$x[k]
    hexes$y <- hexes$y + fr$y[k]
  }
  logo_df <- NULL
  if (logos || gene_logos) {
    gl <- if (logos) .glyphs() else NULL
    wg <- if (gene_logos) .word_glyphs() else NULL
    out <- list()
    for (i in which(fr$n_all > 0)) {
      hb <- h[h$bin == fr$bin[i], ]
      H <- fr$H[i]
      cx <- fr$x[i]
      cy <- fr$y[i]
      if (logos) {
        apex <- c(cx, cy + H)
        lu <- c(cx - H * c30, cy + H / 2)
        ru <- c(cx + H * c30, cy + H / 2)
        for (ch in c("a", "b")) {
          f <- logo_freqs(logo_align(hb[[paste0("c", ch)]], logo_trim))
          pp <- if (ch == "a") .logo_polys(f, lu, apex, glyphs = gl) else
            .logo_polys(f, apex, ru, glyphs = gl)
          if (!is.null(pp)) out[[length(out) + 1]] <- dplyr::mutate(pp, bin = fr$bin[i], chain = ch)
        }
      }
      if (gene_logos) {
        ll <- c(cx - H * c30, cy - H / 2)
        bot <- c(cx, cy - H)
        lr <- c(cx + H * c30, cy - H / 2)
        for (ch in c("a", "b")) {
          ends <- if (ch == "a") list(ll, bot) else list(bot, lr)
          jg <- hb[[paste0("j", ch)]]
          has_j <- !is.null(jg) && any(!is.na(jg))
          parts <- list(list(g = hb[[paste0("v", ch)]], from = 0, to = if (has_j) 0.5 else 1, t = "v"))
          if (has_j) parts[[2]] <- list(g = jg, from = 0.5, to = 1, t = "j")
          for (pt in parts) {
            pp <- .word_logo_polys(gene_words(pt$g), ends[[1]], ends[[2]], glyphs = wg,
                                   from = pt$from, to = pt$to)
            if (!is.null(pp) && nrow(pp)) {
              out[[length(out) + 1]] <- dplyr::mutate(pp, bin = fr$bin[i],
                                                      chain = paste0(pt$t, ch, "_genes"))
            }
          }
        }
      }
    }
    logo_df <- dplyr::bind_rows(out)
    if (nrow(logo_df)) {
      logo_df$group <- as.integer(factor(paste(logo_df$bin, logo_df$chain, logo_df$glyph)))
      logo_df$subgroup <- logo_df$contour
      logo_df <- logo_df[order(logo_df$group, logo_df$contour, logo_df$order), ]
    }
  }

  # ---- frame labels, outside each frame along its spoke ----
  cu <- cos(fr$theta)
  su <- sin(fr$theta)
  gap <- fr$H + ifelse(su > 0.2, ext_up, ifelse(su < -0.2, ext_dn, 0)) + 1
  lab1 <- if (isTRUE(cumulative)) ifelse(fr$bin == 0, "TCRdist 0 (identical)",
                                         paste0("TCRdist \u2264 ", fr$edge)) else
    ifelse(fr$bin == 0, "TCRdist 0 (identical)", paste("TCRdist", fr$label))
  body <- sprintf("%s clone%s \u00b7 %s %s", .fmt(fr$n_all), ifelse(fr$n_all == 1, "", "s"),
                  .fmt(fr$donors), donor_noun)
  tail <- ifelse(fr$n_all > fr$n_shown, paste0("\n", fr$n_shown, " shown"), "")
  labels <- tibble::tibble(bin = fr$bin, x = fr$x + gap * cu, y = fr$y + gap * su,
                           label = paste0(lab1, "\n", body, tail),
                           hjust = ifelse(cu > 0.3, 0, ifelse(cu < -0.3, 1, 0.5)),
                           vjust = ifelse(su > 0.3, 0, ifelse(su < -0.3, 1, 0.5)))

  # ---- rings, ring labels, spokes, box ----
  tt <- seq(0, 2 * pi, length.out = 361)
  rings <- dplyr::bind_rows(lapply(seq_along(ring_r), function(i)
    tibble::tibble(ring = i, x = ring_r[i] * cos(tt), y = ring_r[i] * sin(tt))))
  a_lab <- pi / 2 + pi / nb                   # between the last bin and bin 0
  ring_labels <- tibble::tibble(x = ring_r * cos(a_lab), y = ring_r * sin(a_lab),
                                label = c(utils::head(c(0, bl$hi[-1]), -1),
                                          paste("TCRdist", radius)))
  ring_labels$label <- as.character(ring_labels$label)
  spokes <- tibble::tibble(bin = fr$bin, x = 0, y = 0, xend = fr$x, yend = fr$y)
  box <- .rounded_rect(0, 0, box_w, box_h, 0.9)
  box_text <- tibble::tibble(x = 0, y = (rev(seq_along(box_lines)) - (length(box_lines) + 1) / 2) *
                               lh(.hx$size_label + 0.5), label = box_lines)

  # ---- limits ----
  lab_w <- max(nchar(c(lab1, body))) * cw(.hx$size_label)
  xs <- c(rings$x, fr$x - fr$H, fr$x + fr$H,
          labels$x + ifelse(labels$hjust == 0, lab_w, ifelse(labels$hjust == 1, -lab_w, 0)),
          labels$x - ifelse(labels$hjust == 0.5, lab_w / 2, 0),
          labels$x + ifelse(labels$hjust == 0.5, lab_w / 2, 0))
  lab_h <- 3 * lh(.hx$size_label)
  ys <- c(rings$y, fr$y - fr$H - ext_dn, fr$y + fr$H + ext_up,
          labels$y + ifelse(labels$vjust == 0, lab_h, ifelse(labels$vjust == 1, -lab_h, lab_h / 2)),
          labels$y - ifelse(labels$vjust == 0.5, lab_h / 2, 0))
  xlim <- range(xs) + c(-1, 1)
  ylim <- range(ys) + c(-1, 1)
  m <- hd$meta
  m$xlim <- xlim
  m$ylim <- ylim
  m$radar <- TRUE
  m$radius <- radius
  m$cumulative <- isTRUE(cumulative)
  structure(list(frames = frames, hexes = hexes, logos = logo_df, labels = labels,
                 rings = rings, ring_labels = ring_labels, spokes = spokes, box = box,
                 box_text = box_text, bins = hd$bins, hits = h, levels = hd$levels,
                 pattern_levels = hd$pattern_levels, rows = hd$rows, meta = m),
            class = "hexradar_data")
}

#' @export
print.hexradar_data <- function(x, ...) {
  cat(sprintf("<hexradar_data> 1 query x %d bins, %s hits (%s shown), %d fill levels\n",
              x$meta$n_bins, .fmt(nrow(x$hits)), .fmt(sum(x$hits$shown)), nrow(x$levels)))
  invisible(x)
}

.rounded_rect <- function(cx, cy, w, h, r) {
  r <- min(r, w / 2, h / 2)
  arc <- function(x0, y0, a0) {
    t <- seq(a0, a0 + pi / 2, length.out = 10)
    data.frame(x = x0 + r * cos(t), y = y0 + r * sin(t))
  }
  hw <- w / 2 - r
  hh <- h / 2 - r
  tibble::as_tibble(rbind(arc(cx + hw, cy + hh, 0), arc(cx - hw, cy + hh, pi / 2),
                          arc(cx - hw, cy - hh, pi), arc(cx + hw, cy - hh, 3 * pi / 2)))
}

.hexradar_plot <- function(rd, title = NULL, caption = NULL) {
  m <- rd$meta
  env <- .default_scale_env(rd)
  hexes <- rd$hexes
  names(hexes)[names(hexes) == "level"] <- ".level"
  names(hexes)[names(hexes) == "pattern"] <- ".pattern"
  has_fill <- !is.null(m$fill_label)
  has_pat <- !is.null(m$pattern_label)
  hex_aes <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$group)
  if (has_fill) hex_aes$fill <- rlang::quo(.data$.level)
  if (has_pat) hex_aes$pattern <- rlang::quo(.data$.pattern)
  lw <- .pt2lw
  p <- ggplot2::ggplot(environment = env) +
    ggplot2::geom_path(data = rd$rings, ggplot2::aes(x = .data$x, y = .data$y, group = .data$ring),
                       colour = .hx$frame_line, linewidth = lw(0.6), linetype = "44") +
    ggplot2::geom_text(data = rd$ring_labels, ggplot2::aes(x = .data$x, y = .data$y,
                       label = .data$label), colour = .hx$muted, size = .hx$size_header / ggplot2::.pt,
                       hjust = 1.1, vjust = -0.3) +
    ggplot2::geom_segment(data = rd$spokes, ggplot2::aes(x = .data$x, y = .data$y,
                          xend = .data$xend, yend = .data$yend),
                          colour = "#8f8d85", linewidth = lw(0.9)) +
    ggplot2::geom_polygon(data = rd$frames, ggplot2::aes(x = .data$x, y = .data$y,
                          group = .data$group), fill = .hx$frame_fill, colour = .hx$frame_line,
                          linewidth = lw(.hx$frame_lwd_pt))
  p <- p + if (has_fill) {
    geom_hextile(data = hexes, mapping = hex_aes, colour = .hx$frame_fill,
                 linewidth = lw(.hx$hex_lwd_pt))
  } else {
    geom_hextile(data = hexes, mapping = hex_aes, colour = .hx$frame_fill, fill = .hx$hues[1],
                 linewidth = lw(.hx$hex_lwd_pt))
  }
  if (!is.null(rd$logos) && nrow(rd$logos)) {
    p <- p + ggplot2::geom_polygon(data = rd$logos, ggplot2::aes(x = .data$x, y = .data$y,
                                   group = .data$group, subgroup = .data$subgroup),
                                   fill = .hx$ink, colour = NA)
  }
  leg <- function(lab) if (is.null(lab)) NULL else
    sprintf("%s: clones within TCRdist %s%s", lab, m$radius,
            if (isTRUE(m$cumulative)) " (cumulative)" else "")
  p <- p +
    ggplot2::geom_polygon(data = rd$box, ggplot2::aes(x = .data$x, y = .data$y),
                          fill = "white", colour = .hx$frame_line, linewidth = lw(1)) +
    ggplot2::geom_text(data = rd$box_text, ggplot2::aes(x = .data$x, y = .data$y,
                       label = .data$label), family = "mono", colour = .hx$ink,
                       size = (.hx$size_label + 0.5) / ggplot2::.pt) +
    ggplot2::geom_text(data = rd$labels, ggplot2::aes(x = .data$x, y = .data$y,
                       label = .data$label, hjust = .data$hjust, vjust = .data$vjust),
                       colour = .hx$ink2, size = .hx$size_label / ggplot2::.pt, lineheight = 1.05) +
    ggplot2::coord_fixed(xlim = m$xlim, ylim = m$ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(title = title, caption = caption, fill = leg(m$fill_label),
                  pattern = leg(m$pattern_label)) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, byrow = TRUE),
                    pattern = ggplot2::guide_legend(ncol = 3, byrow = TRUE)) +
    theme_hexrows() +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0.5, colour = .hx$ink,
                                                        size = 8, margin = ggplot2::margin(t = 6)),
                   plot.caption.position = "plot")
  attr(p, "hexrows") <- rd
  p
}
