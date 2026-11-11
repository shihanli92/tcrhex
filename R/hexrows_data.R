# Style constants (Python figure)
.hx <- list(
  page = "#fcfcfb", frame_fill = "#f3f2ee", frame_line = "#c3c2b7",
  ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781", row_line = "#d9d8d0",
  other = "#b4b3ac", unknown = "#d6d5ce",
  hues = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
           "#6250d6", "#e34948"),
  frame_lwd_pt = 0.8, hex_lwd_pt = 0.6,
  logo_offset = 0.35, logo_height = 2.6, pad = 0.7,
  size_label = 7, size_tag = 5.5, size_count = 5.5, size_header = 6.5
)
.pt2lw <- function(pt) pt / (72.27 / 25.4) / 0.75 * (72.27 / 72) # points -> ggplot linewidth

#' The data behind a hexrows figure
#'
#' The pure data step of [hexrows()]: computes clones, paired TCRdist hits,
#' bins, the clones shown in each bin, their honeycomb positions, the CDR3
#' logos and all labels, as tidy tibbles. `hexrows()` only draws them.
#'
#' @param cells A cell table from a reader ([read_screpertoire()],
#'   [read_10x_contigs()], [read_airr()], [cells_from_df()]).
#' @param mapping Aesthetics made with [ggplot2::aes()]: `fill` (any
#'   categorical expression of the cell columns) and optionally `pattern` (a
#'   second categorical expression drawn as a texture). Expressions are
#'   evaluated per cell; every column of the query table is also available as
#'   `.row_<column>` (e.g. `.row_allele`), and `.row` is the row number. The
#'   per-cell values are then summarised per clone with `summarise_fill`.
#' @param queries Query TCRs: `NULL` (default: `top_public(cells, 20)`), a
#'   query table such as from [top_public()], or a character vector
#'   `"TRAV:CDR3a:TRBV:CDR3b"` (optionally named) - see [as_queries()].
#' @param radius Paired TCRdist search radius.
#' @param bin_width Width of the distance bins (bin 0 is `d == 0`).
#' @param max_hex Maximum hexagons (clones) shown per bin.
#' @param logos Draw CDR3 logos on the frames.
#' @param logo_trim Residues trimmed at each CDR3 end before the logo.
#' @param logo_font `NULL` for the built-in DejaVu Sans Mono Bold outlines,
#'   or a font family outlined with `systemfonts::glyph_outline()` (systemfonts
#'   version 1.2 or later).
#' @param summarise_fill Rule turning the values of a clone's cells into one
#'   value, e.g. [summarise_unique()] (default), [summarise_tissue()],
#'   [summarise_most_common()] or your own `function(values)`.
#' @param fold_levels For a non-factor fill: keep the `fold_levels - 1` most
#'   common levels (over all rows) and fold the rest into `"Other"` when there
#'   are more than `fold_levels` levels. Factors keep all their levels (use
#'   forcats to lump them).
#' @param dominant Rows from [top_public()] whose alpha concordance is below
#'   this are muted and tagged "no dominant alpha".
#' @param row_tag Optional character vector, one per query row, shown as the
#'   third (small) label line, replacing the default tag.
#' @param clones,hits Optional precomputed [tcr_clones()] and
#'   [tcrdist_paired()] results (hits must cover `radius`). With
#'   `cells = NULL`, `hits` may be any per-clone table (`query`, `d`, `d_a`,
#'   `d_b`, `donor`, `ca`, `cb`, `cells` and metadata columns); the
#'   aesthetics are then evaluated on its (already per-clone) columns.
#' @param unit_in Inches per hexagon unit, used to place the text (as by
#'   [ggsave_hexrows()]).
#' @param cumulative If `TRUE`, bin k holds every clone with `d` up to the
#'   bin's upper edge (bin 0: `d == 0`; then `d <= 12`, `d <= 24`, ...), so
#'   the hexagon composition, the largest-remainder split, the count and the
#'   CDR3 logos accumulate as the distance grows. Headers read "0", "<=12",
#'   "<=24", .... The legend and right-hand counts still count each clone once
#'   per row.
#' @param organism `"human"` or `"mouse"` (tcrdist3 V-gene loop table). By
#'   default taken from `queries` / `hits` when they come from
#'   [top_public()] / [tcrdist_paired()] with an organism, else `"human"`.
#' @return A list of tibbles: `frames`, `hexes`, `logos`, `counts`, `labels`,
#'   `lines`, `bins` (the per-bin table: `row`, `bin`, `level`, `clones`,
#'   `shown`), `hits` (every hit with its level, display `bin`, its own
#'   distance bin `own_bin` and whether it is shown; with `cumulative = TRUE`
#'   a clone appears in every bin from its own bin onwards),
#'   `rows` (the query table), `levels` (level counts), and `meta`.
#' @export
#' @examples
#' cells <- tcrhex_example_cells()
#' hd <- hexrows_data(cells, ggplot2::aes(fill = disease), queries = top_public(cells, 2))
#' hd$bins
hexrows_data <- function(cells, mapping = ggplot2::aes(), queries = NULL, radius = 108,
                         bin_width = 12, max_hex = 19, logos = TRUE, logo_trim = 1,
                         logo_font = NULL, summarise_fill = summarise_unique,
                         fold_levels = 8, dominant = 0.10, row_tag = NULL,
                         clones = NULL, hits = NULL, unit_in = 0.072,
                         organism = NULL, cumulative = FALSE) {
  organism <- organism %||% attr(queries, "organism") %||% attr(hits, "organism") %||% "human"
  organism <- match.arg(organism, c("human", "mouse"))
  if (is.null(cells)) {
    # clone-level input: precomputed hits carrying the clone metadata
    if (is.null(hits) || is.null(queries)) {
      stop("Without `cells`, give both `queries` and `hits` (one row per query x clone, ",
           "with the clone's metadata columns).", call. = FALSE)
    }
    from_top <- FALSE
    q <- as_queries(queries, organism = organism)
    hits <- tibble::as_tibble(hits)
    if (!".clone_id" %in% names(hits)) hits$.clone_id <- seq_len(nrow(hits))
    if (!"cells" %in% names(hits)) hits$cells <- 1L
  } else {
    cells <- .as_cells(cells)
    if (is.null(clones)) clones <- tcr_clones(cells)
    from_top <- is.null(queries)
    if (from_top) queries <- top_public(cells, 20, organism = organism)
    q <- as_queries(queries, organism = organism)
    if (is.null(hits)) hits <- tcrdist_paired(q, clones, radius = radius, organism = organism)
  }
  q$.row <- seq_len(nrow(q))
  hits <- hits[hits$d <= radius & hits$query %in% q$name, , drop = FALSE]
  hits$.row <- q$.row[match(hits$query, q$name)]
  hits$.bin <- tcrdist_bin(hits$d, bin_width, radius)

  # ---- aesthetics: evaluate per cell, summarise per (row, clone) ----
  fill_q <- mapping$fill
  pat_q <- mapping$pattern
  cc <- if (is.null(cells)) NULL else attr(clones, "cell_clone")
  ev <- .eval_aes(list(fill = fill_q, pattern = pat_q), hits, cells, cc, q,
                  summarise_fill)
  hits$.fill_raw <- ev$fill
  hits$.pattern_raw <- ev$pattern
  lv <- .make_levels(hits$.fill_raw, fold_levels)
  hits$.level <- lv$value
  plv <- if (is.null(pat_q)) NULL else .make_levels(hits$.pattern_raw, Inf)
  hits$.pattern <- if (is.null(plv)) factor("none") else plv$value
  hits$.own_bin <- hits$.bin
  hits_all <- hits   # one row per (row, clone)

  # ---- cumulative: a clone counts in its own bin and every wider one ----
  if (isTRUE(cumulative) && nrow(hits)) {
    last <- max(bin_labels(bin_width, radius)$bin)
    reps <- last - hits$.bin + 1L
    hits <- hits[rep(seq_len(nrow(hits)), reps), , drop = FALSE]
    hits$.bin <- as.integer(sequence(reps, from = hits_all$.bin))
  }

  # ---- choose the shown clones per bin ----
  hits <- hits[order(hits$.row, hits$.bin, as.integer(hits$.level), hits$d,
                     -hits$cells, hits$.clone_id), , drop = FALSE]
  hits$.shown <- FALSE
  key <- paste(hits$.row, hits$.bin)
  idx_by_bin <- split(seq_len(nrow(hits)), factor(key, levels = unique(key)))
  for (ix in idx_by_bin) {
    lev <- hits$.level[ix]
    n <- table(lev)
    n <- n[n > 0]
    k <- allot(stats::setNames(as.integer(n), names(n)), max_hex)
    for (l in names(k)) {
      il <- ix[as.character(lev) == l]
      hits$.shown[il[seq_len(k[[l]])]] <- TRUE
    }
  }
  bins <- dplyr::count(hits, row = .data$.row, bin = .data$.bin, level = .data$.level,
                       name = "clones")
  sh <- dplyr::count(hits[hits$.shown, ], row = .data$.row, bin = .data$.bin,
                     level = .data$.level, name = "shown")
  bins <- dplyr::left_join(bins, sh, by = c("row", "bin", "level"))
  bins$shown[is.na(bins$shown)] <- 0L
  bins <- bins[order(bins$row, bins$bin, as.integer(bins$level)), ]

  # ---- geometry ----
  c30 <- cos(pi / 6)
  Hmax <- (.patch_extent(max_hex) + 1 + .hx$pad) / c30
  bl <- bin_labels(bin_width, radius)
  nb <- nrow(bl)
  nr <- nrow(q)
  logo_out <- (.hx$logo_offset + .hx$logo_height)
  px <- 2 * Hmax * c30 + if (logos) 2 * logo_out * 0.5 + 1.6 else 2.2
  py <- 2 * Hmax + if (logos) logo_out * c30 + 2.2 else 2.2
  grid <- expand.grid(.bin = bl$bin, .row = seq_len(nr))
  grid$cx <- grid$.bin * px
  grid$cy <- -(grid$.row - 1) * py
  tot <- dplyr::count(hits, .row = .data$.row, .bin = .data$.bin, name = "n_all")
  nsh <- dplyr::count(hits[hits$.shown, ], .row = .data$.row, .bin = .data$.bin,
                      name = "n_shown")
  grid <- dplyr::left_join(grid, tot, by = c(".row", ".bin"))
  grid <- dplyr::left_join(grid, nsh, by = c(".row", ".bin"))
  grid$n_all[is.na(grid$n_all)] <- 0L
  grid$n_shown[is.na(grid$n_shown)] <- 0L
  grid$H <- if (logos) Hmax else (vapply(pmax(grid$n_shown, 1), .patch_extent, 1) + 1 + .hx$pad) / c30
  grid$frame_id <- seq_len(nrow(grid))

  fv <- hex_vertices(grid$cx, grid$cy, grid$H)
  frames <- tibble::tibble(row = grid$.row[fv$id], bin = grid$.bin[fv$id],
                           group = fv$id, x = fv$x, y = fv$y)

  # small hexagons
  shown <- hits[hits$.shown, , drop = FALSE]
  shown <- shown[order(shown$.row, shown$.bin, as.integer(shown$.level),
                       as.integer(shown$.pattern), shown$d, -shown$cells,
                       shown$.clone_id), , drop = FALSE]
  shown$slot <- stats::ave(seq_len(nrow(shown)), shown$.row, shown$.bin, FUN = seq_along)
  nshown <- stats::ave(seq_len(nrow(shown)), shown$.row, shown$.bin, FUN = length)
  lay <- lapply(sort(unique(nshown)), function(k) dplyr::mutate(hex_layout(k), k = k))
  lay <- dplyr::bind_rows(lay)
  m <- match(paste(nshown, shown$slot), paste(lay$k, lay$slot))
  gi <- match(paste(shown$.row, shown$.bin), paste(grid$.row, grid$.bin))
  shown$hx <- grid$cx[gi] + lay$x[m]
  shown$hy <- grid$cy[gi] + lay$y[m]
  hv <- hex_vertices(shown$hx, shown$hy, 1)
  hexes <- tibble::tibble(row = shown$.row[hv$id], bin = shown$.bin[hv$id],
                          group = hv$id, .clone_id = shown$.clone_id[hv$id],
                          d = shown$d[hv$id], level = shown$.level[hv$id],
                          pattern = shown$.pattern[hv$id], x = hv$x, y = hv$y)

  # logos
  logo_df <- NULL
  if (logos) {
    gl <- .glyphs(logo_font)
    ne <- grid[grid$n_all > 0, ]
    out <- vector("list", nrow(ne) * 2)
    k <- 0
    for (i in seq_len(nrow(ne))) {
      hb <- hits[hits$.row == ne$.row[i] & hits$.bin == ne$.bin[i], ]
      H <- ne$H[i]
      cx <- ne$cx[i]
      cy <- ne$cy[i]
      apex <- c(cx, cy + H)
      lu <- c(cx - H * c30, cy + H / 2)
      ru <- c(cx + H * c30, cy + H / 2)
      for (ch in c("a", "b")) {
        fr <- logo_freqs(logo_align(hb[[paste0("c", ch)]], logo_trim))
        p <- if (ch == "a") .logo_polys(fr, lu, apex, .hx$logo_offset, .hx$logo_height, glyphs = gl)
             else .logo_polys(fr, apex, ru, .hx$logo_offset, .hx$logo_height, glyphs = gl)
        if (is.null(p)) next
        k <- k + 1
        p$row <- ne$.row[i]
        p$bin <- ne$.bin[i]
        p$chain <- ch
        out[[k]] <- p
      }
    }
    logo_df <- dplyr::bind_rows(out)
    if (nrow(logo_df)) {
      logo_df$group <- as.integer(factor(paste(logo_df$row, logo_df$bin, logo_df$chain,
                                               logo_df$glyph)))
      logo_df$subgroup <- logo_df$contour
      logo_df <- logo_df[order(logo_df$group, logo_df$contour, logo_df$order), ]
    }
  }

  # counts
  counts <- tibble::tibble(row = grid$.row, bin = grid$.bin,
                           x = grid$cx + grid$H * c30 * 0.92,
                           y = grid$cy - grid$H * 0.62, label = as.character(grid$n_all),
                           colour = ifelse(grid$n_all > grid$n_shown, .hx$ink2, .hx$muted))

  # row lines
  lines <- tibble::tibble(row = seq_len(nr), x = 0, xend = max(bl$bin) * px,
                          y = -(seq_len(nr) - 1) * py, yend = -(seq_len(nr) - 1) * py)

  # ---- labels ----
  labs <- .row_labels(q, hits, from_top || !is.null(attr(queries, "source")),
                      dominant, row_tag)
  cw <- function(size) 0.6 * size / 72 / unit_in       # mono char width (units)
  lh <- function(size) 1.25 * size / 72 / unit_in      # line height (units)
  left_chars <- max(nchar(c(labs$line1, labs$line2)), 0)
  tag_w <- max(nchar(labs$tag), 0) * cw(.hx$size_tag)
  left_w <- max(left_chars * cw(.hx$size_label), tag_w)
  x_left <- -Hmax * c30 - (if (logos) logo_out * 0.5 else 0) - 2 - left_w
  ry <- -(labs$.row - 1) * py
  L1 <- lh(.hx$size_label)
  lab_left <- dplyr::bind_rows(
    tibble::tibble(row = labs$.row, part = "line1", x = x_left, y = ry + L1 * 0.9,
                   label = labs$line1, colour = ifelse(labs$muted, .hx$muted, .hx$ink),
                   size = .hx$size_label),
    tibble::tibble(row = labs$.row, part = "line2", x = x_left, y = ry - L1 * 0.1,
                   label = labs$line2, colour = ifelse(labs$muted, .hx$muted, .hx$ink),
                   size = .hx$size_label),
    tibble::tibble(row = labs$.row, part = "tag", x = x_left, y = ry - L1 * 1.15,
                   label = labs$tag, colour = .hx$muted, size = .hx$size_tag))
  lab_left <- lab_left[!is.na(lab_left$label) & lab_left$label != "", ]
  x_right <- max(bl$bin) * px + Hmax * c30 + 1.6 + 4 * cw(.hx$size_count)
  rr <- dplyr::summarise(dplyr::group_by(hits_all, .data$.row), clones = dplyr::n(),
                         donors = dplyr::n_distinct(.data$donor), .groups = "drop")
  rr <- dplyr::left_join(tibble::tibble(.row = seq_len(nr)), rr, by = ".row")
  rr$clones[is.na(rr$clones)] <- 0L
  rr$donors[is.na(rr$donors)] <- 0L
  lab_right <- tibble::tibble(row = rr$.row, part = "right", x = x_right,
                              y = -(rr$.row - 1) * py,
                              label = sprintf("%s clone%s / %s donor%s", .fmt(rr$clones),
                                              ifelse(rr$clones == 1, "", "s"), .fmt(rr$donors),
                                              ifelse(rr$donors == 1, "", "s")),
                              colour = .hx$ink2, size = .hx$size_label)
  right_w <- max(nchar(lab_right$label)) * cw(.hx$size_label)
  top <- Hmax + (if (logos) logo_out * c30 else 0)
  y_head <- top + 1.2
  head_lab <- if (isTRUE(cumulative)) ifelse(bl$bin == 0, "0", paste0("\u2264", bl$hi)) else bl$label
  axis_lab <- if (isTRUE(cumulative)) "cumulative: all clones with paired TCRdist \u2264 the bin edge" else
    "paired TCRdist to the row's query (bins)"
  header <- dplyr::bind_rows(
    tibble::tibble(row = 0L, part = "bin", x = bl$bin * px, y = y_head, label = head_lab,
                   colour = .hx$ink2, size = .hx$size_header),
    tibble::tibble(row = 0L, part = "axis", x = max(bl$bin) * px / 2,
                   y = y_head + lh(.hx$size_header) * 1.3,
                   label = axis_lab,
                   colour = .hx$muted, size = .hx$size_header))
  labels <- dplyr::bind_rows(lab_left, lab_right, header)

  xlim <- c(x_left - 1, x_right + right_w + 1)
  ylim <- c(-(nr - 1) * py - Hmax - 2.2, y_head + lh(.hx$size_header) * 2.2)
  levels_tab <- dplyr::count(hits_all, level = .data$.level, name = "clones", .drop = FALSE)
  pat_tab <- if (is.null(plv)) NULL else dplyr::count(hits_all, level = .data$.pattern,
                                                       name = "clones", .drop = FALSE)
  rows <- q
  rows$tag <- labs$tag
  rows$muted <- labs$muted
  hits_out <- hits[, c("query", ".row", ".clone_id", "d", "d_a", "d_b", ".bin",
                       ".own_bin", ".level", ".pattern", ".shown",
                       intersect(c("donor", "va", "ca", "vb", "cb", "cells"), names(hits)))]
  names(hits_out)[names(hits_out) %in% c(".row", ".bin", ".own_bin", ".level", ".pattern",
                                         ".shown")] <-
    c("row", "bin", "own_bin", "level", "pattern", "shown")
  structure(list(
    frames = frames, hexes = hexes, logos = logo_df, counts = counts, labels = labels,
    lines = lines, bins = bins, hits = hits_out, rows = rows, levels = levels_tab,
    pattern_levels = pat_tab,
    meta = list(Hmax = Hmax, px = px, py = py, xlim = xlim, ylim = ylim, unit_in = unit_in,
                radius = radius, bin_width = bin_width, max_hex = max_hex, logos = logos,
                logo_trim = logo_trim, n_rows = nr, n_bins = nb,
                cumulative = isTRUE(cumulative),
                fill_label = if (is.null(fill_q)) NULL else rlang::as_label(fill_q),
                pattern_label = if (is.null(pat_q)) NULL else rlang::as_label(pat_q),
                folded = lv$folded, fill_is_factor = lv$is_factor)),
    class = "hexrows_data")
}

#' @export
print.hexrows_data <- function(x, ...) {
  cat(sprintf("<hexrows_data> %d rows x %d bins, %s hits (%s shown), %d fill levels\n",
              x$meta$n_rows, x$meta$n_bins, .fmt(nrow(x$hits)), .fmt(sum(x$hits$shown)),
              nrow(x$levels)))
  cat("components:", paste(names(unclass(x)), collapse = ", "), "\n")
  invisible(x)
}

# Evaluate fill / pattern quosures per cell (of hit clones) and summarise per
# (row, clone). Returns list(fill = <vector aligned with hits>, pattern = ...).
.eval_aes <- function(quos, hits, cells, cc, q, summarise_fill) {
  res <- list(fill = rep("clone", nrow(hits)), pattern = rep("none", nrow(hits)))
  todo <- names(quos)[!vapply(quos, is.null, TRUE)]
  if (!length(todo) || !nrow(hits)) return(res)
  # cell-level frame: one row per (hit, cell)
  hk <- tibble::tibble(.hit = seq_len(nrow(hits)), .clone_id = hits$.clone_id,
                       .row = hits$.row)
  if (is.null(cc)) {
    # clone-level only (no cell map): evaluate on the clone columns of hits
    dat <- hits
    hit_of <- seq_len(nrow(hits))
  } else {
    cc <- cc[cc$.clone_id %in% hits$.clone_id, , drop = FALSE]
    ci <- split(match(cc$cell, cells$cell), cc$.clone_id)
    rows_cells <- ci[as.character(hk$.clone_id)]
    hit_of <- rep(hk$.hit, lengths(rows_cells))
    dat <- cells[unlist(rows_cells, use.names = FALSE), , drop = FALSE]
    dat$.row <- hits$.row[hit_of]
    dat$d <- hits$d[hit_of]
  }
  qcols <- setdiff(names(q), ".row")
  for (cn in qcols) dat[[paste0(".row_", cn)]] <- q[[cn]][match(dat$.row, q$.row)]
  ord <- order(hit_of)
  g <- hit_of[ord]
  first <- !duplicated(g)
  for (a in todo) {
    v <- rlang::eval_tidy(quos[[a]], data = dat)
    if (length(v) == 1) v <- rep(v, nrow(dat))
    if (length(v) != nrow(dat)) {
      stop(sprintf("aes(%s = %s) gave %d values for %d cells.", a,
                   rlang::as_label(quos[[a]]), length(v), nrow(dat)), call. = FALSE)
    }
    gi <- match(g, unique(g))
    s <- .summarise_col(v[ord], gi, max(gi), first, summarise_fill)
    out <- s[match(seq_len(nrow(hits)), unique(g))]
    res[[a]] <- out
  }
  res
}

# Level order: factor levels as given; otherwise by clone count (desc, ties
# alphabetical), folded to `fold - 1` + "Other" when more than `fold` levels.
# NA -> "unknown" (last).
.make_levels <- function(x, fold) {
  is_f <- is.factor(x)
  ch <- as.character(x)
  unk <- is.na(ch) | tolower(ch) %in% c("unknown", "na", "")
  ch[unk] <- "unknown"
  folded <- character()
  if (is_f) {
    lv <- levels(x)
    lv <- lv[lv %in% ch]
    lv <- c(setdiff(lv, "unknown"), setdiff(unique(ch), c(lv, "unknown")))
  } else {
    tab <- table(ch[ch != "unknown"])
    lv <- names(tab)[order(-as.integer(tab), names(tab))]
    if (length(lv) > fold) {
      folded <- lv[fold:length(lv)]
      lv <- c(lv[seq_len(fold - 1)], "Other")
      ch[ch %in% folded] <- "Other"
    }
  }
  if (any(ch == "unknown")) lv <- c(lv, "unknown")
  list(value = factor(ch, levels = lv), folded = folded, is_factor = is_f)
}

# Left-hand row labels.
.row_labels <- function(q, hits, top, dominant, row_tag) {
  n <- nrow(q)
  rank <- if ("rank" %in% names(q)) as.character(q$rank) else NULL
  vbw <- max(nchar(q$vb), nchar(q$va))
  pre_w <- if (is.null(rank)) 0 else max(nchar(rank)) + 2
  pre1 <- if (is.null(rank)) "" else formatC(rank, width = -(pre_w))
  pre2 <- strrep(" ", pre_w)
  line1 <- paste0(pre1, formatC(q$vb, width = -(vbw + 2)), q$cb)
  line2 <- paste0(pre2, formatC(q$va, width = -(vbw + 2)), q$ca)
  tag <- rep("", n)
  muted <- rep(FALSE, n)
  has_conc <- "concordance" %in% names(q) && any(!is.na(q$concordance))
  if (has_conc) {
    conc <- q$concordance
    muted <- !is.na(conc) & conc < dominant
    tag <- ifelse(muted, sprintf("no dominant \u03b1 (top \u03b1 in %s%%)", round(100 * conc)),
                  sprintf("\u03b1 in %s%% of paired people", round(100 * conc)))
    tag[is.na(conc)] <- ""
  }
  extra <- NULL
  if ("epitope" %in% names(q)) extra <- q$epitope
  if ("label" %in% names(q)) extra <- ifelse(is.na(extra %||% NA), q$label,
                                             paste(extra, q$label, sep = " \u00b7 "))
  if (!is.null(extra)) {
    extra[is.na(extra)] <- ""
    tag <- ifelse(tag == "", extra, ifelse(extra == "", tag, paste(tag, extra, sep = " \u00b7 ")))
  }
  if (!has_conc && is.null(extra) && !identical(q$name, paste(q$va, q$ca, q$vb, q$cb, sep = ":"))) {
    tag <- q$name
  }
  if (!is.null(row_tag)) {
    if (length(row_tag) != n) stop("`row_tag` must have one value per query row (", n, ").",
                                   call. = FALSE)
    tag <- ifelse(is.na(row_tag), "", as.character(row_tag))
  }
  tibble::tibble(.row = q$.row, line1 = line1, line2 = line2, tag = tag, muted = muted)
}
