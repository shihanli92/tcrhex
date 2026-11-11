#' Summarise the values of one clone
#'
#' The default rule for summarising a cell-level column per clone: a single
#' distinct value gives that value, several give `"multiple"`. Missing values
#' are ignored (all missing gives `NA`). Write your own with the same signature
#' (`function(values)` returning one value), e.g. [summarise_most_common()].
#'
#' @param values The values of one column for the cells of one clone.
#' @return One value.
#' @export
#' @examples
#' summarise_unique(c("blood", "blood"))
#' summarise_unique(c("blood", "lung", NA))
summarise_unique <- function(values) {
  u <- unique(values[!is.na(values)])
  if (length(u) == 0) return(values[NA_integer_][1])
  if (length(u) == 1) return(u)
  if (is.factor(values)) return(factor("multiple"))
  "multiple"
}

#' @rdname summarise_unique
#' @export
#' @examples
#' summarise_most_common(c("blood", "blood", "lung"))
summarise_most_common <- function(values) {
  v <- values[!is.na(values)]
  if (!length(v)) return(values[NA_integer_][1])
  tab <- table(as.character(v))
  names(tab)[which.max(tab)]
}

#' Collapse cells into clones
#'
#' A clone is donor x TRAV + CDR3a + TRBV + CDR3b (amino acid). Only paired
#' cells (both chains) enter. `cells` is the number of cells; every metadata
#' column is summarised per clone by one rule, `summarise` (default:
#' [summarise_unique()], so a single value stays, several give `"multiple"`).
#' Donor-level columns (one value per donor) therefore pass through unchanged.
#'
#' @param cells A cell table from one of the readers.
#' @param summarise Function applied to the values of every metadata column
#'   within a clone; must return one value.
#' @return A tibble with `.clone_id`, `donor`, `va`, `ca`, `vb`, `cb`,
#'   `cells`, then the summarised metadata columns. J genes (`ja`, `jb`), if
#'   present, are summarised by their most common value (a clone is defined
#'   by V + CDR3, so its cells may rarely differ in J).
#' @export
#' @examples
#' cells <- tcrhex_example_cells()
#' tcr_clones(cells)
tcr_clones <- function(cells, summarise = summarise_unique) {
  cells <- .as_cells(cells)
  p <- cells[!is.na(cells$ca) & !is.na(cells$cb), , drop = FALSE]
  key <- paste(p$donor, p$va, p$ca, p$vb, p$cb, sep = "\r")
  ord <- order(p$donor, p$va, p$ca, p$vb, p$cb, method = "radix")
  p <- p[ord, , drop = FALSE]
  key <- key[ord]
  first <- !duplicated(key)
  g <- cumsum(first)
  ng <- sum(first)
  out <- tibble::tibble(.clone_id = seq_len(ng), donor = p$donor[first],
                        va = p$va[first], ca = p$ca[first], vb = p$vb[first],
                        cb = p$cb[first], cells = tabulate(g, ng))
  meta <- setdiff(names(p), c(.tcr_cols))
  for (m in meta) {
    out[[m]] <- if (m %in% c("ja", "jb")) .most_common_by(p[[m]], g, ng) else
      .summarise_col(p[[m]], g, ng, first, summarise)
  }
  attr(out, "cell_clone") <- tibble::tibble(cell = p$cell, .clone_id = g)
  class(out) <- c("tcrhex_clones", setdiff(class(out), "tcrhex_clones"))
  out
}

# Most common non-NA value per group (ties alphabetical), vectorised; used for
# the J genes, which may differ between cells of one V + CDR3 clone.
.most_common_by <- function(x, g, ng) {
  ok <- !is.na(x)
  out <- rep(NA_character_, ng)
  if (!any(ok)) return(out)
  t <- dplyr::count(tibble::tibble(g = g[ok], x = as.character(x[ok])), .data$g, .data$x)
  t <- t[order(t$g, -t$n, t$x, method = "radix"), ]
  t <- t[!duplicated(t$g), ]
  out[t$g] <- t$x
  out
}

# Summarise column x by group g (1..ng, sorted) with rule f; fast path for the
# default rule.
.summarise_col <- function(x, g, ng, first, f) {
  if (identical(f, summarise_unique)) {
    if (is.list(x)) x <- vapply(x, function(z) paste(z, collapse = ";"), "")
    na <- is.na(x)
    # first non-NA value per group
    idx <- seq_along(x)
    idx[na] <- NA
    fi <- vapply(split(idx, g), function(z) { z <- z[!is.na(z)]; if (length(z)) z[1] else NA_integer_ }, 1L)
    if (length(fi) != ng) stop("internal: group mismatch")
    ref <- x[fi[g]]
    diff <- !na & !is.na(ref) & as.character(x) != as.character(ref)
    multi <- tabulate(g[diff], ng) > 0
    val <- x[fi]
    if (any(multi)) {
      if (is.factor(val)) {
        val <- factor(as.character(val), levels = unique(c(levels(x), "multiple")))
      } else {
        val <- as.character(val)
      }
      val[multi] <- "multiple"
    }
    return(val)
  }
  res <- lapply(split(x, g), f)
  if (all(vapply(res, function(z) is.factor(z), TRUE))) {
    lv <- unique(unlist(lapply(res, levels)))
    return(factor(vapply(res, as.character, ""), levels = lv))
  }
  unlist(res, use.names = FALSE)
}

#' Synthetic example cell table
#'
#' A small SYNTHETIC cell table (invented TCRs, donors and metadata, not real
#' data) for examples and tests.
#' @param n_donors Number of synthetic donors.
#' @param seed Random seed.
#' @return A `tcrhex_cells` table.
#' @export
#' @examples
#' tcrhex_example_cells()
tcrhex_example_cells <- function(n_donors = 6, seed = 1) {
  withr_seed <- function(seed, code) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
    set.seed(seed)
    code
  }
  withr_seed(seed, {
    aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
    mut <- function(s, k) {
      ch <- strsplit(s, "")[[1]]
      for (i in seq_len(k)) {
        pos <- sample(4:(length(ch) - 3), 1)
        ch[pos] <- sample(aa, 1)
      }
      paste(ch, collapse = "")
    }
    base <- list(c("TRAV12-2", "CAVNGGSQGNLIF", "TRBV19", "CASSIRSSYEQYF"),
                 c("TRAV27", "CAGAGSQGNLIF", "TRBV6-5", "CASSYSGGSYEQYF"),
                 c("TRAV21", "CAVRDSNYQLIW", "TRBV20-1", "CSARDRTGNGYTF"))
    donors <- sprintf("D%02d", seq_len(n_donors))
    tissues <- c("blood", "PBMC", "lung", "colon", "skin", "tumor", "lymph node")
    rows <- list()
    for (d in donors) {
      for (b in seq_along(base)) {
        nclones <- sample(3:10, 1)
        for (k in seq_len(nclones)) {
          t <- base[[b]]
          if (k > 1) {
            t[2] <- mut(t[2], sample(0:2, 1))
            t[4] <- mut(t[4], sample(0:3, 1))
          }
          ncell <- sample(1:4, 1)
          rows[[length(rows) + 1]] <- data.frame(
            donor = d, va = t[1], ca = t[2], vb = t[3], cb = t[4],
            tissue = sample(tissues, ncell, replace = TRUE),
            stringsAsFactors = FALSE)
        }
      }
    }
    df <- do.call(rbind, rows)
    df$cell <- sprintf("cell%04d", seq_len(nrow(df)))
    dinfo <- data.frame(donor = donors,
                        disease = rep_len(c("healthy", "disease X", "disease Y"), n_donors),
                        hla = rep_len(c("A*02:01;A*24:02;B*08:01", "A*01:01;A*02:01;B*07:02",
                                        "A*03:01;A*11:01;B*35:01", NA), n_donors))
    df <- merge(df, dinfo, by = "donor")
    # synthetic J genes: mostly the motif's J, sometimes another
    jmap_a <- c("TRAV12-2" = "TRAJ42", "TRAV27" = "TRAJ42", "TRAV21" = "TRAJ33")
    jmap_b <- c("TRBV19" = "TRBJ2-7", "TRBV6-5" = "TRBJ2-7", "TRBV20-1" = "TRBJ1-2")
    df$ja <- ifelse(stats::runif(nrow(df)) < 0.8, jmap_a[df$va],
                    sample(c("TRAJ9", "TRAJ49", "TRAJ53"), nrow(df), TRUE))
    df$jb <- ifelse(stats::runif(nrow(df)) < 0.8, jmap_b[df$vb],
                    sample(c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"), nrow(df), TRUE))
    cells_from_df(df, map = c(cell = "cell", donor = "donor", va = "va", ca = "ca", vb = "vb",
                              cb = "cb", ja = "ja", jb = "jb"), verbose = FALSE)
  })
}
