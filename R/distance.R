#' Paired TCRdist from query TCRs to clones
#'
#' Paired TCRdist (alpha + beta) with the tcrdist3 human defaults. Per chain:
#' CDR1 + CDR2 + CDR2.5 (pMHC-facing loop) of the V gene, weight 1, compared
#' position by position on the IMGT alignment, plus the CDR3 with weight 3,
#' trimming 3 residues at the N- and 2 at the C-terminus. Mismatches score
#' `min(4, 4 - BLOSUM62(a, b))`; a length difference costs 4 per residue
#' (times the weight), with the gap placed where it gives the smallest
#' distance (tcrdist3 `fixed_gappos = False`). Distances match tcrdist3's
#' `TCRrep(organism = "human", chains = c("alpha", "beta"))` exactly (see the
#' package tests).
#'
#' Only query x clone distances are computed. Loop distances are
#' precomputed per V gene, CDR3 distances per unique CDR3 in C++.
#'
#' @param queries Query TCRs: a character vector
#'   `"TRAV:CDR3a:TRBV:CDR3b"` (optionally named) or a data frame with
#'   columns `va`, `ca`, `vb`, `cb` (and optionally `name`); see
#'   [as_queries()].
#' @param clones Clone table from [tcr_clones()] (or any data frame with `va`,
#'   `ca`, `vb`, `cb`).
#' @param radius Keep clones with `d <= radius`; `Inf` keeps all.
#' @param organism `"human"` (default) or `"mouse"`: which tcrdist3 V-gene
#'   loop table to use. Mouse distances match tcrdist3's
#'   `TCRrep(organism = "mouse")`.
#' @param verbose Report clones dropped for unknown V genes.
#' @return A tibble of hits, one row per (query, clone) within the radius:
#'   `query` (name), `.clone_id`, `d`, `d_a`, `d_b`, then the clone's `donor`,
#'   `va`, `ca`, `vb`, `cb`, `cells`.
#' @export
#' @examples
#' clones <- tcr_clones(tcrhex_example_cells())
#' tcrdist_paired("TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF", clones, radius = 50)
tcrdist_paired <- function(queries, clones, radius = 108, organism = c("human", "mouse"),
                           verbose = TRUE) {
  organism <- match.arg(organism)
  q <- as_queries(queries, organism = organism)
  if (!".clone_id" %in% names(clones)) clones$.clone_id <- seq_len(nrow(clones))
  ok_a <- .v_known(clones$va, organism)
  ok_b <- .v_known(clones$vb, organism)
  ok <- ok_a & ok_b & !is.na(clones$ca) & !is.na(clones$cb)
  if (any(!ok) && verbose) {
    unk <- sort(unique(c(clones$va[!ok_a], clones$vb[!ok_b])))
    unk <- unk[!is.na(unk)]
    guess <- .guess_organism(c(clones$va, clones$vb))
    hint <- if (!is.na(guess) && guess != organism)
      sprintf(" Most V genes are %s genes: use organism = \"%s\"?", guess, guess) else ""
    warning(sprintf(
      "%s of %s clones dropped: V gene not in the tcrdist3 %s loop table (%s).%s",
      .fmt(sum(!ok)), .fmt(nrow(clones)), organism,
      paste(utils::head(unk, 10), collapse = ", "), hint), call. = FALSE)
  }
  cl <- clones[ok, , drop = FALSE]
  dm <- .tcr_dm()
  ua <- unique(cl$ca)
  ub <- unique(cl$cb)
  ia <- match(cl$ca, ua)
  ib <- match(cl$cb, ub)
  ea <- .encode_seqs(ua)
  eb <- .encode_seqs(ub)
  out <- vector("list", nrow(q))
  for (i in seq_len(nrow(q))) {
    va_d <- .vloop_dist(q$va[i], organism)[cl$va]
    vb_d <- .vloop_dist(q$vb[i], organism)[cl$vb]
    ca_d <- .tcrdist_one_to_many(.encode_seqs(q$ca[i])$codes, ea$codes, ea$offsets,
                                 ea$lengths, dm, 1L, 4L, 3L, 2L, FALSE)
    cb_d <- .tcrdist_one_to_many(.encode_seqs(q$cb[i])$codes, eb$codes, eb$offsets,
                                 eb$lengths, dm, 1L, 4L, 3L, 2L, FALSE)
    d_a <- unname(va_d) + 3L * ca_d[ia]
    d_b <- unname(vb_d) + 3L * cb_d[ib]
    d <- d_a + d_b
    keep <- which(d <= radius)
    out[[i]] <- tibble::tibble(query = q$name[i], .clone_id = cl$.clone_id[keep],
                               d = d[keep], d_a = d_a[keep], d_b = d_b[keep])
  }
  hits <- dplyr::bind_rows(out)
  info <- clones[, intersect(c(".clone_id", "donor", "va", "ca", "vb", "cb", "cells"),
                             names(clones)), drop = FALSE]
  hits <- dplyr::left_join(hits, info, by = ".clone_id")
  attr(hits, "radius") <- radius
  attr(hits, "dropped_unknown_v") <- sum(!ok)
  attr(hits, "organism") <- organism
  hits
}

#' Validate query TCRs
#'
#' Turns query TCRs into a validated query table. V names are normalised and
#' checked against the tcrdist3 loop table; CDR3s must be amino-acid letters
#' starting with C. All problems are reported together before anything is
#' computed.
#'
#' @param queries Character vector `"TRAV:CDR3a:TRBV:CDR3b"` (names, if any,
#'   become query names) or a data frame with `va`, `ca`, `vb`, `cb` and
#'   optionally `name`, `label`, `epitope`, `concordance`, `rank`, ... (kept).
#' @param organism `"human"` or `"mouse"`: V genes are checked against that
#'   organism's tcrdist3 loop table.
#' @return A tibble with `name`, `va`, `ca`, `vb`, `cb` and any extra columns.
#' @export
#' @examples
#' as_queries(c(flu = "TRAV27:CAGAGSQGNLIF:TRBV19:CASSIRSSYEQYF"))
as_queries <- function(queries, organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  if (is.character(queries)) {
    parts <- strsplit(queries, ":", fixed = TRUE)
    if (any(lengths(parts) != 4)) {
      bad <- queries[lengths(parts) != 4]
      stop("Queries must look like 'TRAV27:CAGAGSQGNLIF:TRBV19:CASSIRSSYEQYF'; got: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
    parts <- unname(parts)
    q <- tibble::tibble(va = vapply(parts, `[`, "", 1), ca = vapply(parts, `[`, "", 2),
                        vb = vapply(parts, `[`, "", 3), cb = vapply(parts, `[`, "", 4))
    q$name <- if (!is.null(names(queries)) && all(names(queries) != "")) names(queries) else queries
  } else if (is.data.frame(queries)) {
    q <- tibble::as_tibble(queries)
    miss <- setdiff(c("va", "ca", "vb", "cb"), names(q))
    if (length(miss)) stop("Query table lacks column(s): ", paste(miss, collapse = ", "),
                           call. = FALSE)
    if (!"name" %in% names(q)) q$name <- paste(q$va, q$ca, q$vb, q$cb, sep = ":")
  } else {
    stop("`queries` must be a character vector or a data frame.", call. = FALSE)
  }
  q$va <- normalize_v(q$va)
  q$vb <- normalize_v(q$vb)
  q$ca <- toupper(trimws(q$ca))
  q$cb <- toupper(trimws(q$cb))
  q$name <- as.character(q$name)
  probs <- character()
  ba <- which(!.v_known(q$va, organism))
  bb <- which(!.v_known(q$vb, organism))
  if (length(ba)) probs <- c(probs, sprintf("query %s: unknown %s TRAV '%s'", q$name[ba], organism, q$va[ba]))
  if (length(bb)) probs <- c(probs, sprintf("query %s: unknown %s TRBV '%s'", q$name[bb], organism, q$vb[bb]))
  ca <- which(is.na(q$ca) | !.valid_cdr3(q$ca))
  cb <- which(is.na(q$cb) | !.valid_cdr3(q$cb))
  if (length(ca)) probs <- c(probs, sprintf("query %s: CDR3a '%s' is not an amino-acid sequence starting with C", q$name[ca], q$ca[ca]))
  if (length(cb)) probs <- c(probs, sprintf("query %s: CDR3b '%s' is not an amino-acid sequence starting with C", q$name[cb], q$cb[cb]))
  if (anyDuplicated(q$name)) probs <- c(probs, "query names must be unique")
  if (length(c(ba, bb))) {
    guess <- .guess_organism(c(q$va, q$vb))
    if (!is.na(guess) && guess != organism) {
      probs <- c(probs, sprintf("(these look like %s V genes: use organism = \"%s\")", guess, guess))
    }
  }
  if (length(probs)) stop("Invalid queries:\n  ", paste(probs, collapse = "\n  "), call. = FALSE)
  q[, c("name", "va", "ca", "vb", "cb", setdiff(names(q), c("name", "va", "ca", "vb", "cb")))]
}

#' The most public TCRs as queries
#'
#' Finds the `n` beta chains (TRBV + CDR3b) found in the most donors and pairs
#' each with its dominant alpha. Per donor, the alpha (TRAV + CDR3a) with the
#' most cells among that donor's paired cells with the beta (ties
#' alphabetical); the dominant alpha is the most common per-donor alpha
#' (ties: more cells, then alphabetical). `concordance` = donors whose
#' per-donor alpha is the dominant one / donors with any paired alpha for that
#' beta. Betas without any paired alpha cannot be queries and are skipped.
#' Ties in donor counts are broken by more cells, then alphabetically.
#'
#' Beta-only cells count towards `donors` (a beta is public in a donor even
#' without a captured alpha), not towards `paired_donors`.
#'
#' @param cells A cell table.
#' @param n Number of queries.
#' @param min_donors Minimum number of donors for a beta.
#' @param organism `"human"` or `"mouse"` (V genes must be in its loop table).
#' @return A query table (see [as_queries()]) with `rank`, `donors`,
#'   `paired_donors`, `donors_with_alpha`, `concordance`.
#' @export
#' @examples
#' top_public(tcrhex_example_cells(), n = 3)
top_public <- function(cells, n = 20, min_donors = 1, organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  cells <- .as_cells(cells)
  b <- cells[!is.na(cells$cb) & .v_known(cells$vb, organism), , drop = FALSE]
  b$beta <- paste(b$vb, b$cb, sep = ":")
  per <- dplyr::summarise(dplyr::group_by(b, .data$beta),
                          donors = dplyr::n_distinct(.data$donor), n_cells = dplyr::n(),
                          .groups = "drop")
  per <- per[per$donors >= min_donors, ]
  per <- per[order(-per$donors, -per$n_cells, per$beta, method = "radix"), ]
  p <- b[!is.na(b$ca) & .v_known(b$va, organism) & b$beta %in% per$beta, , drop = FALSE]
  p$alpha <- paste(p$va, p$ca, sep = ":")
  # per donor: the alpha with the most cells (ties alphabetical)
  da <- dplyr::count(p, .data$beta, .data$donor, .data$alpha, name = "k")
  da <- da[order(da$beta, da$donor, -da$k, da$alpha, method = "radix"), ]
  da <- da[!duplicated(da[, c("beta", "donor")]), ]
  tot <- dplyr::count(p, .data$beta, .data$alpha, name = "alpha_cells")
  dom <- dplyr::count(da, .data$beta, .data$alpha, name = "donors_with_alpha")
  dom <- dplyr::left_join(dom, tot, by = c("beta", "alpha"))
  dom <- dom[order(dom$beta, -dom$donors_with_alpha, -dom$alpha_cells, dom$alpha,
                   method = "radix"), ]
  dom <- dom[!duplicated(dom$beta), ]
  pd <- dplyr::count(da, .data$beta, name = "paired_donors")
  res <- dplyr::left_join(per, dom, by = "beta")
  res <- dplyr::left_join(res, pd, by = "beta")
  res <- res[!is.na(res$alpha), , drop = FALSE]
  res <- utils::head(res, n)
  if (nrow(res) < n) {
    message(sprintf("top_public(): only %d beta chains with a paired alpha%s.", nrow(res),
                    if (min_donors > 1) sprintf(" in >= %d donors", min_donors) else ""))
  }
  sa <- strsplit(res$alpha, ":", fixed = TRUE)
  sb <- strsplit(res$beta, ":", fixed = TRUE)
  q <- tibble::tibble(
    name = paste(res$alpha, res$beta, sep = ":"),
    va = vapply(sa, `[`, "", 1), ca = vapply(sa, `[`, "", 2),
    vb = vapply(sb, `[`, "", 1), cb = vapply(sb, `[`, "", 2),
    rank = seq_len(nrow(res)), donors = res$donors, paired_donors = res$paired_donors,
    donors_with_alpha = res$donors_with_alpha,
    concordance = res$donors_with_alpha / res$paired_donors)
  attr(q, "source") <- "top_public"
  attr(q, "organism") <- organism
  q
}

#' The most expanded clones as queries, optionally per group
#'
#' Picks the `n` clones (donor x TRAV + CDR3a x TRBV + CDR3b) with the most
#' cells, overall or within each level of `by` (any column or expression of
#' the cell table, e.g. `by = tissue`, `by = paste(tissue, timepoint)`).
#' Within a group, clones are ranked by cells in that group, then by all of
#' the clone's cells, then alphabetically. The same TCR can be picked in
#' several groups (each gets its own row). Only paired cells with V genes in
#' the organism's loop table count.
#'
#' @param cells A cell table.
#' @param n Clones per group.
#' @param by Optional grouping: an unquoted column or expression evaluated in
#'   `cells`. Cells where it is `NA` are ignored. Groups come in factor-level
#'   order, else alphabetical.
#' @param min_cells Minimum cells (in the group) for a clone to be picked.
#' @inheritParams top_public
#' @return A query table (see [as_queries()]) with `group`, `rank_in_group`,
#'   `cells_in_group`, `cells` (all cells of the clone), `donor` and a `label`
#'   ("top 1 in lung: 35 cells, donor m1"), shown as the row tag by
#'   [hexrows()].
#' @export
#' @examples
#' cells <- tcrhex_example_cells()   # synthetic
#' top_expanded(cells, n = 1, by = tissue)
top_expanded <- function(cells, n = 1, by = NULL, min_cells = 1,
                         organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  cells <- .as_cells(cells)
  by_q <- rlang::enquo(by)
  p <- cells[!is.na(cells$ca) & !is.na(cells$cb) & .v_known(cells$va, organism) &
               .v_known(cells$vb, organism), , drop = FALSE]
  if (rlang::quo_is_null(by_q)) {
    g <- rep("all cells", nrow(p))
    glev <- "all cells"
  } else {
    g <- rlang::eval_tidy(by_q, data = p)
    if (length(g) == 1) g <- rep(g, nrow(p))
    glev <- if (is.factor(g)) levels(g) else sort(unique(as.character(g[!is.na(g)])))
    g <- as.character(g)
  }
  p$.group <- g
  p <- p[!is.na(p$.group), , drop = FALSE]
  tcr <- paste(p$va, p$ca, p$vb, p$cb, sep = ":")
  p$.clone <- paste(p$donor, tcr, sep = "\r")
  tot <- table(p$.clone)
  cnt <- dplyr::count(tibble::tibble(group = p$.group, clone = p$.clone, donor = p$donor,
                                     va = p$va, ca = p$ca, vb = p$vb, cb = p$cb),
                      .data$group, .data$clone, .data$donor, .data$va, .data$ca, .data$vb,
                      .data$cb, name = "cells_in_group")
  cnt$cells <- as.integer(tot[cnt$clone])
  cnt <- cnt[cnt$cells_in_group >= min_cells, , drop = FALSE]
  cnt$tcr <- paste(cnt$va, cnt$ca, cnt$vb, cnt$cb, sep = ":")
  cnt <- cnt[order(match(cnt$group, glev), -cnt$cells_in_group, -cnt$cells, cnt$tcr,
                   cnt$donor, method = "radix"), , drop = FALSE]
  cnt$rank_in_group <- stats::ave(seq_len(nrow(cnt)), cnt$group, FUN = seq_along)
  cnt <- cnt[cnt$rank_in_group <= n, , drop = FALSE]
  q <- tibble::tibble(
    name = if (rlang::quo_is_null(by_q)) sprintf("#%d %s", cnt$rank_in_group, cnt$tcr) else
      sprintf("%s #%d: %s", cnt$group, cnt$rank_in_group, cnt$tcr),
    va = cnt$va, ca = cnt$ca, vb = cnt$vb, cb = cnt$cb,
    group = cnt$group, rank_in_group = cnt$rank_in_group,
    cells_in_group = cnt$cells_in_group, cells = cnt$cells, donor = cnt$donor,
    label = sprintf("top %d in %s: %s cell%s, donor %s", cnt$rank_in_group, cnt$group,
                    .fmt(cnt$cells_in_group), ifelse(cnt$cells_in_group == 1, "", "s"),
                    cnt$donor))
  attr(q, "source") <- "top_expanded"
  attr(q, "organism") <- organism
  q
}
