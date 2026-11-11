# ---------------------------------------------------------------------------
# Readers: every reader returns the same validated cell table
# (cell, donor, va, ca, vb, cb, <metadata ...>) of class "tcrhex_cells".
# ---------------------------------------------------------------------------

.tcr_cols <- c("cell", "donor", "va", "ca", "vb", "cb")

#' Read scRepertoire output into a cell table
#'
#' Accepts the list of data frames returned by `scRepertoire::combineTCR()`, a
#' Seurat or SingleCellExperiment object after `scRepertoire::combineExpression()`,
#' or a plain data frame of that meta data. The `CTgene` and `CTaa` columns are
#' parsed into TRAV / CDR3a / TRBV / CDR3b: chains are separated by `"_"`
#' (alpha first), several chains of one locus by `";"`, and the genes of a chain
#' by `"."` (V.J.C for alpha, V.D.J.C for beta); a missing chain is `"NA"`.
#'
#' scRepertoire does not keep UMI counts in `CTgene` / `CTaa` and lists
#' same-locus chains alphabetically by nucleotide CDR3, so for a cell with two
#' alpha (or beta) chains the first listed chain is used and the cell is
#' counted in the summary (run `combineTCR(filterMulti = TRUE)` to let
#' scRepertoire keep the most-read chain instead).
#'
#' @param x A `combineTCR()` list, a Seurat object, a SingleCellExperiment, or
#'   a data frame with `CTgene` and `CTaa` columns (barcodes in a `barcode`
#'   column or the row names).
#' @param donor The donor (person) of each cell: an unquoted expression
#'   evaluated in the cell data (e.g. `donor = patient` or
#'   `donor = sub("[BL]$", "", sample)`), or a single string for one donor.
#' @param ... Further metadata columns to create, as in [dplyr::mutate()],
#'   e.g. `tissue = ifelse(grepl("B$", sample), "blood", "lung")`.
#' @param verbose Print the summary.
#' @return A cell table (tibble of class `tcrhex_cells`): one row per cell with
#'   at least one productive TRA or TRB chain, columns `cell`, `donor`, `va`,
#'   `ca`, `vb`, `cb`, the J genes `ja`, `jb` when the input has them, and all
#'   metadata columns.
#' @seealso [read_10x_contigs()], [read_airr()], [cells_from_df()]
#' @export
#' @examples
#' if (requireNamespace("scRepertoire", quietly = TRUE)) {
#'   data("contig_list", package = "scRepertoire")
#'   combined <- scRepertoire::combineTCR(contig_list,
#'     samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
#'   cells <- read_screpertoire(combined, donor = sub("[BL]$", "", sample),
#'                              tissue = ifelse(grepl("B$", sample), "B", "L"))
#' }
read_screpertoire <- function(x, donor, ..., verbose = TRUE) {
  df <- .screp_meta(x)
  if (!all(c("CTgene", "CTaa") %in% names(df))) {
    stop("`x` has no CTgene / CTaa columns: is it scRepertoire output?", call. = FALSE)
  }
  df <- df[!is.na(df$CTgene) | !is.na(df$CTaa), , drop = FALSE]
  g <- .split2(df$CTgene)
  a <- .split2(df$CTaa)
  stats <- list(input = nrow(df))
  pick <- function(genes, aas, locus, jpos) {
    gl <- strsplit(ifelse(is.na(genes), "NA", genes), ";", fixed = TRUE)
    al <- strsplit(ifelse(is.na(aas), "NA", aas), ";", fixed = TRUE)
    n <- pmax(lengths(gl), lengths(al))
    stats[[paste0("extra_", locus)]] <<- sum(n > 1 & !is.na(aas) & aas != "NA")
    v <- vapply(gl, function(z) sub("\\..*$", "", z[1]), "")
    j <- vapply(gl, function(z) {
      f <- strsplit(z[1], ".", fixed = TRUE)[[1]]
      if (length(f) >= jpos) f[jpos] else NA_character_
    }, "")
    cd <- vapply(al, function(z) z[1], "")
    v[v %in% c("NA", "")] <- NA
    j[j %in% c("NA", "", "None")] <- NA
    cd[cd %in% c("NA", "")] <- NA
    list(v = v, j = j, cdr3 = cd)
  }
  pa <- pick(g[, 1], a[, 1], "a", 2)   # alpha: V.J.C
  pb <- pick(g[, 2], a[, 2], "b", 3)   # beta:  V.D.J.C
  meta <- df[, setdiff(names(df), c("barcode", "CTgene", "CTaa", "CTnt", "CTstrict",
                                    "TCR1", "TCR2", "cdr3_aa1", "cdr3_aa2",
                                    "cdr3_nt1", "cdr3_nt2")), drop = FALSE]
  cells <- tibble::tibble(cell = df$barcode, va = pa$v, ca = pa$cdr3,
                          vb = pb$v, cb = pb$cdr3, ja = pa$j, jb = pb$j)
  cells <- dplyr::bind_cols(cells, tibble::as_tibble(meta))
  .finish_cells(cells, rlang::enquo(donor), rlang::enquos(...), stats, verbose,
                source = "scRepertoire")
}

# Extract a meta data frame (with a barcode column) from scRepertoire outputs.
.screp_meta <- function(x) {
  if (inherits(x, "Seurat")) {
    df <- x[[]]
    df$barcode <- rownames(df)
  } else if (inherits(x, "SingleCellExperiment") || inherits(x, "SummarizedExperiment")) {
    df <- as.data.frame(SummarizedExperiment::colData(x))
    df$barcode <- rownames(df)
  } else if (is.data.frame(x)) {
    df <- as.data.frame(x)
    if (!"barcode" %in% names(df)) df$barcode <- rownames(df)
  } else if (is.list(x) && all(vapply(x, is.data.frame, TRUE))) {
    nm <- names(x)
    if (is.null(nm)) nm <- as.character(seq_along(x))
    df <- dplyr::bind_rows(lapply(seq_along(x), function(i) {
      d <- as.data.frame(x[[i]])
      if (!"sample" %in% names(d)) d$sample <- nm[i]
      d
    }))
  } else {
    stop("`x` must be a combineTCR() list, a Seurat object, a SingleCellExperiment ",
         "or a data frame.", call. = FALSE)
  }
  df$barcode <- as.character(df$barcode)
  df
}

# Split "alpha_beta" strings into a 2-column matrix ("NA" stays "NA").
.split2 <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA_NA"
  i <- regexpr("_", x, fixed = TRUE)
  a <- ifelse(i > 0, substr(x, 1, i - 1), x)
  b <- ifelse(i > 0, substr(x, i + 1, nchar(x)), "NA")
  cbind(a, b)
}

#' Read Cell Ranger contig annotations into a cell table
#'
#' Reads one or several 10x Cell Ranger `filtered_contig_annotations.csv` files
#' (or data frames with the same columns: `barcode`, `chain`, `v_gene`, `cdr3`,
#' `productive`, `high_confidence`, `umis`, `reads`) and keeps productive,
#' high-confidence TRA / TRB contigs. When a cell has two or more chains of one
#' locus, the one with most UMIs (then reads, then alphabetical) is kept and the
#' cell is counted in the summary.
#'
#' @param path Path(s) to the CSV file(s), or a data frame / named list of data
#'   frames. Names (of `path` or of the list) become a `sample` column; with
#'   several samples, barcodes are prefixed `"<sample>_"`, as scRepertoire does.
#' @param meta Optional metadata data frame: one row per cell (joined by its
#'   barcode column) or one row per sample (`by = "sample"`), or both
#'   (`by = c("barcode", "sample")`).
#' @param by Join columns of `meta`. A `"barcode"` (or `"cell"`, `"cell_id"`)
#'   column is matched against the raw barcode or the `"<sample>_"`-prefixed
#'   one; `"sample"` against the sample names.
#' @inheritParams read_screpertoire
#' @inherit read_screpertoire return
#' @export
#' @examples
#' if (requireNamespace("scRepertoire", quietly = TRUE)) {
#'   data("contig_list", package = "scRepertoire")
#'   contigs <- setNames(contig_list[1:2], c("P17B", "P17L"))
#'   cells <- read_10x_contigs(contigs, donor = "P17",
#'                             tissue = ifelse(sample == "P17B", "B", "L"))
#' }
read_10x_contigs <- function(path, meta = NULL, by = "barcode", donor, ...,
                             verbose = TRUE) {
  tabs <- .read_tables(path, reader = function(f) {
    readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                    col_types = readr::cols(.default = readr::col_character()))
  })
  chains <- dplyr::bind_rows(lapply(names(tabs), function(s) {
    d <- tabs[[s]]
    need <- c("barcode", "chain", "v_gene", "cdr3")
    miss <- setdiff(need, names(d))
    if (length(miss)) stop("contig table '", s, "' lacks column(s): ",
                           paste(miss, collapse = ", "), call. = FALSE)
    tibble::tibble(
      sample = s, barcode = as.character(d$barcode), locus = as.character(d$chain),
      v = d$v_gene, j = if ("j_gene" %in% names(d)) d$j_gene else NA_character_,
      cdr3 = d$cdr3,
      productive = if ("productive" %in% names(d)) .as_flag(d$productive) else TRUE,
      high_confidence = if ("high_confidence" %in% names(d)) .as_flag(d$high_confidence) else TRUE,
      umis = if ("umis" %in% names(d)) suppressWarnings(as.numeric(d$umis)) else NA_real_,
      reads = if ("reads" %in% names(d)) suppressWarnings(as.numeric(d$reads)) else NA_real_)
  }))
  .contigs_to_cells(chains, meta, by, rlang::enquo(donor), rlang::enquos(...),
                    verbose, source = "10x contigs", multi = length(tabs) > 1)
}

#' Read an AIRR rearrangement table into a cell table
#'
#' Reads AIRR rearrangement TSV file(s) (columns `cell_id`, `locus`, `v_call`,
#' `junction_aa`, `productive`, and `umi_count` / `duplicate_count` /
#' `consensus_count`). The junction (CDR3 including the conserved C and F/W,
#' as in Cell Ranger's `cdr3`) is used; `cdr3_aa` is used, with the flanking
#' residues re-added only if present, when `junction_aa` is missing. Chain
#' selection is as in [read_10x_contigs()].
#'
#' @inheritParams read_10x_contigs
#' @inherit read_screpertoire return
#' @export
read_airr <- function(path, meta = NULL, by = "barcode", donor, ..., verbose = TRUE) {
  tabs <- .read_tables(path, reader = function(f) {
    readr::read_tsv(f, show_col_types = FALSE, progress = FALSE, comment = "#",
                    col_types = readr::cols(.default = readr::col_character()))
  })
  chains <- dplyr::bind_rows(lapply(names(tabs), function(s) {
    d <- tabs[[s]]
    need <- c("cell_id", "locus", "v_call")
    miss <- setdiff(need, names(d))
    if (length(miss)) stop("AIRR table '", s, "' lacks column(s): ",
                           paste(miss, collapse = ", "), call. = FALSE)
    cdr3 <- if ("junction_aa" %in% names(d)) d$junction_aa else {
      warning("AIRR table has no junction_aa; using cdr3_aa (without the conserved ",
              "C / F flanks, distances will differ from 10x input).", call. = FALSE)
      d$cdr3_aa
    }
    num <- function(col) if (col %in% names(d)) suppressWarnings(as.numeric(d[[col]])) else rep(NA_real_, nrow(d))
    umis <- num("umi_count")
    reads <- num("duplicate_count")
    reads[is.na(reads)] <- num("consensus_count")[is.na(reads)]
    tibble::tibble(sample = s, barcode = as.character(d$cell_id),
                   locus = toupper(as.character(d$locus)), v = d$v_call,
                   j = if ("j_call" %in% names(d)) d$j_call else NA_character_, cdr3 = cdr3,
                   productive = if ("productive" %in% names(d)) .as_flag(d$productive) else TRUE,
                   high_confidence = TRUE, umis = umis, reads = reads)
  }))
  .contigs_to_cells(chains, meta, by, rlang::enquo(donor), rlang::enquos(...),
                    verbose, source = "AIRR", multi = length(tabs) > 1)
}

.as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES")
}

.read_tables <- function(path, reader) {
  if (is.data.frame(path)) return(list(`1` = tibble::as_tibble(path)))
  if (is.list(path)) {
    nm <- names(path)
    if (is.null(nm)) nm <- as.character(seq_along(path))
    out <- lapply(path, function(p) if (is.data.frame(p)) tibble::as_tibble(p) else reader(p))
    return(setNames(out, nm))
  }
  nm <- names(path)
  if (is.null(nm)) nm <- if (length(path) == 1) "1" else basename(dirname(path))
  if (anyDuplicated(nm)) nm <- make.unique(nm)
  setNames(lapply(path, reader), nm)
}

# Long chain table (sample, barcode, locus, v, cdr3, productive,
# high_confidence, umis, reads) -> cell table.
.contigs_to_cells <- function(chains, meta, by, donor_q, dots, verbose, source, multi) {
  stats <- list(input = dplyr::n_distinct(paste(chains$sample, chains$barcode)))
  chains$cell <- if (multi) paste0(chains$sample, "_", chains$barcode) else chains$barcode
  keep <- chains$productive & chains$high_confidence & chains$locus %in% c("TRA", "TRB") &
    !is.na(chains$cdr3) & !chains$cdr3 %in% c("", "None", "NA")
  chains <- chains[keep, , drop = FALSE]
  chains$umis[is.na(chains$umis)] <- -1
  chains$reads[is.na(chains$reads)] <- -1
  chains <- chains[order(chains$cell, chains$locus, -chains$umis, -chains$reads,
                         chains$v, chains$cdr3), , drop = FALSE]
  key <- paste(chains$cell, chains$locus)
  nper <- table(key)
  first <- !duplicated(key)
  multi_key <- names(nper)[nper > 1]
  stats$extra_a <- sum(grepl(" TRA$", multi_key))
  stats$extra_b <- sum(grepl(" TRB$", multi_key))
  ch <- chains[first, , drop = FALSE]
  a <- ch[ch$locus == "TRA", ]
  b <- ch[ch$locus == "TRB", ]
  cells <- dplyr::distinct(ch[, c("cell", "sample", "barcode")])
  cells <- tibble::tibble(cell = cells$cell, va = a$v[match(cells$cell, a$cell)],
                          ca = a$cdr3[match(cells$cell, a$cell)],
                          vb = b$v[match(cells$cell, b$cell)],
                          cb = b$cdr3[match(cells$cell, b$cell)],
                          ja = a$j[match(cells$cell, a$cell)],
                          jb = b$j[match(cells$cell, b$cell)],
                          sample = cells$sample, barcode = cells$barcode)
  if (!multi && all(cells$sample == "1")) cells$sample <- NULL
  if (!is.null(meta)) {
    meta <- tibble::as_tibble(meta)
    miss <- setdiff(by, names(meta))
    if (length(miss)) stop("`meta` has no column(s) ", paste(miss, collapse = ", "), call. = FALSE)
    bcol <- intersect(by, c("barcode", "cell", "cell_id"))[1]
    others <- setdiff(by, bcol)
    jb <- stats::setNames(others, others)
    if (!is.na(bcol)) {
      meta[[bcol]] <- as.character(meta[[bcol]])
      # match on the prefixed cell id or the raw barcode, whichever matches more
      use_cell <- mean(cells$cell %in% meta[[bcol]]) >= mean(cells$barcode %in% meta[[bcol]])
      jb <- c(stats::setNames(bcol, if (use_cell) "cell" else "barcode"), jb)
    }
    for (o in others) meta[[o]] <- as.character(meta[[o]])
    if (any(!others %in% names(cells))) {
      stop("Join column(s) ", paste(setdiff(others, names(cells)), collapse = ", "),
           " not in the contig data (only `sample` is).", call. = FALSE)
    }
    meta <- meta[, c(unname(jb), setdiff(names(meta), c(unname(jb), names(cells)))), drop = FALSE]
    if (anyDuplicated(meta[, unname(jb), drop = FALSE])) {
      stop("`meta` has duplicated keys in ", paste(unname(jb), collapse = " + "), call. = FALSE)
    }
    cells <- dplyr::left_join(cells, meta, by = jb)
  }
  cells$barcode <- NULL
  .finish_cells(cells, donor_q, dots, stats, verbose, source = source)
}

#' Build a cell table from any data frame
#'
#' The generic reader: maps columns of `df` (one row per cell, one TRA and one
#' TRB chain per row) onto the internal names.
#'
#' @param df Data frame, one row per cell.
#' @param map Named character vector `c(internal = "column in df")` for
#'   `cell`, `donor`, `va`, `ca`, `vb`, `cb`, and optionally the J genes `ja`,
#'   `jb` (used by `hexrows(gene_logos = TRUE)`). Names not given default to
#'   themselves; `cell` defaults to the row number if absent.
#' @inheritParams read_screpertoire
#' @inherit read_screpertoire return
#' @export
#' @examples
#' # synthetic example (not real data)
#' df <- data.frame(barcode = c("c1", "c2", "c3"), patient = c("D1", "D1", "D2"),
#'                  TRAV = c("TRAV12-2*01", "TRAV12-2", "TRAV27"),
#'                  CDR3a = c("CAVNGGSQGNLIF", "CAVNGGSQGNLIF", NA),
#'                  TRBV = "TRBV19", CDR3b = "CASSIRSSYEQYF", tissue = "blood")
#' cells_from_df(df, map = c(cell = "barcode", donor = "patient", va = "TRAV",
#'                           ca = "CDR3a", vb = "TRBV", cb = "CDR3b"))
cells_from_df <- function(df, map = c(cell = "cell", donor = "donor", va = "va",
                                      ca = "ca", vb = "vb", cb = "cb"),
                          ..., verbose = TRUE) {
  df <- tibble::as_tibble(df)
  full <- setNames(.tcr_cols, .tcr_cols)
  full[names(map)] <- map
  if (!full[["cell"]] %in% names(df)) df[[full[["cell"]]]] <- as.character(seq_len(nrow(df)))
  miss <- setdiff(unname(full), names(df))
  if (length(miss)) stop("`df` lacks mapped column(s): ", paste(miss, collapse = ", "),
                         call. = FALSE)
  out <- tibble::tibble(cell = as.character(df[[full[["cell"]]]]),
                        donor = df[[full[["donor"]]]],
                        va = df[[full[["va"]]]], ca = df[[full[["ca"]]]],
                        vb = df[[full[["vb"]]]], cb = df[[full[["cb"]]]])
  jmap <- map[intersect(names(map), c("ja", "jb"))]
  for (jn in names(jmap)) {
    if (!jmap[[jn]] %in% names(df)) stop("`df` lacks mapped column: ", jmap[[jn]], call. = FALSE)
    out[[jn]] <- df[[jmap[[jn]]]]
  }
  rest <- df[, setdiff(names(df), c(unname(full), unname(jmap))), drop = FALSE]
  rest <- rest[, setdiff(names(rest), names(out)), drop = FALSE]
  out <- dplyr::bind_cols(out, rest)
  .finish_cells(out, rlang::quo(donor), rlang::enquos(...),
                list(input = nrow(out)), verbose, source = "data frame")
}

# Common tail: evaluate donor / extra columns, clean and validate chains,
# record the summary, set the class.
.finish_cells <- function(cells, donor_q, dots, stats, verbose, source) {
  cells <- tibble::as_tibble(cells)
  if (!rlang::quo_is_missing(donor_q)) {
    dv <- rlang::eval_tidy(donor_q, data = cells)
    if (length(dv) == 1) dv <- rep(dv, nrow(cells))
    cells$donor <- dv
  }
  if (!"donor" %in% names(cells)) {
    stop("No donor: give `donor =` (a column, an expression or a single string).",
         call. = FALSE)
  }
  if (length(dots)) cells <- dplyr::mutate(cells, !!!dots)
  cells$donor <- as.character(cells$donor)
  cells$va <- normalize_v(cells$va)
  cells$vb <- normalize_v(cells$vb)
  cells$ca <- toupper(trimws(as.character(cells$ca)))
  cells$cb <- toupper(trimws(as.character(cells$cb)))
  cells$ca[cells$ca %in% c("", "NA", "NONE")] <- NA
  cells$cb[cells$cb %in% c("", "NA", "NONE")] <- NA
  bad_a <- !is.na(cells$ca) & !.valid_cdr3(cells$ca)
  bad_b <- !is.na(cells$cb) & !.valid_cdr3(cells$cb)
  stats$invalid_cdr3 <- sum(bad_a | bad_b)
  cells$ca[bad_a] <- NA
  cells$cb[bad_b] <- NA
  cells$va[is.na(cells$ca)] <- NA
  cells$vb[is.na(cells$cb)] <- NA
  for (jn in intersect(c("ja", "jb"), names(cells))) {
    cells[[jn]] <- normalize_v(cells[[jn]])
    cells[[jn]][is.na(cells[[if (jn == "ja") "ca" else "cb"]])] <- NA
  }
  has_a <- !is.na(cells$va) & !is.na(cells$ca)
  has_b <- !is.na(cells$vb) & !is.na(cells$cb)
  no_donor <- is.na(cells$donor) | cells$donor == ""
  stats$no_chain <- sum(!has_a & !has_b)
  stats$no_donor <- sum(no_donor & (has_a | has_b))
  cells <- cells[(has_a | has_b) & !no_donor, , drop = FALSE]
  if (anyDuplicated(cells$cell)) {
    cells$cell <- make.unique(cells$cell, sep = "#")
    stats$duplicated_cell_ids <- TRUE
  }
  first <- c(.tcr_cols, intersect(c("ja", "jb"), names(cells)))
  cells <- cells[, c(first, setdiff(names(cells), first))]
  .new_cells(cells, stats, source, verbose)
}

.new_cells <- function(cells, stats = list(), source = "", verbose = FALSE) {
  pa <- !is.na(cells$ca) & !is.na(cells$cb)
  stats$cells <- nrow(cells)
  stats$paired <- sum(pa)
  stats$alpha_only <- sum(!is.na(cells$ca) & is.na(cells$cb))
  stats$beta_only <- sum(is.na(cells$ca) & !is.na(cells$cb))
  stats$donors <- dplyr::n_distinct(cells$donor)
  stats$source <- source
  attr(cells, "tcrhex_summary") <- stats
  class(cells) <- c("tcrhex_cells", setdiff(class(cells), "tcrhex_cells"))
  if (isTRUE(verbose)) {
    .print_summary(stats)
  }
  cells
}

.valid_cdr3 <- function(x) grepl("^C[ACDEFGHIKLMNPQRSTVWY]+$", x)

.print_summary <- function(s) {
  cat(sprintf("tcrhex cell table (%s): %s cells, %s paired (TRA+TRB), %s donors\n",
              s$source, .fmt(s$cells), .fmt(s$paired), .fmt(s$donors)))
  cat(sprintf("  unpaired, kept for top_public() but not for TCRdist: %s alpha-only, %s beta-only\n",
              .fmt(s$alpha_only), .fmt(s$beta_only)))
  extra <- c(s$extra_a %||% NA, s$extra_b %||% NA)
  if (!all(is.na(extra))) {
    cat(sprintf("  cells with 2+ chains of one locus (one kept): %s TRA, %s TRB\n",
                .fmt(extra[1]), .fmt(extra[2])))
  }
  drop <- c(`no productive TRA/TRB` = s$no_chain %||% 0,
            `invalid CDR3 (chain removed)` = s$invalid_cdr3 %||% 0,
            `no donor` = s$no_donor %||% 0)
  if (any(drop > 0)) {
    cat("  dropped / cleaned:", paste(sprintf("%s %s", .fmt(drop[drop > 0]),
                                            names(drop)[drop > 0]), collapse = "; "), "\n")
  }
  invisible(s)
}

.fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.tcrhex_cells <- function(x, ...) {
  s <- attr(x, "tcrhex_summary")
  if (!is.null(s)) .print_summary(s)
  y <- x
  class(y) <- setdiff(class(y), "tcrhex_cells")
  print(y, ...)
  invisible(x)
}

#' Summary of a cell table
#' @param object A `tcrhex_cells` table.
#' @param ... Unused.
#' @return The summary list (invisibly), printed.
#' @export
summary.tcrhex_cells <- function(object, ...) {
  s <- attr(object, "tcrhex_summary")
  .print_summary(s)
}

# Coerce any input to a cell table (for hexrows() etc.)
.as_cells <- function(x) {
  if (inherits(x, "tcrhex_cells")) return(x)
  if (is.data.frame(x) && all(.tcr_cols[-1] %in% names(x))) {
    if (!"cell" %in% names(x)) x$cell <- as.character(seq_len(nrow(x)))
    return(.finish_cells(x, rlang::quo(donor), list(), list(input = nrow(x)),
                         verbose = FALSE, source = "data frame"))
  }
  stop("Expected a cell table from read_screpertoire(), read_10x_contigs(), ",
       "read_airr() or cells_from_df().", call. = FALSE)
}
