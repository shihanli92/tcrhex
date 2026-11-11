# ---------------------------------------------------------------------------
# Helpers that return new columns / scales; none of them is a plot mode.
# ---------------------------------------------------------------------------

.tissue_levels <- c("blood", "blood + tissue", "lymphoid / bone marrow", "lung / airway",
                    "gut", "skin", "CNS / CSF", "tumour (labelled)", "other tissue",
                    "other / multiple tissues", "unknown")

#' Default keyword rules of tissue_class()
#'
#' A named list (class -> keywords), applied in order to the lower-cased
#' tissue text; the first class with a matching keyword wins.
#' @return A named list of character vectors.
#' @export
#' @examples
#' names(tissue_rules())
tissue_rules <- function() {
  list(
    "tumour (labelled)" = c("tumour", "tumor", "carcinoma", "melanoma", "metasta",
                            "glioma", "lymphoma"),
    "blood" = c("blood", "pbmc"),
    "CNS / CSF" = c("cerebrospinal", "brain", "leptomen", "cortex", "hippocamp",
                    "substantia nigra", "lobe", "spinal cord", "cerebell", "meninges"),
    "lymphoid / bone marrow" = c("lymph node", "tonsil", "spleen", "thymus", "bone marrow",
                                 "adenoid"),
    "lung / airway" = c("lung", "bronch", "trache", "nasal", "airway", "larynx", "pharynx",
                        "alveol", "pleura", "lavage", "sputum"),
    "gut" = c("ileum", "colon", "duodenum", "stomach", "rectum", "intestin", "jejunum",
              "mesenter", "ileocecal", "colorect", "appendix", "esophag", "oesophag",
              "caecum", "cecum", "gastric"),
    "skin" = c("skin", "dermis", "epidermis", "blister"))
}

#' Classify free-text tissue into compartments
#'
#' Keyword rules (first match wins, case-insensitive): see [tissue_rules()].
#' Empty / missing values and "cell culture" give `"unknown"`; anything else
#' unmatched gives `"other tissue"`. Values that already are class labels
#' are kept, so `tissue_class()` can be applied twice.
#'
#' @param x Character vector of tissue descriptions.
#' @param rules Named list of keyword vectors (class -> keywords), in order.
#' @return A factor with the tissue class levels (blood, blood + tissue,
#'   lymphoid / bone marrow, lung / airway, gut, skin, CNS / CSF, tumour
#'   (labelled), other tissue, other / multiple tissues, unknown).
#' @export
#' @examples
#' tissue_class(c("PBMC", "Lung tumor", "colon mucosa", "CSF", "cell culture", "liver"))
tissue_class <- function(x, rules = tissue_rules()) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  known <- as.character(x) %in% c(.tissue_levels, names(rules))
  out[known] <- as.character(x)[known]
  unk <- !known & (is.na(s) | s == "" | grepl("cell culture", s, fixed = TRUE))
  out[unk] <- "unknown"
  for (cl in names(rules)) {
    todo <- is.na(out)
    if (!any(todo)) break
    hit <- Reduce(`|`, lapply(rules[[cl]], function(k) grepl(k, s[todo], fixed = TRUE)))
    out[which(todo)[hit]] <- cl
  }
  out[is.na(out)] <- "other tissue"
  factor(out, levels = unique(c(.tissue_levels, names(rules))))
}

#' Summarise tissue classes of a clone
#'
#' A `summarise_fill` rule for tissue: the values are classed with
#' [tissue_class()] (no-op for class labels); unknown is ignored; one class
#' stays; blood plus any tissue gives `"blood + tissue"`; several non-blood
#' classes give `"other / multiple tissues"`.
#'
#' @param values Tissue values (or classes) of the cells of one clone.
#' @return A one-element factor with the tissue class levels.
#' @export
#' @examples
#' summarise_tissue(c("PBMC", "lung"))
#' summarise_tissue(c("lung", "colon"))
summarise_tissue <- function(values) {
  cl <- as.character(tissue_class(values))
  lv <- levels(tissue_class(character()))
  cl <- unique(cl[cl != "unknown"])
  res <- if (!length(cl)) "unknown" else if (length(cl) == 1) cl else
    if (any(cl %in% c("blood", "blood + tissue"))) "blood + tissue" else
      "other / multiple tissues"
  factor(res, levels = lv)
}

.tissue_fill <- c("blood" = "#0072B2", "blood + tissue" = "#56B4E9",
                  "lymphoid / bone marrow" = "#CC79A7", "lung / airway" = "#E69F00",
                  "gut" = "#E69F00", "skin" = "#E69F00", "CNS / CSF" = "#E69F00",
                  "tumour (labelled)" = "#D55E00", "other tissue" = "#a8a7a0",
                  "other / multiple tissues" = "#a8a7a0", "unknown" = "#d6d5ce")
.tissue_pattern <- c("blood" = "none", "blood + tissue" = "none",
                     "lymphoid / bone marrow" = "none", "lung / airway" = "none",
                     "gut" = "stripe", "skin" = "crosshatch", "CNS / CSF" = "circle",
                     "tumour (labelled)" = "none", "other tissue" = "none",
                     "other / multiple tissues" = "none", "unknown" = "none")

#' Tissue colour and pattern scales
#'
#' Okabe-Ito based colours (colour-blind checked): blood #0072B2; blood +
#' tissue #56B4E9; lymphoid / bone marrow #CC79A7; lung / airway #E69F00
#' plain; gut #E69F00 + stripes; skin #E69F00 + crosshatch; CNS / CSF #E69F00
#' + dots; tumour (labelled) #D55E00; other (/ multiple) tissues #a8a7a0;
#' unknown #d6d5ce. Map the same expression to `fill` and `pattern` (e.g.
#' `aes(fill = tissue_class(tissue), pattern = tissue_class(tissue))`) and
#' add both scales to get one combined legend. On a [hexrows()] plot the
#' legend shows clone counts.
#'
#' @param ... Passed to [ggplot2::scale_fill_manual()] /
#'   [ggplot2::discrete_scale()].
#' @param counts Show clone counts in the legend (hexrows plots only).
#' @return A scale (added to a plot with `+`).
#' @export
scale_fill_tissue <- function(..., counts = TRUE) {
  structure(list(aesthetic = "fill", args = list(...), counts = counts),
            class = "tcrhex_scale")
}

#' @rdname scale_fill_tissue
#' @export
scale_pattern_tissue <- function(..., counts = TRUE) {
  structure(list(aesthetic = "pattern", args = list(...), counts = counts),
            class = "tcrhex_scale")
}

.make_tissue_scale <- function(aesthetic, args, counts = NULL) {
  lab <- if (is.null(counts)) ggplot2::waiver() else .count_labels(counts)
  if (aesthetic == "fill") {
    do.call(ggplot2::scale_fill_manual,
            c(list(values = .tissue_fill, labels = lab, na.value = .tissue_fill[["unknown"]],
                   drop = TRUE), args))
  } else {
    do.call(ggplot2::scale_discrete_manual,
            c(list(aesthetics = "pattern", values = .tissue_pattern, labels = lab,
                   na.value = "none", drop = TRUE), args))
  }
}

#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.tcrhex_scale <- function(object, plot, ...) {
  hd <- attr(plot, "hexrows")
  counts <- NULL
  if (isTRUE(object$counts) && !is.null(hd)) {
    tab <- if (object$aesthetic == "pattern" && !is.null(hd$pattern_levels))
      hd$pattern_levels else hd$levels
    counts <- stats::setNames(tab$clones, as.character(tab$level))
  }
  sc <- .make_tissue_scale(object$aesthetic, object$args, counts)
  out <- ggplot2::ggplot_add(sc, plot, ...)
  if (!is.null(hd)) attr(out, "hexrows") <- hd
  out
}

# ---- HLA ----

# Parse HLA strings into 1-field alleles, e.g. "A02:01;A*24:02;HLA-B*08:01"
# -> c("A*02", "A*24", "B*08"). Returns a list.
.hla_parse <- function(hla) {
  lapply(as.character(hla), function(s) {
    if (is.na(s) || trimws(s) == "") return(NA_character_)
    tok <- trimws(strsplit(s, "[;,| ]+")[[1]])
    tok <- tok[tok != ""]
    a <- .hla_1field(tok)
    unique(a[!is.na(a)])
  })
}

.hla_1field <- function(x) {
  x <- toupper(sub("^HLA-", "", toupper(trimws(x))))
  m <- regmatches(x, regexec(
    "^(A|B|C|E|F|G|DRB[1-9]|DQA1|DQB1|DPA1|DPB1|DRA)\\*?([0-9]{1,3})", x))
  vapply(m, function(z) if (length(z) == 3) sprintf("%s*%02d", z[2], as.integer(z[3])) else NA_character_, "")
}

#' HLA carrier status
#'
#' `hla` is a donor's HLA typing string (e.g. `"A02:01;A24:02;B*08:01"`;
#' separators `;`, `,`, `|` or space), reduced to 1-field alleles; `allele`
#' is compared at 1-field resolution (`"A*02"`, `"A02"` and `"A*02:01"` are
#' the same). Vectorised over both, so a per-row allele works:
#' `aes(fill = hla_carrier(hla, .row_allele))`.
#'
#' @param hla Character vector of HLA typing strings (NA / empty = not typed).
#' @param allele Allele(s), recycled to `length(hla)`.
#' @return A factor with levels `"carrier"`, `"non-carrier"`, `"not typed"`
#'   (`NA` where `allele` is `NA`).
#' @export
#' @examples
#' hla_carrier(c("A02:01;A24:02;B*08:01", "A*01:01;B*07:02", NA), "A*02")
hla_carrier <- function(hla, allele) {
  n <- length(hla)
  al <- rep_len(.hla_1field(as.character(allele)), n)
  al[is.na(rep_len(allele, n))] <- NA
  ps <- .hla_parse(hla)
  out <- vapply(seq_len(n), function(i) {
    if (is.na(al[i])) return(NA_character_)
    p <- ps[[i]]
    if (length(p) == 1 && is.na(p)) return("not typed")
    if (al[i] %in% p) "carrier" else "non-carrier"
  }, "")
  factor(out, levels = c("carrier", "non-carrier", "not typed"))
}

#' The HLA allele most enriched near each query
#'
#' For each query, takes the typed donors with a clone within `radius`
#' (paired TCRdist) and tests every 1-field allele of `loci` they carry: `k`
#' carriers out of `n` typed donors, against the allele's carrier frequency
#' `p0` among all typed donors of `clones` (or `background`), with a
#' one-sided binomial test `P(X >= k)`. The allele with the smallest p-value
#' among those with `k >= min_carriers` is returned (ties: larger k, then
#' alphabetical).
#'
#' @param hits Output of [tcrdist_paired()] (radius >= `radius`).
#' @param clones Clone table with the donor HLA column.
#' @param radius Distance defining the neighbourhood.
#' @param loci HLA loci considered.
#' @param hla Name of the HLA column.
#' @param background Optional named numeric vector of carrier frequencies
#'   (names like `"A*02"`), replacing the frequencies of `clones`.
#' @param min_carriers Minimum carriers `k`.
#' @param queries Optional query table / names giving the output order (one
#'   row per query, `NA` allele where none qualifies).
#' @return A tibble: `query`, `allele`, `k`, `n`, `p0`, `p`, `text`.
#' @export
#' @examples
#' cells <- tcrhex_example_cells()
#' cl <- tcr_clones(cells)
#' q <- top_public(cells, 3)
#' hla_enriched(tcrdist_paired(q, cl, 24), cl, queries = q)
hla_enriched <- function(hits, clones, radius = 24, loci = c("A", "B", "C"),
                         hla = "hla", background = NULL, min_carriers = 3,
                         queries = NULL) {
  if (!hla %in% names(clones)) stop("`clones` has no column '", hla, "'.", call. = FALSE)
  dh <- clones[!duplicated(clones$donor), c("donor", hla)]
  dh$alleles <- .hla_parse(dh[[hla]])
  typed <- !vapply(dh$alleles, function(a) length(a) == 1 && is.na(a), TRUE)
  dh <- dh[typed, ]
  locus_of <- function(a) sub("\\*.*$", "", a)
  dh$alleles <- lapply(dh$alleles, function(a) a[locus_of(a) %in% loci])
  if (is.null(background)) {
    all_al <- unlist(dh$alleles)
    bg <- table(all_al) / nrow(dh)
    background <- stats::setNames(as.numeric(bg), names(bg))
  }
  h <- hits[hits$d <= radius, , drop = FALSE]
  qn <- if (is.null(queries)) unique(hits$query) else
    if (is.data.frame(queries)) {
      if ("name" %in% names(queries)) queries$name else paste(queries$va, queries$ca, queries$vb, queries$cb, sep = ":")
    } else if (!is.null(names(queries)) && all(names(queries) != "")) names(queries) else queries
  res <- lapply(qn, function(qq) {
    don <- unique(h$donor[h$query == qq])
    al <- dh$alleles[dh$donor %in% don]
    n <- length(al)
    none <- tibble::tibble(query = qq, allele = NA_character_, k = NA_integer_, n = n,
                           p0 = NA_real_, p = NA_real_, text = NA_character_)
    if (n == 0) return(none)
    k <- table(unlist(al))
    k <- k[k >= min_carriers & names(k) %in% names(background)]
    if (!length(k)) return(none)
    p0 <- background[names(k)]
    p <- stats::pbinom(as.integer(k) - 1, n, p0, lower.tail = FALSE)
    o <- order(p, -as.integer(k), names(k))[1]
    a <- names(k)[o]
    ko <- as.integer(k[o])
    p0o <- unname(p0[o])
    po <- unname(p[o])
    tibble::tibble(query = qq, allele = a, k = ko, n = n, p0 = p0o, p = po,
                   text = sprintf("%s in %d/%d typed donors within %s (background %s%%), p = %s",
                                  a, ko, n, radius, round(100 * p0o), format(signif(po, 2))))
  })
  dplyr::bind_rows(res)
}
