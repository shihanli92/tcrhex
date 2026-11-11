#' V-gene loop table used by paired TCRdist
#'
#' Human TRAV / TRBV genes with their IMGT-aligned CDR1, CDR2 and CDR2.5
#' (pMHC-facing) loops, derived from tcrdist3's `alphabeta_gammadelta_db.tsv`
#' (MIT License, Copyright (c) 2020 Koshlan Mayer-Blackwell). Genes are keyed
#' without allele; the `*01` allele's loops are used (as when tcrdist3 is given
#' `"TRAV27*01"`). `other_alleles_differ` flags genes whose other alleles have
#' different loops.
#'
#' @return A tibble with columns `gene`, `chain`, `allele`, `cdr1`, `cdr2`,
#'   `cdr25`, `other_alleles_differ`.
#' @export
#' @examples
#' head(vgene_loops())
vgene_loops <- function() {
  if (is.null(.tcrhex_cache$vloops)) {
    f <- system.file("extdata", "tcrdist_vloops_human.tsv", package = "tcrhex")
    .tcrhex_cache$vloops <- readr::read_tsv(f, col_types = "cccccccl", progress = FALSE)
  }
  .tcrhex_cache$vloops
}

.tcrhex_cache <- new.env(parent = emptyenv())

#' Normalise V gene names
#'
#' Strips allele suffixes (`"TRBV19*01"` -> `"TRBV19"`), takes the first call
#' of a comma-separated list, removes an `"HLA-"`-style prefix and restores the
#' IMGT slash in the TRAV/DV genes as written by Cell Ranger
#' (`"TRAV29DV5"` -> `"TRAV29/DV5"`). Empty strings, `"None"` and `"NA"` become
#' `NA`.
#'
#' @param x Character vector of V gene calls.
#' @return Character vector of IMGT gene names.
#' @export
#' @examples
#' normalize_v(c("TRBV19*01", "TRAV29DV5", "TRBV6-5*01,TRBV6-6*01", "None"))
normalize_v <- function(x) {
  x <- as.character(x)
  x <- trimws(sub(",.*$", "", x))
  x <- sub("\\*.*$", "", x)
  x <- toupper(x)
  x <- sub("^([A-Z]+[0-9]+(-[0-9]+)?)(DV[0-9]+)$", "\\1/\\3", x)
  x <- sub("^([A-Z]+[0-9]+(-[0-9]+)?)(OR[0-9]+-[0-9]+)$", "\\1/\\3", x)
  x[x %in% c("", "NONE", "NA")] <- NA_character_
  x
}

# 21 x 21 tcrdist substitution distance matrix: min(4, 4 - BLOSUM62) for the 20
# amino acids in parasail order, 0 on the diagonal; row/col 21 = any other
# symbol (gap '.', X, ...), distance 0 to everything (as in pwseqdist).
.tcr_dm <- function() {
  if (!is.null(.tcrhex_cache$dm)) return(.tcrhex_cache$dm)
  aa <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
  # BLOSUM62 for the 20 amino acids (ARNDCQEGHILKMFPSTWYV order)
  b62 <- c(
    4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,
   -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,
   -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,
   -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,
    0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,
   -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,
   -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,
    0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,
   -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,
   -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,
   -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,
   -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,
   -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,
   -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,
   -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,
    1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,
    0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,
   -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,
   -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,
    0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4)
  B <- matrix(b62, 20, 20, byrow = TRUE, dimnames = list(aa, aa))
  D <- 4L - B
  D[D > 4] <- 4
  diag(D) <- 0
  dm <- matrix(0L, 21, 21)
  dm[1:20, 1:20] <- as.integer(D)
  storage.mode(dm) <- "integer"
  .tcrhex_cache$dm <- dm
  dm
}

# Distance from query V gene to every gene of the same chain, via the loops
# (weight 1 each, no trimming, fixed gap position) -> named integer vector.
.vloop_dist <- function(qv, chain) {
  key <- paste0("vd_", chain, "_", qv)
  if (!is.null(.tcrhex_cache[[key]])) return(.tcrhex_cache[[key]])
  vt <- vgene_loops()
  vt <- vt[vt$chain == chain, ]
  qi <- match(qv, vt$gene)
  dm <- .tcr_dm()
  tot <- integer(nrow(vt))
  for (loop in c("cdr1", "cdr2", "cdr25")) {
    enc <- .encode_seqs(vt[[loop]])
    q <- .encode_seqs(vt[[loop]][qi])$codes
    tot <- tot + .tcrdist_one_to_many(q, enc$codes, enc$offsets, enc$lengths, dm,
                                      1L, 4L, 0L, 0L, TRUE)
  }
  out <- setNames(tot, vt$gene)
  .tcrhex_cache[[key]] <- out
  out
}

# Is a V gene in the loop table for the given chain ("TRA"/"TRB")?
.v_known <- function(v, chain) {
  vt <- vgene_loops()
  v %in% vt$gene[vt$chain == chain]
}
