# Figure regression against the HPC Python output. The bundle holds research
# data and is never shipped: put it in data-raw/bundle/ (untracked) or point
# TCRHEX_BUNDLE at it. Skipped when absent.
bundle_dir <- function() {
  d <- Sys.getenv("TCRHEX_BUNDLE", "")
  if (d == "") d <- file.path(test_path(), "..", "..", "data-raw", "bundle")
  if (dir.exists(file.path(d, "hits"))) normalizePath(d) else NULL
}

read_bundle <- function(d) {
  q <- utils::read.delim(file.path(d, "queries", "top20_queries_CD8mem.tsv"),
                         colClasses = "character")
  q$rank <- as.integer(q$rank)
  q$concordance <- as.numeric(q$concordance)
  files <- sort(list.files(file.path(d, "hits"), "^rank[0-9]+_hits\\.tsv$", full.names = TRUE))
  hits <- do.call(rbind, lapply(seq_along(files), function(i) {
    h <- utils::read.delim(files[i], comment.char = "#", colClasses = "character")
    h$query <- q$name[i]
    h
  }))
  for (v in c("d", "d_a", "d_b", "cells")) hits[[v]] <- as.integer(hits[[v]])
  hits$donor <- hits$person_id
  list(q = q, hits = hits)
}

test_that("bundle: our d_a / d_b equal tcrdist3's on every hit", {
  d <- bundle_dir()
  skip_if(is.null(d), "HPC bundle not present (data-raw/bundle or TCRHEX_BUNDLE)")
  b <- read_bundle(d)
  cl <- unique(b$hits[, c("va", "ca", "vb", "cb")])
  cl$.clone_id <- seq_len(nrow(cl))
  h <- tcrdist_paired(b$q[, c("name", "va", "ca", "vb", "cb")], cl, radius = Inf)
  b$hits$.clone_id <- match(paste(b$hits$va, b$hits$ca, b$hits$vb, b$hits$cb),
                            paste(cl$va, cl$ca, cl$vb, cl$cb))
  m <- merge(b$hits, h, by.x = c("query", ".clone_id"), by.y = c("query", ".clone_id"),
             suffixes = c("", ".r"))
  bad <- m[m$d_a != m$d_a.r | m$d_b != m$d_b.r, ]
  if (nrow(bad)) print(utils::head(bad[, c("query", "va", "ca", "vb", "cb", "d_a", "d_a.r",
                                           "d_b", "d_b.r")], 20))
  expect_equal(nrow(m), nrow(b$hits))
  expect_equal(nrow(bad), 0)
})

compare_table <- function(hd, ref_file) {
  ref <- utils::read.delim(ref_file, colClasses = "character")
  ours <- hd$bins
  ours <- data.frame(row = as.character(ours$row), bin = as.character(ours$bin),
                     level = as.character(ours$level), clones = as.character(ours$clones),
                     shown = as.character(ours$shown))
  ref <- ref[, c("row", "bin", "level", "clones", "shown")]
  key <- function(x) x[order(as.integer(x$row), as.integer(x$bin), x$level), ]
  a <- key(ours); r <- key(ref)
  rownames(a) <- rownames(r) <- NULL
  if (!identical(a, r)) {
    diff <- merge(a, r, by = c("row", "bin", "level"), all = TRUE, suffixes = c(".R", ".py"))
    diff <- diff[is.na(diff$clones.R) | is.na(diff$clones.py) |
                   diff$clones.R != diff$clones.py | diff$shown.R != diff$shown.py, ]
    print(utils::head(diff, 30))
  }
  expect_identical(a, r)
}

test_that("bundle: per-bin tables equal the Python ones (tissue, disease)", {
  d <- bundle_dir()
  skip_if(is.null(d), "HPC bundle not present (data-raw/bundle or TCRHEX_BUNDLE)")
  b <- read_bundle(d)
  hd <- hexrows_data(NULL, ggplot2::aes(fill = tissue_class(tissue_class)), queries = b$q,
                     hits = b$hits, summarise_fill = summarise_tissue, logos = FALSE)
  compare_table(hd, file.path(d, "reference", "top20_rows_CD8mem_r108_tissue_table.tsv"))
  hd <- hexrows_data(NULL, ggplot2::aes(fill = disease), queries = b$q, hits = b$hits,
                     logos = FALSE)
  compare_table(hd, file.path(d, "reference", "top20_rows_CD8mem_r108_disease_table.tsv"))
})
