test_that("10x, AIRR and scRepertoire fixtures give identical cell tables", {
  meta <- read.csv(fx("meta_synthetic.csv"))
  a <- read_10x_contigs(fx("contigs_10x_synthetic.csv"), meta = meta, donor = patient,
                        verbose = FALSE)
  b <- read_airr(fx("airr_synthetic.tsv"), meta = meta, donor = patient, verbose = FALSE)
  s <- read.delim(fx("screp_synthetic.tsv"), colClasses = "character")
  s <- merge(s[, setdiff(names(s), "seurat_clusters")], meta, by = "barcode")
  c <- read_screpertoire(s, donor = patient, verbose = FALSE)
  expect_identical(core(a), core(b))
  expect_identical(core(a), core(c))
  expect_equal(nrow(a), 6)
  # two TRA in AAAG-1: the one with most UMIs is kept
  expect_equal(a$va[a$cell == "AAAG-1"], "TRAV27")
  # nonproductive TRA -> beta-only cell
  expect_true(is.na(a$ca[a$cell == "AACC-1"]))
  # V normalisation: Cell Ranger "TRAV29DV5", allele suffix
  expect_equal(a$va[a$cell == "AACG-1"], "TRAV29/DV5")
  expect_equal(a$vb[a$cell == "AACG-1"], "TRBV7-2")
  s1 <- attr(a, "tcrhex_summary")
  expect_equal(s1$paired, 4)
  expect_equal(s1$beta_only, 2)
  expect_equal(s1$extra_a, 1)
  expect_equal(s1$extra_b, 1)
  expect_equal(attr(c, "tcrhex_summary")$extra_a, 1)
  # all metadata columns kept
  expect_true("seurat_clusters" %in% names(a))
})

test_that("readers print a summary", {
  expect_output(read_10x_contigs(fx("contigs_10x_synthetic.csv"), donor = "P1"),
                "6 cells, 4 paired")
})

test_that("cells_from_df maps columns and strips alleles", {
  df <- data.frame(bc = c("x", "y", "z"), who = c("D1", "D1", NA),
                   a_v = "TRAV12-2*01", a_c = c("CAVNGGSQGNLIF", "notacdr3", "CAVNGGSQGNLIF"),
                   b_v = "TRBV19*01", b_c = "CASSIRSSYEQYF", site = "blood")
  x <- cells_from_df(df, map = c(cell = "bc", donor = "who", va = "a_v", ca = "a_c",
                                 vb = "b_v", cb = "b_c"), verbose = FALSE)
  expect_equal(x$va, c("TRAV12-2", NA))
  expect_equal(nrow(x), 2)               # z has no donor
  expect_equal(attr(x, "tcrhex_summary")$invalid_cdr3, 1)
  expect_true("site" %in% names(x))
})

test_that("normalize_v handles 10x / AIRR spellings", {
  expect_equal(normalize_v(c("TRBV19*01", "TRAV29DV5", "TRAV38-2DV8", "TRBV6-5*01,TRBV6-6*01",
                             "None", "")),
               c("TRBV19", "TRAV29/DV5", "TRAV38-2/DV8", "TRBV6-5", NA, NA))
})

test_that("scRepertoire objects are read (public example data)", {
  skip_if_not_installed("scRepertoire")
  data("contig_list", package = "scRepertoire", envir = environment())
  comb <- scRepertoire::combineTCR(contig_list[1:2], samples = c("P17B", "P17L"))
  x <- read_screpertoire(comb, donor = "P17", verbose = FALSE)
  expect_s3_class(x, "tcrhex_cells")
  expect_true(all(x$donor == "P17"))
  y <- read_10x_contigs(setNames(contig_list[1:2], c("P17B", "P17L")), donor = "P17",
                        verbose = FALSE)
  # identical on cells with one chain per locus
  single <- unlist(lapply(comb, function(d) d$barcode[!grepl(";", d$CTaa)]))
  a <- as.data.frame(x[x$cell %in% single, 1:6])
  b <- as.data.frame(y[match(a$cell, y$cell), 1:6])
  rownames(a) <- rownames(b) <- NULL
  attr(a, "tcrhex_summary") <- attr(b, "tcrhex_summary") <- NULL
  expect_identical(a, b)
})
