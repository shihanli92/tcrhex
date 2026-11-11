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
  # mouse spellings
  expect_equal(normalize_v(c("TRAV15D-1-DV6D-1", "TRAV15D-1/DV6D-1*01", "TRAV16DDV11",
                             "TRAV13-4DV7", "TRAV14D-1", "TRAV6N-7")),
               c("TRAV15D-1/DV6D-1", "TRAV15D-1/DV6D-1", "TRAV16D/DV11", "TRAV13-4/DV7",
                 "TRAV14D-1", "TRAV6N-7"))
  expect_true(all(tcrhex:::.v_known(c("TRAV15D-1/DV6D-1", "TRAV16D/DV11", "TRAV13-4/DV7"),
                                    "mouse")))
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

test_that("Seurat meta.data with AIRR-style _vj / _vdj columns (synthetic mouse cells)", {
  md <- data.frame(
    barcode = c("AAAC-1", "AAAG-1", "AAAT-1", "AACC-1"),
    mouse_id = c("m1", "m1", "m2", "m2"),
    v_call_vj = c(NA, "TRAV6N-7", "TRAV14D-1", "TRAV19"),
    junction_aa_vj = c(NA, "CALGPGTGGYKVVF", "CAARMSNYNVLYF", "CAAGSQGTGSKLSF"),
    v_call_vdj = c(NA, "TRBV1", "TRBV20", NA),
    junction_aa_vdj = c(NA, "CTCRTGTGTYEQYF", "CGARDRGKSLEEVFF", NA),
    umis_vj = c(NA, "5", "2", "5"), cluster = c("0", "1", "1", "2"))
  md <- dplyr::group_by(md, v_call_vj, junction_aa_vj)   # grouped, as after group_by()
  x <- cells_from_df(md, map = c(cell = "barcode", donor = "mouse_id", va = "v_call_vj",
                                 ca = "junction_aa_vj", vb = "v_call_vdj",
                                 cb = "junction_aa_vdj"), verbose = FALSE)
  s <- attr(x, "tcrhex_summary")
  expect_equal(c(s$cells, s$paired, s$alpha_only), c(3L, 2L, 1L))
  expect_equal(attr(x, "tcrhex_summary")$no_chain, 1L)
  expect_true(all(c("umis_vj", "cluster") %in% names(x)))
  expect_false(dplyr::is_grouped_df(x))
  cl <- tcr_clones(x)
  q <- top_public(x, 5, organism = "mouse")
  expect_equal(attr(q, "organism"), "mouse")
  hd <- hexrows_data(x, ggplot2::aes(fill = cluster), queries = q, logos = FALSE)
  expect_equal(sum(hd$bins$clones[hd$bins$bin == 0]), 2)
})
