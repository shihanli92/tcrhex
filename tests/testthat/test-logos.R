test_that("logo alignment: modal, longer and shorter sequences", {
  s <- c("CASSLGQAYEQYF", "CASSLGQAYEQYF", "CASSLGAYEQYF", "CASSLGQGAYEQYF")
  a <- logo_align(s, trim = 1)
  expect_true(all(nchar(a) == 11))                 # modal 13 - 2
  expect_equal(a[1], "ASSLGQAYEQY")
  # shorter (10): gap after ceil(10/2) = 5
  expect_equal(a[3], "ASSLG-AYEQY")
  # longer (12): keep first ceil(11/2) = 6 and last floor(11/2) = 5
  expect_equal(a[4], "ASSLGQAYEQY")
  f <- logo_freqs(a)
  col6 <- f[f$pos == 6, ]
  expect_equal(sum(col6$freq), 0.75)                # gap column is shorter
  expect_equal(col6$letter[nrow(col6)], "Q")        # most frequent on top
})

test_that("logo polygons are built from glyph outlines", {
  g <- tcrhex:::.glyphs()
  expect_true(all(c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "-") %in% g$letter))
  expect_true(all(c(LETTERS, 0:9, "-", "/", "|") %in% g$letter))
  expect_true(all(g$x >= 0 & g$x <= 1 & g$y >= 0 & g$y <= 1))
  p <- tcrhex:::.logo_polys(logo_freqs(logo_align(c("CASSF", "CASSF"), 1)), c(0, 0), c(3, 0))
  # edge along x: letters stand on the baseline 0.35 above, full height 2.6
  expect_equal(min(p$y), 0.35, tolerance = 1e-6)
  expect_equal(max(p$y), 0.35 + 2.6, tolerance = 1e-6)
  expect_equal(range(p$x), c(0.04, 3 - 0.04), tolerance = 1e-6)
})

test_that("gene_words: V and J counted independently, threshold and cap", {
  w <- gene_words(c(rep("TRAV12-2", 6), rep("TRAV21", 3), "TRAV27"), min_freq = 0.2)
  expect_equal(w$word, c("AV12-2", "AV21"))
  expect_equal(w$freq, c(0.6, 0.3))
  j <- gene_words(c(rep("TRAJ42", 6), rep(NA, 4)))
  expect_equal(j$freq, 0.6)                 # missing J: shorter stack
  expect_equal(nrow(gene_words(paste0("TRBV", 1:9), 0, max_words = 4)), 4)
})

test_that("word logos hang outside the edge, upright, most frequent nearest", {
  w <- tibble::tibble(word = c("AV12-2", "BV3"), freq = c(0.75, 0.25))
  # horizontal edge from (0, 0) to (6, 0): outside = below (right of the direction)
  p <- tcrhex:::.word_logo_polys(w, c(0, 0), c(6, 0), offset = 0.35, height = 2.6)
  expect_lte(max(p$y), -0.35 + 0.1)                       # digits overshoot cap height a bit
  expect_gte(min(p$y), -0.35 - 2.6 - 0.1)
  first <- p[p$glyph <= 6, ]                               # "AV12-2" (6 glyphs)
  second <- p[p$glyph > 6, ]
  expect_gt(mean(first$y), mean(second$y))                 # nearest the edge
  expect_equal(range(p$x), c(0.24, 5.76), tolerance = 0.25)
})

test_that("hexrows gene logos sit on the two lower edges", {
  cells <- tcrhex_example_cells()
  hd <- hexrows_data(cells, queries = top_public(cells, 1), gene_logos = TRUE)
  expect_true(all(c("a", "b", "va_genes", "ja_genes", "vb_genes", "jb_genes") %in%
                    hd$logos$chain))
  g <- hd$logos[hd$logos$chain %in% c("va_genes", "ja_genes") & hd$logos$bin == 0, ]
  H <- hd$meta$Hmax
  # below the lower-left edge: under the line from (-H cos30, -H/2) to (0, -H)
  edge_y <- -H / 2 - (g$x + H * cos(pi / 6)) * tan(pi / 6)
  expect_true(all(g$y <= edge_y + 0.1))
  # V on the first half of the edge (left), J on the second
  p0 <- c(-H * cos(pi / 6), -H / 2)
  u <- (g$x - p0[1]) * cos(pi / 6) + (g$y - p0[2]) * (-sin(pi / 6))   # along the edge
  L <- H
  expect_lte(max(u[g$chain == "va_genes"]), 0.5 * L + 1e-6)
  expect_gte(min(u[g$chain == "ja_genes"]), 0.5 * L - 1e-6)
  expect_true(all(g$x <= 0.5 + 1e-6))
  expect_true(hd$meta$gene_logos)
  p <- hexrows(cells, ggplot2::aes(fill = disease), queries = top_public(cells, 2),
               gene_logos = TRUE)
  expect_silent(ggplot2::ggplot_build(p))
})
