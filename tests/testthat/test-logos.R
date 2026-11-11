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
  expect_setequal(unique(g$letter), c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "-"))
  expect_true(all(g$x >= 0 & g$x <= 1 & g$y >= 0 & g$y <= 1))
  p <- tcrhex:::.logo_polys(logo_freqs(logo_align(c("CASSF", "CASSF"), 1)), c(0, 0), c(3, 0))
  # edge along x: letters stand on the baseline 0.35 above, full height 2.6
  expect_equal(min(p$y), 0.35, tolerance = 1e-6)
  expect_equal(max(p$y), 0.35 + 2.6, tolerance = 1e-6)
  expect_equal(range(p$x), c(0.04, 3 - 0.04), tolerance = 1e-6)
})
