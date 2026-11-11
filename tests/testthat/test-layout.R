test_that("hex_layout: k = 1, 7, 19 positions and banding order", {
  l1 <- hex_layout(1)
  expect_equal(c(l1$x, l1$y), c(0, 0))
  l7 <- hex_layout(7)
  expect_equal(nrow(l7), 7)
  r <- sqrt(l7$x^2 + l7$y^2)
  expect_equal(sum(abs(r) < 1e-9), 1)
  expect_true(all(abs(r[abs(r) > 1e-9] - sqrt(3)) < 1e-9))
  # rows top to bottom, left to right within a row
  expect_equal(round(l7$y, 6), c(1.5, 1.5, 0, 0, 0, -1.5, -1.5))
  expect_true(all(diff(l7$x[1:2]) > 0), all(diff(l7$x[3:5]) > 0))
  l19 <- hex_layout(19)
  expect_equal(as.numeric(table(round(l19$y, 6))), c(3, 4, 5, 4, 3))
  expect_equal(round(max(sqrt(l19$x^2 + l19$y^2)), 6), round(2 * sqrt(3), 6))
  # every pair of centres at least the lattice spacing apart
  d <- as.matrix(dist(cbind(l19$x, l19$y)))
  expect_gte(min(d[upper.tri(d)]), sqrt(3) - 1e-9)
  # y non-increasing in slot order (banding)
  expect_true(all(diff(round(l19$y, 6)) <= 0))
})

test_that("hex_layout ties broken by angle", {
  l2 <- hex_layout(2)
  # second slot: of the 6 neighbours, smallest atan2. y = -1.5 * r gives -0 for
  # r = 0 and atan2(-0, x < 0) = -pi (as in Python), so the left neighbour wins.
  expect_equal(nrow(l2), 2)
  nb <- l2[abs(l2$x) + abs(l2$y) > 0, ]
  expect_equal(c(round(nb$x, 6), nb$y), c(round(-sqrt(3), 6), 0))
  l3 <- hex_layout(3)  # next: -120 degrees
  expect_true(any(abs(l3$x + sqrt(3) / 2) < 1e-9 & l3$y == -1.5))
})

test_that("allot: largest remainder sums to k and keeps proportions", {
  expect_equal(allot(c(a = 30, b = 10, c = 2), 19), c(a = 14L, b = 4L, c = 1L))
  for (i in 1:200) {
    n <- sample(1:50, sample(1:6, 1), replace = TRUE)
    k <- sample(1:25, 1)
    a <- allot(n, k)
    expect_equal(sum(a), min(k, sum(n)))
    expect_true(all(a <= n))
    if (sum(n) > k) expect_true(all(abs(a - k * n / sum(n)) < 1))
  }
  expect_equal(allot(c(3, 2), 19), c(3L, 2L))
})

test_that("binning: 0 alone, 1-12, last bin capped at the radius", {
  expect_equal(tcrdist_bin(c(0, 1, 12, 13, 24, 25, 107, 108, 109), 12, 108),
               c(0L, 1L, 1L, 2L, 2L, 3L, 9L, 9L, NA))
  b <- bin_labels(12, 108)
  expect_equal(nrow(b), 10)
  expect_equal(b$label[1:3], c("0", "1\u201312", "13\u201324"))
  expect_equal(b$label[10], "97\u2013108")
  b2 <- bin_labels(12, 100)
  expect_equal(b2$label[nrow(b2)], "97\u2013100")
})
