cells <- tcrhex_example_cells(n_donors = 20)   # synthetic
q1 <- "TRAV27:CAGAGSQGNLIF:TRBV6-5:CASSYSGGSYEQYF"

test_that("hexradar: one query, frames on their rings, nothing overlaps", {
  rd <- hexradar_data(cells, ggplot2::aes(fill = disease), query = q1, radius = 48)
  expect_s3_class(rd, "hexradar_data")
  f <- rd$frames[!duplicated(rd$frames$group), ]
  ctr <- dplyr::summarise(dplyr::group_by(rd$frames, .data$bin), x = mean(.data$x),
                          y = mean(.data$y), H = max(sqrt((.data$x - mean(.data$x))^2 +
                                                          (.data$y - mean(.data$y))^2)))
  expect_equal(nrow(ctr), 5)                                  # 0, 1-12 ... 37-48
  ring_r <- sort(unique(round(sqrt(rd$rings$x^2 + rd$rings$y^2), 6)))
  r <- sqrt(ctr$x^2 + ctr$y^2)
  expect_equal(r, ring_r, tolerance = 1e-6)                   # frame k sits on ring k
  d <- as.matrix(dist(cbind(ctr$x, ctr$y)))
  for (i in 1:4) for (j in (i + 1):5) expect_gt(d[i, j], ctr$H[i] + ctr$H[j])
  # circle centre (default): frames clear the circle
  cr <- max(sqrt(rd$box$x^2 + rd$box$y^2))
  expect_equal(diff(range(sqrt(rd$box$x^2 + rd$box$y^2))), 0, tolerance = 1e-9)
  expect_true(all(r - ctr$H > cr))
  expect_equal(rd$box_text$label, c("query", "TRAV27", "CAGAGSQGNLIF", "TRBV6-5", "CASSYSGGSYEQYF"))
  # the circle holds every text line (mono width estimate)
  cw <- 0.6 * 7.5 / 72 / 0.072
  lh <- 1.3 * 7.5 / 72 / 0.072
  corner <- sqrt((nchar(rd$box_text$label) * cw / 2)^2 + (abs(rd$box_text$y) + lh / 2)^2)
  expect_true(all(corner < cr))
  # box option keeps the rounded rectangle, one line per chain
  rb <- hexradar_data(cells, ggplot2::aes(fill = disease), query = q1, centre = "box")
  expect_length(rb$box_text$label, 3)
  bw <- diff(range(rb$box$x)); bh <- diff(range(rb$box$y))
  cb <- dplyr::summarise(dplyr::group_by(rb$frames, .data$bin), x = mean(.data$x), y = mean(.data$y),
                         H = max(sqrt((.data$x - mean(.data$x))^2 + (.data$y - mean(.data$y))^2)))
  expect_true(all(abs(cb$x) - cb$H * cos(pi / 6) > bw / 2 | abs(cb$y) - cb$H > bh / 2))
  # bins / shown follow the row rules
  per <- aggregate(cbind(clones, shown) ~ bin, data = rd$bins, FUN = sum)
  expect_true(all(per$shown == pmin(per$clones, 50)))
  expect_match(rd$labels$label[1], "TCRdist 0 \\(identical\\)")
  expect_match(rd$labels$label[2], "clones · \\d+ donors")
})

test_that("hexradar: errors for several queries; options build", {
  expect_error(hexradar(cells, query = top_public(cells, 2)), "one query")
  p <- hexradar(cells, ggplot2::aes(fill = tissue_class(tissue), pattern = tissue_class(tissue)),
                query = top_expanded(cells, 1), summarise_fill = summarise_tissue,
                cumulative = TRUE, logos = TRUE, gene_logos = TRUE, max_hex = 19,
                donor_noun = "people") + scale_fill_tissue() + scale_pattern_tissue()
  expect_silent(ggplot2::ggplot_build(p))
  rd <- attr(p, "hexrows")
  expect_true(all(c("a", "b", "va_genes", "vb_genes") %in% rd$logos$chain))
  expect_match(rd$labels$label[2], "TCRdist ≤ 12")
  expect_match(rd$labels$label[2], "people")
  f <- tempfile(fileext = ".png")
  ggsave_hexrows(f, p)
  expect_true(file.exists(f))
})

test_that("hexradar: many 12-wide bins stay compact; empty bins get short labels", {
  rd <- hexradar_data(cells, ggplot2::aes(fill = disease), query = q1, radius = 240)
  expect_equal(rd$meta$n_bins, 21)
  ctr <- dplyr::summarise(dplyr::group_by(rd$frames, .data$bin), x = mean(.data$x),
                          y = mean(.data$y), H = max(sqrt((.data$x - mean(.data$x))^2 +
                                                          (.data$y - mean(.data$y))^2)))
  d <- as.matrix(dist(cbind(ctr$x, ctr$y)))
  hh <- outer(ctr$H, ctr$H, "+")
  expect_true(all(d[upper.tri(d)] > hh[upper.tri(hh)]))       # no two frames overlap
  expect_lt(diff(rd$meta$xlim) * 0.072, 20)                     # page stays reasonable
  lab_empty <- rd$labels$label[rd$labels$empty]
  if (length(lab_empty)) expect_true(all(grepl("· none$", lab_empty)))
  # the ring-label ray avoids every spoke
  a <- atan2(rd$ring_labels$y[1], rd$ring_labels$x[1])
  th <- atan2(ctr$y, ctr$x)
  dd <- abs(((th - a + pi) %% (2 * pi)) - pi)
  expect_gt(min(dd), 0.05)
})
