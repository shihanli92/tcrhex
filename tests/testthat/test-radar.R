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
  bw <- diff(range(rd$box$x)); bh <- diff(range(rd$box$y))
  expect_true(all(abs(ctr$x) - ctr$H * cos(pi / 6) > bw / 2 | abs(ctr$y) - ctr$H > bh / 2))
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
