cells <- tcrhex_example_cells()   # synthetic

test_that("hexrows_data: per-bin table invariants", {
  q <- top_public(cells, 3)
  hd <- hexrows_data(cells, ggplot2::aes(fill = disease), queries = q, max_hex = 5)
  b <- hd$bins
  expect_named(b, c("row", "bin", "level", "clones", "shown"))
  per <- aggregate(cbind(clones, shown) ~ row + bin, data = b, FUN = sum)
  expect_true(all(per$shown == pmin(per$clones, 5)))
  # composition kept: per bin, shown = largest-remainder allotment
  for (i in seq_len(nrow(per))) {
    s <- b[b$row == per$row[i] & b$bin == per$bin[i], ]
    expect_equal(unname(s$shown), unname(allot(setNames(s$clones, s$level), 5)))
  }
  # within a level, the closest clones are shown
  h <- hd$hits
  for (k in split(h, paste(h$row, h$bin, h$level))) {
    if (any(!k$shown) && any(k$shown)) expect_lte(max(k$d[k$shown]), min(k$d[!k$shown]))
  }
  expect_equal(sum(b$clones), nrow(hd$hits))
  expect_true(all(hd$hits$d[hd$hits$bin == 0] == 0))
  # hexagons: 6 vertices each, inside their frame
  expect_equal(nrow(hd$hexes), 6 * sum(b$shown))
  expect_equal(levels(hd$hits$level), c("disease X", "disease Y", "healthy"))
})

test_that("fill evaluated per cell then summarised; factor levels give band order", {
  q <- top_public(cells, 2)
  hd <- hexrows_data(cells, ggplot2::aes(fill = tissue), queries = q, logos = FALSE)
  expect_true("multiple" %in% levels(hd$hits$level))
  hd2 <- hexrows_data(cells, ggplot2::aes(fill = tissue_class(tissue)), queries = q,
                      summarise_fill = summarise_tissue, logos = FALSE)
  expect_true(all(levels(hd2$hits$level) %in% levels(tissue_class(character()))))
  expect_true("blood + tissue" %in% hd2$hits$level)
  # expression and per-row variables
  qq <- q
  qq$allele <- c("A*02", "B*08")
  hd3 <- hexrows_data(cells, ggplot2::aes(fill = hla_carrier(hla, .row_allele)),
                      queries = qq, logos = FALSE)
  expect_setequal(levels(hd3$hits$level), c("carrier", "non-carrier", "not typed"))
  # within a bin, shown hexagons are banded by level order (top rows first)
  hx <- hd$hexes[!duplicated(hd$hexes$group), ]
  for (k in split(hx, paste(hx$row, hx$bin))) {
    expect_false(is.unsorted(as.integer(k$level)))
  }
})

test_that("default folding to max 8 levels for non-factor fills", {
  x <- tcrhex:::.make_levels(c(letters[1:10], rep("a", 5), NA), 8)
  expect_equal(levels(x$value), c("a", letters[2:7], "Other", "unknown"))
  y <- tcrhex:::.make_levels(factor(letters[1:10], levels = rev(letters[1:10])), 8)
  expect_equal(levels(y$value), rev(letters[1:10]))
})

test_that("hexrows returns a ggplot that takes scales and builds", {
  q <- top_public(cells, 2)
  p <- hexrows(cells, ggplot2::aes(fill = disease), queries = q)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_true(length(b$data) >= 5)
  p2 <- p + ggplot2::scale_fill_manual(values = c("disease X" = "red", "disease Y" = "blue",
                                                  healthy = "grey"))
  b2 <- ggplot2::ggplot_build(p2)
  fills <- unique(unlist(lapply(b2$data, function(d) d$fill)))
  expect_true(all(c("red", "blue", "grey") %in% fills))
  p3 <- hexrows(cells, ggplot2::aes(fill = tissue_class(tissue), pattern = tissue_class(tissue)),
                queries = "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF",
                summarise_fill = summarise_tissue) + scale_fill_tissue() + scale_pattern_tissue()
  expect_silent(ggplot2::ggplot_build(p3))
  f <- tempfile(fileext = ".png")
  ggsave_hexrows(f, p3)
  expect_true(file.exists(f))
})

test_that("row tags: top_public concordance, muted rows, user row_tag", {
  q <- top_public(cells, 3)
  q$concordance[2] <- 0.05
  hd <- hexrows_data(cells, ggplot2::aes(fill = disease), queries = q, logos = FALSE)
  expect_match(hd$rows$tag[1], "of paired people")
  expect_match(hd$rows$tag[2], "no dominant")
  expect_true(hd$rows$muted[2])
  hd2 <- hexrows_data(cells, queries = q, row_tag = c("x", "y", "z"), logos = FALSE)
  expect_equal(hd2$rows$tag, c("x", "y", "z"))
  hd3 <- hexrows_data(cells, queries = c(myTCR = "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF"),
                      logos = FALSE)
  expect_equal(hd3$rows$tag, "myTCR")
})

test_that("default fill scale: fixed hues in level order, greys for Other / unknown", {
  cl <- tcrhex_example_cells(n_donors = 12)
  cl$grp <- c(letters[1:11], NA)[match(cl$donor, sort(unique(cl$donor)))]
  p <- hexrows(cl, ggplot2::aes(fill = grp), queries = top_public(cl, 2), logos = FALSE)
  hd <- attr(p, "hexrows")
  lv <- levels(hd$hits$level)
  expect_equal(utils::tail(lv, 2), c("Other", "unknown"))
  expect_length(lv, 9)
  b <- ggplot2::ggplot_build(p)
  hexdat <- b$data[[3]]
  map <- unique(data.frame(level = as.character(hd$hexes$level[!duplicated(hd$hexes$group)]),
                           fill = hexdat$fill[!duplicated(hexdat$group)]))
  expect_equal(map$fill[map$level == "Other"][1], "#b4b3ac")
  expect_equal(map$fill[map$level == "unknown"][1], "#d6d5ce")
  expect_equal(map$fill[map$level == lv[1]][1], "#2a78d6")
  lab <- b$plot$scales$get_scales("fill")$get_labels()
  expect_match(lab[1], "\\(\\d+\\)$")
})
