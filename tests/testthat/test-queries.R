# SYNTHETIC fixture (invented TCRs and donors)
tp_cells <- function() {
  A1 <- c("TRAV12-2", "CAVNGGSQGNLIF"); A2 <- c("TRAV27", "CAGAGSQGNLIF")
  A3 <- c("TRAV21", "CAVRDSNYQLIW"); A4 <- c("TRAV5", "CAESTGKLIF")
  B1 <- c("TRBV19", "CASSIRSSYEQYF"); B2 <- c("TRBV20-1", "CSARDRTGNGYTF")
  B3 <- c("TRBV9", "CASSVGGSSYNEQFF")
  row <- function(d, a, b) data.frame(donor = d, va = a[1], ca = a[2], vb = b[1], cb = b[2])
  na <- c(NA, NA)
  df <- rbind(row("D1", A1, B1), row("D1", A1, B1), row("D1", A2, B1),
              row("D2", A1, B1), row("D2", A2, B1),
              row("D3", na, B1),
              row("D1", A3, B2),
              row("D2", A4, B2), row("D2", A4, B2), row("D2", A4, B2),
              row("D4", na, B3))
  df$cell <- paste0("c", seq_len(nrow(df)))
  cells_from_df(df, verbose = FALSE)
}

test_that("top_public: donors, dominant alpha, concordance", {
  expect_message(q <- top_public(tp_cells(), n = 3), "only 2 beta")
  expect_equal(q$vb, c("TRBV19", "TRBV20-1"))
  expect_equal(q$va, c("TRAV12-2", "TRAV5"))
  expect_equal(q$donors, c(3L, 2L))
  expect_equal(q$paired_donors, c(2L, 2L))
  expect_equal(q$donors_with_alpha, c(2L, 1L))
  expect_equal(q$concordance, c(1, 0.5))
  expect_equal(q$rank, 1:2)
  expect_equal(q$name[1], "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF")
})

test_that("tcr_clones summarises every metadata column by one rule", {
  df <- data.frame(cell = paste0("c", 1:4), donor = c("D1", "D1", "D1", "D2"),
                   va = "TRAV27", ca = "CAGAGSQGNLIF", vb = "TRBV19", cb = "CASSIRSSYEQYF",
                   tissue = c("blood", "lung", "blood", "lung"), age = c(40, 40, 40, 55),
                   batch = c("b1", "b1", NA, "b2"))
  cl <- tcr_clones(cells_from_df(df, verbose = FALSE))
  expect_equal(cl$cells, c(3L, 1L))
  expect_equal(cl$tissue, c("multiple", "lung"))
  expect_equal(cl$age, c(40, 55))
  expect_equal(cl$batch, c("b1", "b2"))
  cl2 <- tcr_clones(cells_from_df(df, verbose = FALSE), summarise = summarise_most_common)
  expect_equal(cl2$tissue, c("blood", "lung"))
})

test_that("top_expanded: most cells per group, ties, labels", {
  # SYNTHETIC cells
  row <- function(d, tissue, va, ca, k) data.frame(donor = d, tissue = tissue, va = va, ca = ca,
                                                  vb = "TRBV19", cb = "CASSIRSSYEQYF")[rep(1, k), ]
  df <- rbind(row("D1", "lung", "TRAV27", "CAGAGSQGNLIF", 3),
              row("D1", "blood", "TRAV27", "CAGAGSQGNLIF", 5),
              row("D2", "lung", "TRAV12-2", "CAVNGGSQGNLIF", 4),
              row("D2", "blood", "TRAV21", "CAVRDSNYQLIW", 1),
              row("D3", NA, "TRAV21", "CAVRDSNYQLIW", 9))
  df$cell <- paste0("c", seq_len(nrow(df)))
  x <- cells_from_df(df, verbose = FALSE)
  q <- top_expanded(x, n = 1, by = tissue)
  expect_equal(q$group, c("blood", "lung"))
  expect_equal(q$va, c("TRAV27", "TRAV12-2"))
  expect_equal(q$cells_in_group, c(5L, 4L))
  expect_equal(q$cells, c(8L, 4L))
  expect_match(q$label[2], "top 1 in lung: 4 cells, donor D2")
  # the same TCR can be top in several groups; names stay unique
  q2 <- top_expanded(x, n = 2, by = tissue)
  expect_equal(nrow(q2), 4)
  expect_false(anyDuplicated(q2$name) > 0)
  expect_equal(q2$va[q2$group == "lung"], c("TRAV12-2", "TRAV27"))
  # overall (no by): D3's 9-cell clone wins, NA tissue is fine without `by`
  expect_equal(top_expanded(x, n = 1)$cells, 9L)
  # factor levels give the group order; works as hexrows queries with tags
  x$tissue <- factor(x$tissue, levels = c("lung", "blood"))
  q3 <- top_expanded(x, n = 1, by = tissue)
  expect_equal(q3$group, c("lung", "blood"))
  hd <- hexrows_data(x, ggplot2::aes(fill = tissue), queries = q3, logos = FALSE)
  expect_equal(hd$rows$tag, q3$label)
})
