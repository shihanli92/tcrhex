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
