test_that("tissue_class keyword rules", {
  x <- c("PBMC", "Peripheral blood", "lung tumor", "Colon mucosa", "CSF",
         "cerebrospinal fluid", "Tonsil", "BAL lavage", "skin blister", "", NA,
         "cell culture", "liver", "frontal lobe", "blood + tissue")
  expect_equal(as.character(tissue_class(x)),
               c("blood", "blood", "tumour (labelled)", "gut", "other tissue",
                 "CNS / CSF", "lymphoid / bone marrow", "lung / airway", "skin",
                 "unknown", "unknown", "unknown", "other tissue", "CNS / CSF",
                 "blood + tissue"))
  expect_s3_class(tissue_class("blood"), "factor")
})

test_that("summarise_tissue", {
  expect_equal(as.character(summarise_tissue(c("PBMC", "blood"))), "blood")
  expect_equal(as.character(summarise_tissue(c("PBMC", "lung"))), "blood + tissue")
  expect_equal(as.character(summarise_tissue(c("lung", "colon"))), "other / multiple tissues")
  expect_equal(as.character(summarise_tissue(c("lung", NA))), "lung / airway")
  expect_equal(as.character(summarise_tissue(c(NA, ""))), "unknown")
})

test_that("hla_carrier at 1-field resolution, vectorised over alleles", {
  h <- c("A02:01;A24:02;B*08:01", "A*01:01,B*07:02", NA, "", "HLA-A*02:06|C*07:01")
  expect_equal(as.character(hla_carrier(h, "A*02")),
               c("carrier", "non-carrier", "not typed", "not typed", "carrier"))
  expect_equal(as.character(hla_carrier(h, "A02:01")), as.character(hla_carrier(h, "A*02")))
  expect_equal(as.character(hla_carrier(h[1:2], c("B*08", "B*07"))), c("carrier", "carrier"))
  expect_true(is.na(hla_carrier(h[1], NA)))
})

test_that("hla_enriched picks the over-represented allele (synthetic fixture)", {
  # 10 typed donors; A*24 in 2 of 10 overall, but in all 4 neighbourhood donors
  hla <- c(rep("A*24:02;A*01:01;B*07:02", 4), rep("A*02:01;A*01:01;B*07:02", 6))
  clones <- data.frame(.clone_id = 1:10, donor = paste0("D", 1:10), hla = hla)
  clones$hla[5:10] <- c("A*02:01;A*03:01;B*08:01", "A*02:01;A*03:01;B*08:01",
                        "A*11:01;A*01:01;B*07:02", "A*02:01;A*03:01;B*44:02",
                        "A*02:01;A*11:01;B*07:02", "A*02:01;A*03:01;B*35:01")
  hits <- data.frame(query = "q1", .clone_id = 1:5, d = c(0, 12, 20, 24, 90),
                     donor = paste0("D", 1:5))
  e <- hla_enriched(hits, clones, radius = 24)
  expect_equal(e$allele, "A*24")
  expect_equal(c(e$k, e$n), c(4L, 4L))
  expect_equal(e$p0, 0.4)
  expect_equal(e$p, 0.4^4)
  # min_carriers
  e2 <- hla_enriched(hits, clones, radius = 12)
  expect_true(is.na(e2$allele))
  expect_match(e$text, "A\\*24 in 4/4")
})
