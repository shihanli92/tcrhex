q1 <- "TRAV27:CAGAGSQGNLIF:TRBV19:CASSIRSSYEQYF"
clone <- function(va, ca, vb, cb) data.frame(.clone_id = 1L, donor = "D", va = va, ca = ca,
                                             vb = vb, cb = cb, cells = 1L)

test_that("hand-checked paired TCRdist cases", {
  # identical
  h <- tcrdist_paired(q1, clone("TRAV27", "CAGAGSQGNLIF", "TRBV19", "CASSIRSSYEQYF"), Inf)
  expect_equal(c(h$d, h$d_a, h$d_b), c(0, 0, 0))
  # one CDR3b substitution inside the trimmed region: I -> L (BLOSUM62 2 -> 4-2 = 2), x3
  h <- tcrdist_paired(q1, clone("TRAV27", "CAGAGSQGNLIF", "TRBV19", "CASSLRSSYEQYF"), Inf)
  expect_equal(h$d_b, 6)
  # substitution in the trimmed N-terminal 3 residues counts 0
  h <- tcrdist_paired(q1, clone("TRAV27", "CAGAGSQGNLIF", "TRBV19", "CATSIRSSYEQYF"), Inf)
  expect_equal(h$d_b, 0)
  # one residue shorter: gap 4 x weight 3 = 12, best gap position -> no mismatches
  h <- tcrdist_paired(q1, clone("TRAV27", "CAGAGSQGNLIF", "TRBV19", "CASSIRSYEQYF"), Inf)
  expect_equal(h$d_b, 12)
  # V-only difference: loops only (no CDR3 term)
  h <- tcrdist_paired(q1, clone("TRAV27", "CAGAGSQGNLIF", "TRBV6-5", "CASSIRSSYEQYF"), Inf)
  v <- tcrhex:::.vloop_dist("TRBV19", "TRB")[["TRBV6-5"]]
  expect_equal(h$d_b, v)
  expect_gt(v, 0)
  vt <- vgene_loops()
  a <- vt[vt$gene == "TRBV19", ]
  b <- vt[vt$gene == "TRBV6-5", ]
  dm <- tcrhex:::.tcr_dm()
  aa <- strsplit("ARNDCQEGHILKMFPSTWYV", "")[[1]]
  idx <- function(s) { i <- match(strsplit(s, "")[[1]], aa); i[is.na(i)] <- 21; i }
  manual <- sum(vapply(c("cdr1", "cdr2", "cdr25"), function(l)
    sum(dm[cbind(idx(a[[l]]), idx(b[[l]]))]), 1))
  expect_equal(v, manual)
})

test_that("substitution matrix is min(4, 4 - BLOSUM62)", {
  dm <- tcrhex:::.tcr_dm()
  expect_equal(dim(dm), c(21L, 21L))
  expect_true(all(diag(dm)[1:20] == 0))
  expect_true(all(dm[21, ] == 0))
  expect_equal(dm[1, 16], 3)   # A-S: BLOSUM 1
  expect_equal(dm[10, 20], 1)  # I-V: BLOSUM 3
  expect_true(isSymmetric(dm))
})

check_tcrdist3 <- function(file) {
  r <- read.delim(gzfile(file), colClasses = c(rep("character", 8), "integer", "integer"))
  r$q <- paste(r$q_va, r$q_ca, r$q_vb, r$q_cb, sep = ":")
  cl <- unique(r[, c("c_va", "c_ca", "c_vb", "c_cb")])
  names(cl) <- c("va", "ca", "vb", "cb")
  cl$.clone_id <- seq_len(nrow(cl))
  r$.clone_id <- match(paste(r$c_va, r$c_ca, r$c_vb, r$c_cb),
                       paste(cl$va, cl$ca, cl$vb, cl$cb))
  h <- tcrdist_paired(unique(r$q), cl, radius = Inf)
  m <- merge(r, h, by.x = c("q", ".clone_id"), by.y = c("query", ".clone_id"),
             suffixes = c(".tcrdist3", ""))
  expect_equal(nrow(m), nrow(r))
  bad <- m[m$d_a != m$d_a.tcrdist3 | m$d_b != m$d_b.tcrdist3, ]
  if (nrow(bad)) {
    print(utils::head(bad[, c("q", "c_va", "c_ca", "c_vb", "c_cb", "d_a.tcrdist3", "d_a",
                              "d_b.tcrdist3", "d_b")], 20))
  }
  expect_equal(nrow(bad), 0)
}

test_that("d_a / d_b match tcrdist3 exactly (public scRepertoire clones)", {
  check_tcrdist3(fx("tcrdist3_pub.tsv.gz"))
})

test_that("d_a / d_b match tcrdist3 exactly (synthetic gap stress set)", {
  check_tcrdist3(fx("tcrdist3_syn.tsv.gz"))
})

test_that("radius filter and unknown V genes", {
  cl <- rbind(clone("TRAV27", "CAGAGSQGNLIF", "TRBV19", "CASSIRSSYEQYF"),
              clone("TRAV27", "CAGAGSQGNLIF", "TRBVX", "CASSIRSSYEQYF"),
              clone("TRAV1-1", "CWWWWWWWWWWWF", "TRBV2", "CWWWWWWWWWWWWWWWF"))
  cl$.clone_id <- 1:3
  expect_warning(h <- tcrdist_paired(q1, cl, 108), "1 of 3 clones dropped")
  expect_equal(h$.clone_id, 1L)
  expect_equal(attr(h, "dropped_unknown_v"), 1)
})

test_that("queries are validated before computing", {
  expect_error(as_queries("TRAV27:CAGAGSQGNLIF:TRBV19"), "must look like")
  expect_error(as_queries("TRAVX:CAGAGSQGNLIF:TRBV19:CASSIRSSYEQYF"), "unknown TRAV")
  expect_error(as_queries("TRAV27:AGAGSQGNLIF:TRBV19:CASS1RSSYEQYF"), "CDR3a .* not")
  q <- as_queries(c(flu = "TRAV27*01:cagagsqgnlif:TRBV19:CASSIRSSYEQYF"))
  expect_equal(q$name, "flu")
  expect_equal(q$va, "TRAV27")
  expect_equal(q$ca, "CAGAGSQGNLIF")
  df <- data.frame(name = "x", va = "TRAV27", ca = "CAGAGSQGNLIF", vb = "TRBV19",
                   cb = "CASSIRSSYEQYF", epitope = "flu M1")
  expect_equal(as_queries(df)$epitope, "flu M1")
})
