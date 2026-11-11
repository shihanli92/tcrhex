# Build the tcrdist3 regression fixtures from the reference distances written by
# data-raw/tcrdist3_reference.py (run in a Python env with tcrdist3 0.3).
#   pub_*: clones from scRepertoire's public example data `contig_list`
#          (read_10x_contigs -> tcr_clones), queries = top_public(cells, 20).
#   syn_*: SYNTHETIC random CDR3s / mutants (not real TCRs) to stress gaps.
# Each fixture row holds both TCRs and tcrdist3's d_a / d_b.
args <- commandArgs(TRUE)  # dir with pub_q/pub_c/pub_ref/syn_q/syn_c/syn_ref .tsv
dir <- args[1]
out <- "tests/testthat/fixtures"
build <- function(prefix, keep_all = FALSE, n_random = 4000) {
  q <- read.delim(file.path(dir, paste0(prefix, "_q.tsv")), colClasses = "character")
  c <- read.delim(file.path(dir, paste0(prefix, "_c.tsv")), colClasses = "character")
  r <- read.delim(file.path(dir, paste0(prefix, "_ref.tsv")),
                  colClasses = c("character", "character", "integer", "integer"))
  set.seed(1)
  if (!keep_all) {
    near <- which(r$d_a + r$d_b <= 150)
    far <- sample(setdiff(seq_len(nrow(r)), near), n_random)
    r <- r[sort(c(near, far)), ]
  }
  qi <- match(r$qid, q$id)
  ci <- match(r$cid, c$id)
  d <- data.frame(q_va = q$va[qi], q_ca = q$ca[qi], q_vb = q$vb[qi], q_cb = q$cb[qi],
                  c_va = c$va[ci], c_ca = c$ca[ci], c_vb = c$vb[ci], c_cb = c$cb[ci],
                  d_a = r$d_a, d_b = r$d_b)
  f <- gzfile(file.path(out, paste0("tcrdist3_", prefix, ".tsv.gz")), "w")
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  close(f)
  cat(prefix, nrow(d), "pairs\n")
}
build("pub")
build("syn", keep_all = FALSE, n_random = 6000)
