# Render the HPC bundle with tcrhex, next to the Python reference PNGs.
# The bundle is research data: keep it in data-raw/bundle/ (git-ignored, not
# in the built package). Usage, from the package folder:
#   Rscript data-raw/render_bundle.R [bundle_dir]
# Writes <bundle>/tcrhex_top20_rows_CD8mem_r108_{tissue,disease,hla_auto}.png
# and the matching _table.tsv files, and reports differences from the Python
# per-bin tables.
suppressMessages({library(tcrhex); library(ggplot2)})
args <- commandArgs(TRUE)
d <- if (length(args)) args[1] else "data-raw/bundle"

q <- read.delim(file.path(d, "queries", "top20_queries_CD8mem.tsv"), colClasses = "character")
q$rank <- as.integer(q$rank)
q$concordance <- as.numeric(q$concordance)
files <- sort(list.files(file.path(d, "hits"), "^rank[0-9]+_hits\\.tsv$", full.names = TRUE))
hits <- do.call(rbind, lapply(seq_along(files), function(i) {
  h <- read.delim(files[i], comment.char = "#", colClasses = "character")
  h$query <- q$name[i]
  h
}))
for (v in c("d", "d_a", "d_b", "cells")) hits[[v]] <- as.integer(hits[[v]])
hits$donor <- hits$person_id
hits$.clone_id <- seq_len(nrow(hits))

out <- function(tag, p) {
  ggsave_hexrows(file.path(d, sprintf("tcrhex_top20_rows_CD8mem_r108_%s.png", tag)), p)
  hd <- attr(p, "hexrows")
  write.table(hd$bins, file.path(d, sprintf("tcrhex_top20_rows_CD8mem_r108_%s_table.tsv", tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ref <- file.path(d, "reference", sprintf("top20_rows_CD8mem_r108_%s_table.tsv", tag))
  if (file.exists(ref)) {
    r <- read.delim(ref)
    a <- as.data.frame(hd$bins)
    a$level <- as.character(a$level)
    m <- merge(a, r[, c("row", "bin", "level", "clones", "shown")],
               by = c("row", "bin", "level"), all = TRUE, suffixes = c(".R", ".py"))
    bad <- m[is.na(m$clones.R) | is.na(m$clones.py) | m$clones.R != m$clones.py |
               m$shown.R != m$shown.py, ]
    cat(sprintf("%s: %d per-bin rows, %d differ from Python\n", tag, nrow(m), nrow(bad)))
    if (nrow(bad)) print(utils::head(bad, 20))
  }
}

out("tissue", hexrows(NULL, aes(fill = tissue_class(tissue_class), pattern = tissue_class(tissue_class)),
                      queries = q, hits = hits, summarise_fill = summarise_tissue) +
  scale_fill_tissue() + scale_pattern_tissue())
out("disease", hexrows(NULL, aes(fill = disease), queries = q, hits = hits))

# HLA auto: per-row enriched allele among donors within 24, colour by carriers
clones <- hits[!duplicated(hits$donor), c("donor", "hla_1field")]
enr <- hla_enriched(hits, clones, radius = 24, hla = "hla_1field", queries = q)
q$allele <- enr$allele
out("hla_auto", hexrows(NULL, aes(fill = hla_carrier(hla_1field, .row_allele)),
                        queries = q, hits = hits, row_tag = enr$text))
