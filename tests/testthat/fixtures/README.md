# Test fixtures

* `contigs_10x_synthetic.csv`, `airr_synthetic.tsv`, `screp_synthetic.tsv`,
  `meta_synthetic.csv`: SYNTHETIC cells (invented TCRs / donors, not real
  data), the same 6 cells in the three reader formats.
* `tcrdist3_pub.tsv.gz`: query x clone pairs with tcrdist3 0.3's d_a / d_b,
  clones from scRepertoire's public example data `contig_list` (MIT);
  pairs with d <= 150 plus 4000 random pairs.
* `tcrdist3_syn.tsv.gz`: SYNTHETIC random CDR3s (length 1-30) and indel mutants,
  with tcrdist3's d_a / d_b, to exercise the gap placement.
  Built by `data-raw/tcrdist3_reference.py` + `data-raw/make_tcrdist3_fixtures.R`.
