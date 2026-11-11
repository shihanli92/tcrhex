# README figures (public scRepertoire example data + SYNTHETIC cells)
suppressMessages({library(tcrhex); library(scRepertoire); library(ggplot2)})
data("contig_list"); data("scRep_example")
samples <- c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
combined <- combineTCR(contig_list, samples = samples)
cells <- read_screpertoire(combined, donor = sub("[BL]$", "", sample))
p <- hexrows(cells, aes(fill = sample), queries = top_public(cells, 10))
ggsave_hexrows("man/figures/ex1_screp_sample.png", p, dpi = 110)
meta <- data.frame(sample = samples, patient = sub("[BL]$", "", samples),
                   tissue = ifelse(grepl("B$", samples), "peripheral blood", "lung"))
cells <- read_10x_contigs(setNames(contig_list, samples), meta = meta, by = "sample",
                          donor = patient)
p <- hexrows(cells, aes(fill = tissue), queries = top_public(cells, 10)) +
  scale_fill_manual(values = c("peripheral blood" = "#0072B2", lung = "#E69F00",
                               multiple = "#56B4E9"))
ggsave_hexrows("man/figures/ex2_10x_tissue.png", p, dpi = 110)
cells <- tcrhex_example_cells(n_donors = 40)
p <- hexrows(cells, aes(fill = hla_carrier(hla, "A*02")),
             queries = c("my TCR" = "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF"),
             title = "SYNTHETIC example data: HLA-A*02 carriers") +
  scale_fill_manual(values = c(carrier = "#D55E00", "non-carrier" = "#0072B2",
                               "not typed" = "#d6d5ce"))
ggsave_hexrows("man/figures/ex4_synthetic_hla.png", p, dpi = 110)
p <- hexrows(cells, aes(fill = tissue_class(tissue), pattern = tissue_class(tissue)),
             queries = top_public(cells, 3), summarise_fill = summarise_tissue,
             title = "SYNTHETIC example data: tissue classes") +
  scale_fill_tissue() + scale_pattern_tissue()
ggsave_hexrows("man/figures/ex3_synthetic_tissue.png", p, dpi = 110)
q1 <- top_public(cells, 1)
p <- hexradar(cells, aes(fill = hla_carrier(hla, "A*02")), query = q1, donor_noun = "people",
              title = "SYNTHETIC example data: hexradar()") +
  scale_fill_manual(values = c(carrier = "#2a78d6", "non-carrier" = "#b4b3ac",
                               "not typed" = "#e4e3dc"))
ggsave_hexrows("man/figures/ex5_synthetic_radar.png", p, dpi = 110)
