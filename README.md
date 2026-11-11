# tcrhex

Paired TCRdist **hex rows** from single-cell TCR data, drawn with ggplot2.

One row per query TCR (TRAV + CDR3α, TRBV + CDR3β). Along the row, one
hexagonal frame per paired-TCRdist bin (0, 1–12, 13–24, … up to the search
radius). Inside each frame, one small hexagon per clone (donor × TRAV + CDR3α
× TRBV + CDR3β), up to 19, coloured by **any** categorical column of your
data. The frames also carry frequency logos of the bin's CDR3α and CDR3β.

* **Input**: ordinary single-cell TCR metadata: scRepertoire (`combineTCR()`
  lists, or Seurat / SingleCellExperiment objects after `combineExpression()`),
  10x Cell Ranger `filtered_contig_annotations.csv`, AIRR rearrangement TSV, or
  any data frame.
* **Paired TCRdist**: reproduces tcrdist3's human defaults exactly (tested
  against tcrdist3 0.3 on 308,380 query × clone pairs). Only queries are
  compared against all clones, in C++: 20 queries × 500,000 clones take
  about 5 s.
* **ggplot-like**: `hexrows(cells, aes(fill = <any column or expression>))`
  returns a ggplot, so `+ scale_fill_*()`, `+ labs()` and `+ theme()` work.

## Install

```r
# from the package folder (or a tarball)
install.packages("remotes")
remotes::install_local("tcrhex")          # or: R CMD INSTALL tcrhex_0.1.0.tar.gz
# optional, for the examples:
BiocManager::install("scRepertoire")
```

Requires a C++ compiler (Rtools on Windows, Xcode command-line tools on macOS).

## Examples

All examples use public data (scRepertoire's bundled `contig_list` /
`scRep_example`: 4 patients, peripheral blood and lung) or the package's
**synthetic** example cells.

### 1. scRepertoire output, `aes(fill = )` any column

```r
library(tcrhex)
library(scRepertoire)
library(ggplot2)

data("contig_list")
samples <- c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
combined <- combineTCR(contig_list, samples = samples)

cells <- read_screpertoire(combined, donor = sub("[BL]$", "", sample))
#> tcrhex cell table (scRepertoire): 27,185 cells, 22,526 paired (TRA+TRB), 4 donors
#>   unpaired, kept for top_public() but not for TCRdist: 798 alpha-only, 3,861 beta-only
#>   cells with 2+ chains of one locus (one kept): 2,965 TRA, 2,477 TRB
p <- hexrows(cells, aes(fill = sample), queries = top_public(cells, 10))
ggsave_hexrows("samples.png", p)
```

A Seurat or SingleCellExperiment object after `combineExpression()` works the
same way, and all of its meta data columns are available:
`read_screpertoire(sobj, donor = patient)` followed by
`hexrows(cells, aes(fill = seurat_clusters))`. A clone whose cells fall in
several clusters is coloured `"multiple"`. Pass
`summarise_fill = summarise_most_common` to use the most frequent cluster
instead.

![](man/figures/ex1_screp_sample.png)

### 2. 10x contigs + a metadata table, with `scale_fill_manual()`

```r
files <- c(P17B = "P17B/filtered_contig_annotations.csv",
           P17L = "P17L/filtered_contig_annotations.csv")  # ... one per sample
# (with the bundled example: files <- setNames(contig_list, samples))
meta <- data.frame(sample  = samples,
                   patient = sub("[BL]$", "", samples),
                   tissue  = ifelse(grepl("B$", samples), "peripheral blood", "lung"))

cells <- read_10x_contigs(files, meta = meta, by = "sample", donor = patient)
#> tcrhex cell table (10x contigs): 27,185 cells, 22,526 paired (TRA+TRB), 4 donors
#>   unpaired, kept for top_public() but not for TCRdist: 798 alpha-only, 3,861 beta-only
#>   cells with 2+ chains of one locus (one kept): 2,965 TRA, 2,477 TRB

p <- hexrows(cells, aes(fill = tissue), queries = top_public(cells, 10)) +
  scale_fill_manual(values = c("peripheral blood" = "#0072B2", lung = "#E69F00",
                               multiple = "#56B4E9"))
ggsave_hexrows("tissue.png", p)
```

![](man/figures/ex2_10x_tissue.png)

Free-text tissues can be grouped into compartments with `tissue_class()`. Use
`summarise_tissue` to turn a blood + lung clone into "blood + tissue", and the
colour-blind-safe tissue colour and pattern scales:

```r
hexrows(cells, aes(fill = tissue_class(tissue), pattern = tissue_class(tissue)),
        summarise_fill = summarise_tissue) +
  scale_fill_tissue() + scale_pattern_tissue()
```

### 3. Your own query TCR, coloured by HLA carrier status

```r
cells <- tcrhex_example_cells(n_donors = 40)   # SYNTHETIC cells with an `hla` column
p <- hexrows(cells, aes(fill = hla_carrier(hla, "A*02")),
             queries = c("my TCR" = "TRAV12-2:CAVNGGSQGNLIF:TRBV19:CASSIRSSYEQYF")) +
  scale_fill_manual(values = c(carrier = "#D55E00", "non-carrier" = "#0072B2",
                               "not typed" = "#d6d5ce"))
```

![](man/figures/ex4_synthetic_hla.png)

`hla_enriched()` finds, for each row, the allele most over-represented among
the donors near the query. Map it with
`aes(fill = hla_carrier(hla, .row_allele))`: every column of the query table
is available as `.row_<column>`. See the vignette.

## The pieces

| function | what it does |
|---|---|
| `read_screpertoire()`, `read_10x_contigs()`, `read_airr()`, `cells_from_df()` | one validated cell table, with a summary of dropped / unpaired / multi-chain cells |
| `tcr_clones()` | donor × TRAV+CDR3α × TRBV+CDR3β clones; every metadata column summarised by one rule |
| `tcrdist_paired()` | query × clone paired TCRdist (`d`, `d_a`, `d_b`) within a radius |
| `top_public()` | the n beta chains in most donors, each with its dominant alpha and concordance |
| `hexrows()` / `hexrows_data()` | the figure / the tidy data behind it (incl. the per-bin table) |
| `tissue_class()`, `summarise_tissue()`, `scale_fill_tissue()`, `scale_pattern_tissue()` | tissue compartments |
| `hla_carrier()`, `hla_enriched()` | HLA carrier status, per-row enriched allele |
| `ggsave_hexrows()` | saves with a page sized to rows × bins (0.072 in per hexagon unit) |

## Licences

MIT. The V-gene loop table is derived from tcrdist3's
`alphabeta_gammadelta_db.tsv` (MIT, © 2020 Koshlan Mayer-Blackwell). The logo
letters are traced from DejaVu Sans Mono Bold (Bitstream Vera licence). See
`inst/COPYRIGHTS`.
