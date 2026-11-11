fx <- function(...) test_path("fixtures", ...)
core <- function(x) {
  x <- as.data.frame(x)[, c("cell", "donor", "va", "ca", "vb", "cb", "tissue")]
  x <- x[order(x$cell), ]
  rownames(x) <- NULL
  x
}
