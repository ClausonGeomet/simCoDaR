# Adapted from CoDA R/composition_study.R at
# c2184f44994a7c3f2d4f36ec74c749e0267a4127; GPL-3.
# Origin: Research Team; history: Matt Clauson. See inst/PROVENANCE.md.
.ilr_basis <- function(d) {
  helmert <- stats::contr.helmert(d)
  t(sweep(helmert, 2L, sqrt(colSums(helmert^2)), "/"))
}
.to_ilr <- function(x, basis) {
  if (any(!is.finite(x)) || any(x <= 0) || ncol(x) != ncol(basis))
    stop("ILR input must be finite, positive and match the basis.")
  z <- log(x) %*% t(basis)
  if (any(!is.finite(z))) stop("Nonfinite ILR center.")
  z
}
.from_ilr <- function(z, basis) {
  if (any(!is.finite(z)) || ncol(z) != nrow(basis))
    stop("ILR coordinates must be finite and match the basis.")
  logx <- z %*% basis
  if (any(!is.finite(logx))) stop("Nonfinite inverse ILR calculation; reduce spread.")
  logx <- sweep(logx, 1L, apply(logx, 1L, max), "-")
  if (any(!is.finite(logx))) stop("Nonfinite stabilized ILR calculation.")
  x <- exp(pmax(logx, -700))
  x / rowSums(x)
}
