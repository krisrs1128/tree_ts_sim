
###############################################################################
# Simulations involving multivariate AR processes. Our interest here is in
# simulating time series on trees by enforcing MAR structure between parents
# and children.
###############################################################################

mar <- function(x0, el, A, eps) {
  stopifnot(length(dim(x0)) == 2)

  # root appears as parent, but not child
  root <- setdiff(el[, 1], el[, 2])
  X <- list()
  X[[root]] <- x0

  # proceed down tree
  parents <- unique(el[, 1])
  for (i in seq_along(parents)) {
    children <- el[which(el[, 1] == parents[i]), 2]

    for (j in seq_along(children)) {
      X[[children[j]]] <- A[[children[j]]] %*% X[[parents[i]]] + eps[[children[j]]]
    }
  }
  X
}
