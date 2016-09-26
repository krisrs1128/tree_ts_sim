
###############################################################################
# Simulations involving multivariate AR processes. Our interest here is in
# simulating time series on trees by enforcing MAR structure between parents
# and children.
###############################################################################

#' @title MAR process, given all parameters
#' @description This propagates a MAR process given the mixing matrices A and
#' innovations eps. While this is more memory intensive than computing eps on
#' the fly, this lets us experiment with different noise structures more easily.
#' @param x0 [p x k matrix] A vector representing the time series for the root
#' node.
#' @param el [N x 2 matrix] A matrix whose rows are  edges between parents and
#' children.
#' @param A [length V list of p x p matrices] A list whose names are the names
#' of vertices in el, and whose values are the mixing matrices associated with
#' each node.
#' @param eps [length V list of p x k vectors] A list whose names are the names
#' of vertices in el, and whose values are the innovation matrices associated
#' with each node.
#' @return X [list of matrices] A list giving the simulated MAR process, whose
#' v^th element is the simulated matrix for node v.#'
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
