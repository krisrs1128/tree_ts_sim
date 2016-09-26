
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


#' @title Normalize mixing matrices so children sum to identity
#' @description The constraint that keeps the children in a MAR process sum to
#' the parent (at each time point) is that the sum of the mixing matrices for
#' sibling nodes is the identity. Given a list of matrices that don't respect
#' this constraint, normalize them so that they do.
#' @param el [N x 2 matrix] A matrix whose rows are  edges between parents and
#' children.
#' @param A0 [list of p x p matrices] A list of matrices with names equal to the
#' entries in el.
#' @return A [list of p x p matrices] A new list of matrices that respect the
#' "sum to parent" constraint.
#' @examples
#' el <- matrix(
#'   c("1", "2", "1", "3", "2", "4", "2", "5", "3", "6", "3", "7"),
#'   ncol = 2,
#'   byrow = TRUE
#' )
#' A0 <- lapply(2:7, function(i) matrix(rnorm(100), 10, 10))
#' names(A0) <- 2:7
#' normalize_tree_mixing(el, A0)
normalize_tree_mixing <- function(el, A0) {
  A <- list()
  parents <- unique(el[, 1])

  for (i in seq_along(parents)) {
    children <- el[which(el[, 1] == parents[i]), 2]
    cur_A  <- A0[children]
    A <- c(A, sum_to_identity(cur_A))
  }
  A
}

#' @title Generate mixing matrices respecting tree structure
#' @description This wraps normalize_tree_mixing to generate random gaussian
#' matrices and modify them to sum to I.
#' @param el [N x 2 matrix] A matrix whose rows are  edges between parents and
#' children.
#' @param p [int] The dimension of each (square) A matrix to generate.
#' @param mean [numeric] The mean of the gaussian noise to generate in the first
#' step.
#' @param sd [numeric] The SD of the gaussian noise to generate in the first
#' step.
#' @return A [list of p x p matrices] A list of random matrices respecting the
#' "sum to I" across siblings constraint.
gaussian_tree_mixing <- function(el, p = 10, mean = 0, sd = 1) {
  root <- setdiff(el[, 1], el[, 2])
  v <- setdiff(unique(el), root)

  A0 <- lapply(v, function(i) { rnorm(p * p, p, mean, sd) })
  names(A0) <- v
  normalize_tree_mixing(el, A0)
}

#' @title Make list of matrices to sum to I
#' @description Given a list of initial matrices, generate a final list whose
#' sum is the identity. We use the cycling approach,
#' A_{1}(1) = I - A_{2}(0) - A_{3}(0) ...
#' then,
#' A_{2}(1) = I - A_{1}(1)- A_{3)(0)...
#' and hope it converges to something whose sum is the identity.
#' @param X_list [list of p x p matrices] Matrices to use as initialization to
#' this algorithm.
#' @param n_iter The number of cycles to run. Probably should check for
#' convergence, but this is cheaper.
#' @return X_list [list of p x p matrices] Matrices whose sum is the identity.
sum_to_identity <- function(X_list, n_iter = 10) {
  p <- nrow(X_list[[1]])
  for (iter in seq_len(n_iter)) {
    for (i in seq_along(X_list)) {
      X_list[[i]] <- diag(p) - Reduce("+", X_list[-i])
    }
  }
  X_list
}
