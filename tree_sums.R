
#' @title Sum of Gaussian Processes
#' @description Each node is associated with a GP, and each leaf is the sum of
#' its ancestors.
#' @param
#' @return
#' @importFrom simData gp_data
#' @examples
#' el <- matrix(
#'   c("1", "2", "1", "3", "2", "4", "2", "5", "3", "6", "3", "7"),
#'   ncol = 2,
#'   byrow = TRUE
#' )
#' times <- sort(runif(100))
#' tree_data <- tree_sum(gp_data, el, times = times, bandwidth = 0.1)
#' mtree_data <- melt(tree_data)
#' mtree_data$times <- times[mtree_data$Var1]
#' ggplot(mtree_data) +
#'   geom_point(aes(x = times, y = value, col = L1))
tree_sum <- function(node_fun, el, ...) {
  args <- list(...)
  x0 <- node_fun(args$times, args$bandwidth)
  p <- length(x0)

  nodes <- unique(as.character(el))
  eps <- replicate(
    length(nodes),
    node_fun(args$times, args$bandwidth),
    simplify = FALSE
  )
  names(eps) <- nodes

  A <- replicate(
    length(nodes),
    diag(p),
    simplify = F
  )
  names(A) <- nodes

  mar(x0, el, A, eps)
}
