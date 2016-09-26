context("mar")

test_that("Generates sum given mixture of identities.", {
  el <- matrix(
    c("1", "2", "1", "3", "2", "4", "2", "5", "3", "6", "3", "7"),
    ncol = 2,
    byrow = TRUE
  )

  v <- unique(as.character(el))
  A <- vector(length = length(v), mode = "list")
  eps <- vector(length = length(v), mode = "list")
  names(A)  <- v
  names(eps) <- v

  p <- 10
  x0 <- matrix(rnorm(p), p, 1)

  for (i in seq_along(v)) {
    A[[v[i]]] <- .5 * diag(p) # binary tree
    eps[[v[i]]] <- rep(0, p)
  }

  X <- list(
    "1" = x0,
    "2" = 0.5 * x0,
    "3" = 0.5 * x0,
    "4" = 0.25 * x0,
    "5" = 0.25 * x0,
    "6" = 0.25 * x0,
    "7" = 0.25 * x0
  )

  expect_equal(mar(x0, el, A, eps), X)
})
