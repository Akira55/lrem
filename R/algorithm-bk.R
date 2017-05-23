
#' LRE solution method based on Blanchard and Kahn (1980, ECTA)
#'
#' This function solves for a linear policy function for the Linear Rational
#' Expectations model of \deqn{(x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}{
#' (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}, where x and y are predetermined and
#' non-predetermined variables, respectively.
#'
#' @param A Square matrix
#' @param nx The number of predetermined variables, \code{nx} is required
#'   by the algorithm.
#'
#' @return List of two functions (g, h), passed to \code{\link{simulate}}
#' @export
#'
lre_auto_bk <- function(A, nx) {

  q <- schur(A)
  ns <- num_stable(q$W)

  bk_condition(nx, ns, silent = FALSE)

  n <- dim(A)[1]
  stbl <- 1:ns
  pre <- 1:nx
  npr <- (nx + 1):n

  Q_1s <- q$Q[pre, stbl]
  Q_2s <- q$Q[npr, stbl]
  A_11 <- A[pre, pre]
  A_12 <- A[pre, npr]

  g <- function(x) Q_2s %*% solve(Q_1s, x)
  h <- function(x) A_11 %*% x + A_12 %*% g(x)

  list(g, h)
}
