
#' LRE solution method based on Klein (2000, JEDC)
#'
#' This function solves for a linear policy function for the Linear Rational
#' Expectations model of \deqn{E (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}{E
#' (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}, where x and y are predetermined and
#' non-predetermined variables, respectively.
#'
#' @param A,E Square matrices of same size. \deqn{zE - A}{zE-A} is a regular
#'   pencil
#' @param nx The number of predetermined variables, \code{nx} is required
#'   by the algorithm.
#'
#' @return List of two functions (g, h), passed to \code{\link{simulate}}
#' @export
#'
#' @references Klein (2000, JEDC) \url{https://doi.org/10.1016/S0165-1889(99)00045-7}
#'
lre_auto_klein <- function(A, E, nx) {

  q <- qz(A, E)
  ns <- num_stable(q$ALPHA, q$BETA)

  bk_condition(nx, ns, silent = FALSE)

  n <- dim(A)[1]
  stbl <- 1:ns
  pre <- 1:nx
  npr <- (nx + 1):n

  S_ss <- q$S[stbl, stbl]
  T_ss <- q$T[stbl, stbl]
  Z_1s <- q$Z[pre, stbl]
  Z_2s <- q$Z[npr, stbl]

  g <- function(x) Z_2s %*% solve(Z_1s, x)
  h <- function(x) Z_1s %*% solve(S_ss, T_ss) %*% solve(Z_1s, x)

  list(g = g, h = h)
}

