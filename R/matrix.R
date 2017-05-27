
#' Schur decomposition
#'
#' Real Schur decomposition. Weakly stable eigenvalues (in discrete-time) are
#' collected at the upper-left corner.
#'
#' @param A Sqaure matrix.
#' @return Schur decomposition of \code{A}.
#'
schur <- function(A) {
  q <- QZ::qz(A)

  # Eigenvals in the Closed Unit Disk are collected at the upper-left corner
  QZ::qz.dtrsen(q$T, q$Q, abs(q$W) <= 1)
}


#' QZ decomposition
#'
#' Real QZ decomposition. Weakly stable generalized eigenvalues (in
#' discrete-time) are collected at the upper-left corner. A pencil (zE - A) is
#' converted to (zS - T).
#'
#' @param A,E Square matrices of same size.
#' @return QZ decomposition. See Note for additional information.
#'
#' @note The matrix labels for the returned values follow those of Klein (2000)
#'   \url{https://doi.org/10.1016/S0165-1889(99)00045-7}
#'   and Kenji Sato's lecture slides
#'   \url{https://www.kenjisato.jp/teaching/ed/2017/}. In particular, \code{S}
#'   and \code{T} returned by this function have reversed roles compared to
#'   \code{QZ::qz}
#'
qz <- function(A, E) {
  q <- QZ::qz(A, E)

  # Eigenvals in the Closed Unit Disk are collected at the upper-left corner
  cud <- abs(q$ALPHA) < abs(q$BETA)
  ret <- QZ::qz.dtgsen(q$S, q$T, q$Q, q$Z, select = cud)

  s <- ret$S
  t <- ret$T
  ret$S <- t
  ret$T <- s

  return(ret)
}
