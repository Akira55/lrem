
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
#' discrete-time) are collected at the upper-left corner.
#'
#' @param A,B Sqaure matrices of same size. (A, B) is a regular pencil.
#' @return QZ decomposition.
#'
qz <- function(A, B) {
  q <- QZ::qz(A, B)

  # Eigenvals in the Closed Unit Disk are collected at the upper-left corner
  cud <- abs(q$ALPHA) < abs(q$BETA)
  ret <- QZ::qz.dtgsen(q$S, q$T, q$Q, q$Z, select = cud)
  return(ret)
}
