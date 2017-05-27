

# Number of predetermined variables
num_predet <- function(A, nx, ny, x0) {
  if (!is.null(nx)) {
    return(nx)
  } else if (!is.null(ny)) {
    return(dim(A)[1] - ny)
  } else if (!is.null(x0)) {
    return(length(x0))
  } else {
    stop("Cannot determine the number of predetermined variables.")
  }
}


# Number of stable eigenvalues
num_stable <- function(alpha, beta = NULL) {
  if (is.null(beta)) {
    return(sum(abs(alpha) <= 1))
  } else {
    return(sum(abs(alpha) <= abs(beta)))
  }
}


#' Check Blanchard-Khan condition
#'
#' @param nx Number of predetermined variables.
#' @param ns Number of stable eigenvalues.
#' @param silent Boolean. If set FALSE, violation of BK condition causes an exception.
#'
#' @return Boolean. TRUE if nx = ns.
#' @export
#'
bk_condition <- function(nx, ns, silent = TRUE) {

  # If silent is FALSE, throw an error and stop execution.
  raise <- if (silent) {
    message
  } else {
    stop
  }

  if (nx == ns) {
    return(TRUE)
  } else if (nx > ns) {
      raise("nx > ns: Solution generally non-existent!")
  } else {
      raise("nx < ns: Indeterminacy case.")
  }
  FALSE
}


# Construct extended coefficients
extend_coeff <- function(A, E = NULL, B, Phi) {

  new_dim <- dim(A) + dim(Phi)

  if (is.null(E)) {
    E_ <- NULL
  } else {
    E_ <- matrix(0, new_dim[1], new_dim[2])
    E_[1:dim(Phi)[1], 1:dim(Phi)[2]] <- diag(1, ncol(Phi))
    E_[(dim(Phi)[1] + 1):nrow(E_), (dim(Phi)[2] + 1):ncol(E_)] <- E
  }

  A_ <- matrix(0, new_dim[1], new_dim[2])
  A_[1:dim(Phi)[1], 1:dim(Phi)[2]] <- Phi
  A_[(dim(Phi)[1] + 1):nrow(A_), 1:dim(Phi)[2]] <- B
  A_[(dim(Phi)[1] + 1):nrow(A_), (dim(Phi)[2] + 1):ncol(A_)] <- A

  return(list(A = A_, E = E_))
}


#' LRE solution method under AR inputs
#'
#' This function solves for a linear policy function for the Linear Rational
#' Expectations model of \deqn{E (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t}) +
#' Bu_{t}}{E (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t}) + Bu_{t}}, where x and y are
#' predetermined and non-predetermined variables, respectively. The input term
#' follows AR process.
#'
#' @param A,E Square matrices of same size. \deqn{zE - A}{zE-A} is a regular
#'   pencil. E = NULL is understood as E = I.
#' @param B coefficient matrix for input term
#' @param Phi matrix. Inputs are assumed to follow VAR process u_{t+1} = Phi %*%
#'   u_{t} + epsilon_{t+1}
#' @param nx,ny,x0 The number of predetermined variables, \code{nx}, is required
#'   by the algorithm. It is passed by \code{nx} or computed from the number of
#'   the non-predetermined \code{ny}, or the initial vector \code{x0}
#'
#' @return List of two functions (g, h), passed to \code{\link{simulate}}
#' @export
#'
lre_ar <- function(A, E = NULL, B, Phi, nx = NULL, ny = NULL, x0 = NULL) {

  nx <- num_predet(A, nx, ny, x0) + ncol(Phi)

  coeff <- extend_coeff(A, E, B, Phi)
  A_ <- coeff$A
  E_ <- coeff$E

  if (is.null(E_)) {
    lre_auto_bk(A_, nx = nx)
  } else {
    # Check if A_ and E_ have same size
    stopifnot(all(dim(A_) == dim(E_)))

    lre_auto_klein(A_, E_, nx = nx)
  }
}


#' LRE solution method
#'
#' This function solves for a linear policy function for the Linear Rational
#' Expectations model of \deqn{E (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}{E
#' (x_{t+1}, y_{t+1}) = A (x_{t}, y_{t})}, where x and y are predetermined and
#' non-predetermined variables, respectively. When E is NULL, it is understood
#' as E = I (identity matrix) and Schur decomposition of A is used to compute
#' the solution. When E is not NULL, the generalized Schur decomposition of (A,
#' E) is used.
#'
#' @param A,E Square matrices of same size. \deqn{zE - A}{zE-A} is a regular
#'   pencil. E = NULL is understood as E = I.
#' @param nx,ny,x0 The number of predetermined variables, \code{nx}, is required
#'   by the algorithm. It is passed by \code{nx} or computed from the number of the
#'   non-predetermined \code{ny}, or the initial vector \code{x0}
#'
#' @return List of two functions (g, h), passed to \code{\link{simulate}}
#' @export
#'
#' @references Klein (2000, JEDC)
#'   \url{https://doi.org/10.1016/S0165-1889(99)00045-7}
#'
#' @seealso \code{\link{lre_auto_klein}}, \code{\link{lre_auto_bk}}
#'
lre_auto <- function(A, E = NULL, nx = NULL, ny = NULL, x0 = NULL) {

  # Check if A is square
  stopifnot(dim(A)[1] == dim(A)[2])

  if (is.null(E)) {
    lre_auto_bk(A, nx = num_predet(A, nx, ny, x0))
  } else {
    # Check if A and E have same size
    stopifnot(all(dim(A) == dim(E)))

    lre_auto_klein(A, E, nx = num_predet(A, nx, ny, x0))
  }
}
