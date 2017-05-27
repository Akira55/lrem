
#' Simulate state space model
#'
#' This function simulates a state space model with a state equation, x_{t+1} =
#' h(x_{t}) + e_{t+1}, and output equation, y_{t} = g(x_{t}), from initial
#' condition x_{0}. In macroeconomics, x is said to be the vector of
#' predetermined variables and y the vector of non-predetermined variables.
#'
#' @param g maps predetermined to non-predetermined
#' @param h maps predetermined forward in time
#' @param x0 initial values for predetermined variables
#' @param t integer, length of simulation
#' @param e matrix, shock sequence
#'
#' @return matrix comprising of simulation output
#' @export
#'
#' @seealso \code{\link{impulse}}
simulate <- function(g = NULL, h, x0, t = NULL, e = NULL) {
  if (is.null(g)) {
    simulate_state_eqn(h, x0, t, e)
  } else {
    simulate_ssr(g, h, x0, t, e)
  }
}


#' @rdname simulate
simulate_state_eqn <- function(h, x0, t, e) {

  # Neaten shock sequence
  . <- adjust_input(e, t, x0)
  e <- .$e
  t <- .$t

  # Initial Condition
  n <- length(x0)
  out <- matrix(0, t, n)
  out[1, ] <- x0

  for (i in 1:(t - 1)) {
    out[i + 1, ] <- h(out[i, ]) + e[i, ]
  }
  out
}

#' @rdname simulate
simulate_ssr <- function(g, h, x0, t, e) {
  # Neaten shock sequence
  . <- adjust_input(e, t, x0)
  e <- .$e
  t <- .$t

  # Predetermined variables
  n1 <- length(x0)
  pre <- 1:n1

  # Non-predetermined variables
  n2 <- length(g(x0))
  npr <- (n1 + 1):(n1 + n2)

  # Initial conditions
  out <- matrix(0, t, n1 + n2)
  out[1, pre] <- x0

  # Populate simulation output
  out[1, npr] <- g(x0)
  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i, pre]) + e[i, ]
    out[i + 1, npr] <- g(out[i + 1, pre])
  }
  out
}


#' Compute impulse response
#'
#' @param g maps predetermined to non-predetermined
#' @param h maps predetermined forward in time
#' @param x0 initial values for predetermined variables. For linear case, they
#'   are zero.
#' @param t integer, length of simulation
#' @param e1 impulse shock
#'
#' @return matrix comprising of simulation output
#' @export
#'
#' @seealso \code{\link{simulate}}
#'
impulse <- function(g, h, x0, t, e1) {
  e <- matrix(e1, nrow = 1)
  simulate(g, h, x0, t, e)
}


# Modify shock sequence, e  -----

adjust_input <- function(e, t, x0) {
  . <- adjust_input_length(e, t)

  e <- adjust_input_width(.$e, x0)
  list(e = e, t = .$t)
}


# Modify e so as to satisfy nrow(e) == t - 1
adjust_input_length <- function(e, t) {

  # If both e and t are NULL, then raise an error
  if (is.null(e) && is.null(t)) stop("Neither e nor t is provided.")

  # If t is NULL, then nrow(e) + 1 is used as t
  if (is.null(t)) return(list(e = e, t = nrow(e) + 1))

  # If simulation length, t, is longer than nrow(e),
  # fill missing inputs with zero
  if (nrow(e) < t - 1) {
    enew <- matrix(0, t - 1, ncol(e))
    enew[1:nrow(e), ] <- e
    return(list(e = enew, t = t))
  }
}


# Modify e so as to satisfy ncol(e) == length(x0)
adjust_input_width <- function(e, x0) {

  # This is not supposed to happen.
  if (is.null(e)) {
    return(matrix(0, 1, length(x0)))
  }

  nx <- length(x0)
  ne <- ncol(e)

  if (nx < ne) {
    message("Shock vector has more components than state vector.
             Information will be lost.")
    enew <- e[, 1:nx]
  } else if (nx > ne) {
    message("Shock vector has fewer components than state vector.")
    enew <- cbind(e, matrix(0, nrow(e), nx - ne))
  } else {
    enew <- e
  }
  enew
}



