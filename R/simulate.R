
#' Simulate state space model
#'
#' This function simulates a state space model with a state equation,
#' x_{t+1} = h(x_{t}), and output equation, y_{t} = g(x_{t}), from initial condition
#' x_{0}. In macroeconomics, x is said to be the vector of predetermined variables and
#' y the vector of non-predetermined variables.
#'
#' @param g maps predetermined to non-predetermined
#' @param h maps predetermined forward in time
#' @param x0 initial values for predetermined variables
#' @param t length of simulation
#'
#' @return matrix comprising of simulation output
#' @export
#'
simulate <- function(g = NULL, h, x0, t) {
  if (is.null(g)) {
    simulate_ssr_no_output(h, x0, t)
  } else {
    simulate_ssr(g, h, x0, t)
  }
}


#' @rdname simulate
simulate_ssr_no_output <- function(h, x0, t) {
  n <- length(x0)
  out <- matrix(0, t, n)
  for (i in 1:(t - 1)) {
    out[i + 1, ] <- h(out[i, ])
  }
  out
}

#' @rdname simulate
simulate_ssr <- function(g, h, x0, t) {
  n1 <- length(x0)     # Number of predetermined variables
  n2 <- length(g(x0))  # Number of non-predetermined variables

  pre <- 1:n1
  npr <- (n1 + 1):(n1 + n2)

  out <- matrix(0, t, n1 + n2)  # Zero matrix for simulation output

  out[1, pre] <- x0     # Initial Condition
  out[1, npr] <- g(x0)  # Eq. (2.1)

  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i, pre])         # Eq. (2.2)
    out[i + 1, npr] <- g(out[i + 1, pre])  # Eq. (2.1)
  }
  out
}
