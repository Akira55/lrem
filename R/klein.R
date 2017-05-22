
qz_klein <- function(A,E) {
  #' QZ decompose of pair matrixes A and E
  #'
  #' @param A a legurar matrix that is the same size as E
  #' @param E a legurar matrix that is the same size as E
  #'
  #' @importFrom QZ qz
  #' @importFrom QZ qz.dtgsen
  #'
  # Args:
  #    A and E: a pair of the same size matrixes
  # Returns
  #  Omega_x, Omega_u, Omega_y, Psi_x, Psi_y: coefficient matrixes
  #
  # QZ decomposition
  ret0 <- qz(A, E)
  #Decomposition elements
  alpha = ret0$ALPHA
  beta = ret0$BETA
  S = ret0$S
  T = ret0$T
  Q = ret0$Q
  Z = ret0$Z
  # Making eignvalue according to the order in Klein(2000)
  select  <- abs(alpha) < abs(beta)
  ret <- qz.dtgsen(S, T, Q, Z, select = select)
  eig <- ret$ALPHAR / ret$BETA

  decom <- list(S="matrix", T="matrix", Q="matrix", Z="matrix", eig <- "array")
  decom$S <- ret$S
  decom$T <- ret$T
  decom$Q <- ret$Q
  decom$Z <- ret$Z
  decom$eig <- eig

  return(decom)
}

klein_coefficients_from <- function(decom){
  #' Making parameters accrding Klein(200)
  #'
  #' @param decom a list object made by qz_klein()
  #' @export
  #
  # Return:
  #   coefficients: a list object which contains Omega_x, Omega_u, Omega_y, Psi_x and Psi_y
  #
  # making block matrixes of each decomposit element
  Z <- make_blocks_of(decom$Z)
  T <- make_blocks_of(decom$S)#'講義ノートに合わせるためにS,Tを交換しておく。
  S <- make_blocks_of(decom$T)
  B <- matrix(c(0. , 0. , 0. , 1.0 ), nrow=4, ncol=1)
  # Making block matrixes of C
  C <- t(decom$Q) %*% B
  nrowC <- nrow(C)
  rnC_half <- nrow(C)/2
  Cs <- matrix(C[1:rnC_half,])
  Cu <- matrix(C[(rnC_half+1):nrowC,])
  # To define parameters
  Z1s <- Z$a11
  Z1u <- Z$a12
  Z2s <- Z$a21
  Z2u <- Z$a22
  Tss <- T$a11
  Tsu <- T$a12
  Tuu <- T$a22
  Sss <- S$a11
  Ssu <- S$a12
  Suu <- S$a22
  iSss <- solve(Sss)
  iTuu <- solve(Tuu)
  iZ1s <- solve(Z1s)

  Omega_x <- Z1s %*% iSss %*% Tss %*% iZ1s
  Omega_u <- Z1s %*% iSss %*% (Cs - (Tsu - Tss %*% iZ1s %*% Z1u) %*% iTuu %*% Cu)
  Omega_y <-  (Z1u - Z1s %*% iSss %*% Ssu) + Z1s %*% iSss %*% (Tsu - Tss %*%  iZ1s  %*%  Z1u) %*%  iTuu  %*%  Suu
  Psi_x <- Z2s %*%  iZ1s
  Psi_y <- Z2u - Z2s %*% iZ1s %*% Z1u

  coefficients <- list(Omega_x ="matrix", Omega_u = "matrix" , Omega_y = "matrix", Psi_x = "matrix", Psi_y = "matrix")
  coefficients$Omega_x <- Omega_x
  coefficients$Omega_u <- Omega_u
  coefficients$Omega_y <- Omega_y
  coefficients$Psi_x <-  Psi_x
  coefficients$Psi_y <- Psi_y
  coefficients$Tuu <- Tuu
  coefficients$iTuu <- iTuu
  coefficients$Cu <- Cu
  coefficients$Suu <- Suu

  return(coefficients)
}



make_blocks_of <- function(A) {
  #' making block matrix of a square matrix A
  #'
  #' @param A a square matrix A
  #'
  # ATTENTION: We do not require a matrix A to be square to make blocks of it
  #            But in this time we just consider a square matrix A
  # Args:
  #   A: a matrix
  #
  # Returns:
  #   block: a list that contains block matrixes of A
  #
  # Error handling
  if(nrow(A) != ncol(A)) {
    return("this matrix is not square")
  }

  # take a mat size to make block matrix
  size <- nrow(A)
  n_first <- 0:(size/2)
  n_latter <- (size/2+1):size

  # make block as follow
  # A =([a11 a12]
  #     [a21 a22])

  a11 <- A[n_first, n_first]
  a21 <- A[n_latter, n_first]
  a12 <- A[n_first, n_latter]
  a22 <- A[n_latter, n_latter]

  block <- list(a11 = "matrix", a21 = "matrix", a12 = "matrix", a22 = "matrix")
  block$a11 <- a11
  block$a21 <- a21
  block$a12 <- a12
  block$a22 <- a22

  return(block)
}



yu <- function(u, decom) {
  #' caluculating y
  #'
  #' @param u an array that contains shocks
  #' @param decom a list object made by qz_klein()
  #'
  # Returns:
  #   y: an array that contains a  caluculation result
  coeff <- klein_coefficients_from(decom)

  # a block matrix right downward of QZ decomposed matrix S
  Suu <- coeff$Suu

  # a block matrix right downward of QZ decomposed matrix T
  iTuu <- coeff$iTuu
  Cu <- coeff$Cu

  # caluculate y
  length_n <- dim(u)
  y <- array(1:length_n*2, dim = c(2, length_n+1))
  u0 <- reverse_1dim_array(u)
  y[,1] <- matrix(c(0,0), c(2,1)) #'the assumption that y_{T+1} = 0

  # y_T1 <- iTuu %*% Suu %*% y_T1
  # y1   <- matrix(c(0,0), c(2,1)) #at first, pervious value y1 = 0
  for (i in 1:length_n) {
    yplus <- iTuu %*% Suu %*% y[,i]
    y1 <- yplus - iTuu %*% Cu * u0[i]
    y[,i+1] <- y1
  }
  y <- reverse_2dim_array(y)
  return(y)
}


reverse_1dim_array <- function(original) {
  #' returns reversed array of original array
  #'
  #' @param original an array to be reversed
  #'
  # Returns:
  #   reversed: an array which contains elements in the reversed order of original
  #
  n_size <- dim(original)
  reversed <-as.array(1:n_size)
  for (i in 1:n_size) {
    reversed[n_size-i+1] <- original[i]
  }
  return(reversed)
}

reverse_2dim_array <- function(original) {
  #' returns reversed array of original array
  #'
  #' @param original an array to be reversed
  #'
  # Returns:
  #   reversed: an array which contains elements in the reversed order of original
  #
  n_size <- dim(original)
  reversed <-array(1:n_size[2]*n_size[1], dim = c(n_size[1], n_size[2]))
  for (i in 1:n_size[2]) {
    reversed[,n_size[2]-i+1] <- original[,i]
  }
  return(reversed)
}


h <- function(x1, u, coefficients, y) {
  #' h a function that computes the predetermined of the next period
  #' @param x1 the predetermied variable
  #' @param u shocks
  #' @param coefficients coeffients matrixes that is provid by klein_coefficients_from()
  #' @param y an valu that from an array computed by yu()
  # Returns:
  #    x1: the non-predetermined

  Omega_x <- coefficients$Omega_x
  Omega_u <- coefficients$Omega_u
  Omega_y <- coefficients$Omega_y


  #x1[,i+1] <- Omega_x %*% x1[,i] + Omega_u * u[i] + Omega_y %*% y[,i+1]
  x1 <- Omega_x %*% x1 + Omega_u * u + Omega_y %*% y
  return(x1)
}


g <- function(x1, coefficients, y) {
  #' a function that computes the non-predetermined from the predetermined
  #'
  #' @param x1 the predetermied variable
  #' @param coefficients coeffients matrixes that is provid by klein_coefficients_from()
  #' @param y an valu that from an array computed by yu()
  #'
  # Returns:
  #    x1: the non-predetermined

  Psi_x <- coefficients$Psi_x
  Psi_y <- coefficients$Psi_y

  #x2[,i+1] <- Psi_x %*% x1[,i] + Psi_y %*% y[,i]
  x2 <- Psi_x %*% x1 + Psi_y %*% y
  return(x2)
}


simulate <- function(A, E, g, h, x0, t, u) {
  #' simulate a model using function that decides pre and non-pre determinant variables
  #'
  #' @param A a legular matrix that is the same size as E
  #' @param E a legurar matrix that is the same size as E
  #' @param g a function that computes the non-predetermined from the predetermined
  #' @param h a function that computes the predetermined of the next period
  #' @param x0 an initial condition
  #' @param t a number of times to caluculate
  #' @param u an array that contains shocks
  #'
  #' @export
  #'
  # Return:
  #     out: an matrix that contains simulation result

  #関数g,h に必要
  decom <- qz_klein(A,E)
  coefficients <- klein_coefficients_from(decom)
  y <- yu(u, decom)

  n1 <- length(x0)
  n2 <- length(g(x0, coefficients, y[,1]))

  #the number of pre and non-predeterminant variables
  pre <- 1:n1
  npr <- (n1 + 1):(n1 + n2)
  # an array to store them
  out <- matrix(0, t, n1 + n2)  # Zero matrix for simulation output

  out[1, pre] <- x0     # Initial Condition
  out[1, npr] <- g(x0, coefficients, y[,1])  # Eq. (2.1)

  for (i in 1:(t - 1)) {
    out[i + 1, pre] <- h(out[i,pre], u[i], coefficients, y[, i+1])
    out[i + 1, npr] <- g(out[i+1,npr], coefficients, y[,i])
  }
  return(out)
}
