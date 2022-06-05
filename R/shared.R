#' Given the model setting, generate the transformation matrices.
#' @export
T_generator <- function(k0, k1, nlabel){
  dg <- k0 + k1
  K <- nlabel * k0 + k1
  ntopic <- K
  Ts <- array(rep(0, dg * ntopic * nlabel), dim = c(nlabel, ntopic, dg))
  for (i in 1:nlabel){
    if(k0 > 0){
      for (j in 1:k0) Ts[i, (i - 1) * k0 + j, j] = 1
    }
    if(k1 > 0){
      for (j in 1:k1) Ts[i, nlabel* k0 + j, k0 + j] = 1
    }
  }
  
  Ts
}

#' Helpful function to summarize the sampling results.
#' @export
trans_to_matrix <- function(x, d1, d2, byfeature = TRUE){
  if (is.matrix(x)) return(x)
  if (length(x) == 1) return(matrix(rep(x, d1 * d2), nrow = d1))
  if (byfeature) return(matrix(rep(x, d1), nrow = d1, byrow = TRUE))
  return(matrix(rep(x, d2), nrow = d1, byrow = FALSE))
}

#' Generate transformed memberships
#' @export
U <- function(G, Y, Ts){
  N <- length(Y)
  K <- dim(Ts)[2]
  U <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N){
    U[i,] <- Ts[Y[i], ,] %*% G[i,]
  }
  U
}


