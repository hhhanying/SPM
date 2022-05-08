# Used to generate array of transformation matrix T
T_generator <- function(k0, k1, nlabel){
  dg <- k0 + k1
  K <- nlabel * k0 + k1
  ntopic <- K
  Ts <- array(rep(0, dg * ntopic * nlabel), dim = c(nlabel, ntopic, dg))
  for (i in 1:nlabel){
    for (j in 1:k0) Ts[i, (i - 1) * k0 + j, j] = 1
    for (j in 1:k1) Ts[i, nlabel* k0 + j, k0 + j] = 1
  }
  
  Ts
}

# To accommodate more complicated cases, we allow hyperparameters for topic parameters to be different among dimensions / topics
# When the hyperparameter is a constant (same for all) or a vector (different only among dimensions / topics), this function is used to transform it to a matrix to adapt to the stan sampling
trans_to_matrix <- function(x, d1, d2, byfeature = TRUE){
  if (is.matrix(x)) return(x)
  if (length(x) == 1) return(matrix(rep(x, d1 * d2), nrow = d1))
  if (byfeature) return(matrix(rep(x, d1), nrow = d1, byrow = TRUE))
  return(matrix(rep(x, d2), nrow = d1, byrow = FALSE))
}




