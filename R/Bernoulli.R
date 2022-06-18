#' Generate data from given parameters.
#'
#' @param P: the probability for each topic.
#' @param logits: the logit for each topic. If not provided, it will be calculated from P.
#' @param w: the probability for each class. Default to be NULL, i.e., uniform distribution.
#' @return A list res containing: 1) X: the document; 2) Y: the labels; 3) G: the untransformed memberships; 4) U: the transformed memberships; 5) PX: the success probability for each data point and 6) the logit of the success probability for each data point. 
#' @export
document_generator_Bernoulli <- function(a, rho, Ts, N, P = NULL, logits = NULL, w = NULL, seed = NULL){
  if (!is.null(seed)){ # if seed is provided, set seed
    set.seed(seed)
  }  
  
  if (!is.null(logits)){ # if logits are not provided, calculate it from P
    logits <- log(P / (1 - P))
  }   
  
  # get necessary parameters
  nlabel <- dim(Ts)[1]
  d <- dim(logits)[2]
  K <- dim(Ts)[2]
  
  Y <- sample(x = nlabel, size = N, replace = TRUE, prob = w) # sample labels
  G <- gtools::rdirichlet(n = N, alpha = a * rho) # sample memberships
  U <- matrix(NA, nrow = N, ncol = K) # transform memberships
  for (i in 1:N){
    U[i,] <- Ts[Y[i], ,] %*% G[i,]
  }  
  
  logitX <- U %*% logits # get parameters for generated distributions
  PX <- exp(logitX) / (1 + exp(logitX))
  
  X <- rbinom(n = N * d, size = 1, prob = as.vector(PX))
  X <- matrix(X, nrow = N, byrow = FALSE)
  
  list(X = X, Y = Y, G = G, U = U, PX = PX, logitX = logitX) 
}