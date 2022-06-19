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
  
  if (is.null(logits)){ # if logits are not provided, calculate it from P
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

#' Simulate data from given memberships
#' 
#' @return A matrix X.
#' @export
SPM_simulator_Bernoulli <- function(P = NULL, logits = NULL, Ts = NULL, G = NULL, Y = NULL, U = NULL, seed = NULL){
  if (!is.null(seed)){ # if seed is provided, set seed
    set.seed(seed)
  }  
  
  if (is.null(logits)){ # if logits are not provided, calculate it from P
    logits <- log(P / (1 - P))
  } 
  
  if(is.null(U)){ # if transformed memberships are not given, calculate it from the untransformed memberships and labels
    K <- dim(Ts)[2]
    N <- length(Y)
    
    U <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N){
      U[i,] <- Ts[Y[i], ,] %*% G[i,]
    }    
  } else{
    N <- dim(U)[1]
  }
  
  d <- dim(logits)[2]
  
  logitX <- U %*% logits # get parameters for generated distributions
  PX <- exp(logitX) / (1 + exp(logitX))
  
  X <- rbinom(n = N * d, size = 1, prob = as.vector(PX))
  X <- matrix(X, nrow = N, byrow = FALSE)
  
  X
}

#' Get the parameters for the generated distribution for each data point
#' 
#' @export
get_parameters_Bernoulli <- function(P = NULL, logits = NULL, G = NULL, Ts = NULL, Y = NULL, U = NULL){
  if (is.null(logits)){ # if logits are not provided, calculate it from P
    logits <- log(P / (1 - P))
  } 
  
  if (is.null(U)){ # if given membership need to be transformed
    N <- length(Y)
    K <- dim(logits)[1]
    U <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N){
      U[i,] <- Ts[Y[i], ,] %*% G[i,]
    }
  } 
  
  logitX <- U %*% logits 
  PX <- exp(logitX) / (1 + exp(logitX))
  
  list(U = U, PX = PX, logitX = logitX)
}

#' Fit the SPM model for Bernoulli distribution.
#' 
#' @param save_trace: default to be FALSE. If TRUE will also return the trace.
#' @export
SPM_training_Normal <- function(X, Y, Ts, b, alpha, alpha_p, beta_p, VI = FALSE, ntrace = 1000, nchain = 1, nskip = 2, seed = 1, save_trace = FALSE){
  SPM_Bernoulli_stancode <- "
  data {
    int<lower=0> dg;  // dim(membership)
    int<lower=0> N;  // size of training set
    int<lower=0> nlabel;
    int<lower=0> ntopic;
    int<lower=0> d; //dim of data
    
    int<lower=0, upper=1> X[N, d];
    int<lower=1> Y[N];         // label, start from 1 
    
    matrix[ntopic, dg] T[nlabel];
    
    real<lower=0> b;
    vector<lower=0>[dg] alpha;
    
    matrix[ntopic, d] alpha_p;
    matrix[ntopic, d] beta_p;
  }
  parameters {
    real<lower=0> a;
    simplex[dg] rho;
    
    simplex[dg] G[N];
    
    matrix<lower=0, upper=1>[ntopic, d] P;
  }
  model{
    matrix[ntopic, d] logits;
    vector[ntopic] U[N];
    row_vector[d] logit_X[N];
    
    a ~ exponential(b);
    rho ~ dirichlet(alpha);
  
    // topics
    for (i in 1:ntopic){
      P[i] ~ beta(alpha_p[i], beta_p[i]);
      logits[i] = logit(P[i]);
    }
    
    for(i in 1:N){
      G[i] ~ dirichlet(a * rho);
      U[i] = T[Y[i]] * G[i];
      logit_X[i] = U[i]' * logits;
      X[i] ~ bernoulli_logit(logit_X[i]);
    }
  }
  "
  
  N <- length(Y)
  nlabel <- dim(Ts)[1]
  dg <- dim(Ts)[3]
  ntopic <- dim(Ts)[2]
  d <- dim(X)[2]
  K <- ntopic
  
  dat_fit <- list(
    dg = dg, N = N, nlabel = nlabel, ntopic = ntopic, d = d,
    X = X, Y = Y, T = Ts,
    b = b, alpha = alpha,
    alpha_p = trans_to_matrix(alpha_p, ntopic, d), 
    beta_p = trans_to_matrix(beta_p, ntopic, d)
  )
  
  if(VI){
    model <- rstan::stan_model(model_code = SPM_Bernoulli_stancode)
    fit_train <- rstan::vb(model, data = dat_fit, seed = seed)
  } else{
    fit_train <-rstan::stan(model_code = SPM_Bernoulli_stancode,
                            data = dat_fit,
                            chains = nchain,
                            iter = ntrace,
                            seed = seed)    
  }
  
  trace <- as.matrix(fit_train)
  
  res <- list()
  
  if(save_trace){
    res[["trace"]] <- trace
  }
  
  # save results
  nsample <- dim(trace)[1]
  nsave <- nsample %/% nskip
  index_save <- (1:nsave) * nskip
  trace <- trace[index_save,]
  
  paras <- c("G", "P")
  d1s <- c(N, K)
  d2s <- c(dg, d)
  for (i in 1:length(paras)){
    para <- paras[i]
    d1 <- d1s[i]
    d2 <- d2s[i]
    tem <- matrix(NA, nrow = d1, ncol = d2)
    for (i1 in 1:d1){
      for(i2 in 1:d2){
        varname <- sprintf("%s[%s,%s]", para, i1, i2)
        tem[i1, i2] <- mean(trace[, varname])
      }
    } 
    res[[para]] <- tem
  }  
  
  res[["a"]] <- mean(trace[, "a"])
  
  tem <- rep(0, dg)
  for (i in 1:dg){
    varname <- sprintf("rho[%s]", i)
    tem[i] <- mean(trace[, varname])
  }
  res[["rho"]] <- tem
  
  res
}

#' Predict the labels from given topics.
#' 
#' Use MCMC to calculate the posterior probability for labels and choose the MAP.
#' 
#' @param w: the prior distribution of the labels. Default to be NULL, meaning uniform distribution.
#' @param nsample: the sample size to draw for estimating the posterior probability.
#' @export
SPM_predicting_Bernoulli <- function(X, a, rho, Ts, P = NULL, logits = NULL, w = NULL, nsample = 1000, seed = NULL){
  if (!is.null(seed)){ # if seed is provided, set seed
    set.seed(seed)
  }  
  
  if (is.null(logits)){ # if logits are not provided, calculate it from P
    logits <- log(P / (1 - P))
  } 
  
  if (is.null(w)) { # if no information about the prior of labels, treat it as uniform
    w <- rep(1, nlabel) 
  }
  
  nlabel <- dim(Ts)[1]
  N <- dim(X)[1]  
  
  probs <- matrix(NA, nrow = N, ncol = nlabel) # stores the posterior probability of labels
  Y <- rep(0, N)
  
  for(i in 1:N){ # for each data point
    G <- gtools::rdirichlet(n = nsample, alpha = a * rho) # draw G
    
    for(y in 1:nlabel){ # calculate log posterior for each possible label
      U <- G %*% t(Ts[y,,]) # nsample * K
      x <- matrix(rep(X[i,], nsample), nrow = nsample, byrow = TRUE) # each row is X[i]
      logitX <- U %*% logits 
      logP <- x * logitX - log(1 + exp(logitX)) # logP[i1, i2] = log P(X[i, i2] | G = G[i1], Y = y)
      logPX <- apply(logP, 1, sum) # logPX[j] = log P(X| G = G[j])
      probs[i, y] <- matrixStats::logSumExp(logPX) + log(w[y]) # should actually - log(nsample)
    }
    
    Y[i] <- which.max(probs[i,]) # get the posterior estimate 
  }
  
  list(posterior = probs, labels = Y) # return posterior distributions and the labels
}

#' Estimate the memberships.
#' 
#' @export
SPM_membership_Bernoulli <- function(X, Y, a, rho, Ts, P = NULL, logits = NULL, VI = FALSE, ntrace = 1000, nchain = 2, nskip = 2, seed = 1, save_trace = FALSE){
  if (is.null(logits)){ # if logits are not provided, calculate it from P
    logits <- log(P / (1 - P))
  } 
  SPM_Bernoulli_membership_stancode <- "
  data {
    int<lower=0> dg;  // dim(membership)
    int<lower=0> N;  // size of training set
    int<lower=0> nlabel;
    int<lower=0> ntopic;
    int<lower=0> d; //dim of data
    
    int<lower=0, upper=1> X[N, d];
    int<lower=1> Y[N];         // label, start from 1 
    
    matrix[ntopic, dg] T[nlabel];
    
    real<lower=0> a;
    vector[dg] rho;
    
    matrix[ntopic, d] logits;
  }
  parameters {
    simplex[dg] G[N];
  }
  model{
    for (i in 1:N){
      vector[ntopic] u;
      row_vector[d] logit_X;
      
      G[i] ~ dirichlet(a * rho);
      
      u = T[Y[i]] * G[i];
      logit_X = u' * logits;
      
      X[i] ~ bernoulli_logit(logit_X);
    }
  }
  "
  
  dg <- dim(Ts)[3]
  N <- length(Y)
  nlabel <- dim(Ts)[1]
  ntopic <- dim(Ts)[2]
  d <- dim(X)[2]
  K <- ntopic
  
  if(is.matrix(rho)){
    rho <- as.vector(rho)
  }
  
  dat_fit <- list(
    dg = dg, N = N, nlabel = nlabel, ntopic = ntopic, d = d,
    X = X, Y = Y,
    T = Ts,
    a = a, rho = rho,
    logits = logits
  )
  
  if (VI){
    model <- rstan::stan_model(model_code = SPM_Bernoulli_membership_stancode)
    fit_estimate <-rstan::vb(model, data = dat_fit, seed = seed)
  }else{
    fit_estimate <-rstan::stan(
      model_code = SPM_Bernoulli_membership_stancode,
      data = dat_fit,
      chains = nchain,
      iter = ntrace,
      seed = seed)
  }
  
  trace <- as.matrix(fit_estimate)
  
  res <- list()
  
  if(save_trace){
    res[["trace"]] <- trace
  }
  
  # save results
  nsample <- dim(trace)[1]
  nsave <- nsample %/% nskip
  index_save <- (1:nsave) * nskip
  trace <- trace[index_save,]  
  
  G <- matrix(NA, nrow = N, ncol = dg)
  for (i1 in 1:N){
    for(i2 in 1:dg){
      varname <- sprintf("G[%s,%s]", i1, i2)
      G[i1, i2] <- mean(trace[, varname])
    }
  } 
  
  res[["G"]] <- G
  
  res
}


