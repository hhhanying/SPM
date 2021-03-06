#' Simulate data from given parameters.
#'
#' @param a: dispersion parameter for membership
#' @param rho: weight parameter for membership
#' @return A list res containing: 1) X: the document; 2) Y: the labels; 3) G: the untransformed memberships and 4) U: the transformed memberships. 
#' @export
document_generator_Normal <- function(a, rho, Ts, Lambda, Tau, N, w = NULL, seed = NULL){
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  nlabel <- dim(Ts)[1]
  d <- dim(Tau)[2]
  K <- dim(Ts)[2]
  
  Y <- sample(x = nlabel, size = N, replace = TRUE, prob = w)
  G <- gtools::rdirichlet(n = N, alpha = a * rho)
  #browser()
  U <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N){
    U[i,] <- Ts[Y[i], ,] %*% G[i,]
  }
  
  LambdaX <- U %*% Lambda
  TauX <- U %*% Tau
  SigmaX <- sqrt(1 / LambdaX) 
  MuX <- TauX / LambdaX
  #browser()
  X <- matrix(NA, nrow = N, ncol = d)
  for(i in 1:N){
    X[i,] <- rnorm(n = d, mean = MuX[i,], sd = SigmaX[i,])
  }
  
  list(X = X, Y = Y, G = G, U = U)  
}

#' Given the estimated topics and memberships, simulate data
#' @return X: the simulated data
#' @export
SPM_simulator_Normal <- function(Lambda = NULL, Tau = NULL, Mu = NULL, S = NULL, Ts = NULL, G = NULL, Y = NULL, U = NULL, seed = NULL){
  # if random seed is provided
  if (!is.null(seed)){
    set.seed(seed)
  }
  # if membership is untransformed
  if (is.null(U)){
    K <- dim(Ts)[2]
    N <- length(Y)
    
    U <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N){
      U[i,] <- Ts[Y[i], ,] %*% G[i,]
    }
    
  } else{
    N <- dim(U)[1]
  }
  # what form of parameters are given
  if (is.null(Lambda)){
    Lambda <- 1 / S
  }
  if (is.null(Tau)){
    Tau <- Mu * Lambda
  }
  
  d <- dim(Tau)[2]
  
  LambdaX <- U %*% Lambda
  TauX <- U %*% Tau
  SigmaX <- sqrt(1 / LambdaX) 
  MuX <- TauX / LambdaX
  
  X <- matrix(NA, nrow = N, ncol = d)
  for(i in 1:N){
    X[i,] <- rnorm(n = d, mean = MuX[i,], sd = SigmaX[i,])
  }
  
  X 
}
#' Fit a SPM model on the training set given the hyperparameters.
#' @return A list containing all estimated parameters.
#' @export
SPM_training_Normal <- function(X, Y, Ts, b, alpha, mu_Mu, sigma2_Mu, alpha_Lambda, beta_Lambda, VI = FALSE, ntrace = 1000, nchain = 2, nskip = 2, seed = 1){
  SPM_Normal_stancode <-"
  data {
    int<lower=0> dg;  // dim(membership)
    int<lower=0> N;  // size of training set
    int<lower=0> nlabel;
    int<lower=0> ntopic;
    int<lower=0> d; //dim of data
    
    matrix[N, d] X;
    int<lower=1> Y[N];         // label, start from 1 
    
    matrix[ntopic, dg] T[nlabel];
    
    real<lower=0> b;
    vector<lower=0>[dg] alpha;
    matrix[ntopic, d] mu_Mu;
    matrix[ntopic, d] sigma2_Mu;
    matrix[ntopic, d] alpha_Lambda;
    matrix[ntopic, d] beta_Lambda;
    
  }
  parameters {
    real<lower=0> a;
    simplex[dg] rho;
    
    simplex[dg] G[N];
    
    matrix[ntopic, d] Mu;
    matrix<lower=0>[ntopic, d] Lambda;

  }
  model{
    matrix[ntopic, d] Tau;
    vector[ntopic] U[N];
    row_vector[d] lambda_X[N];
    row_vector[d] tau_X[N];
    
    a ~ exponential(b);
    rho ~ dirichlet(alpha);

    // topics
    for (i in 1:ntopic){
      Lambda[i] ~ gamma(alpha_Lambda[i], beta_Lambda[i]);
      Mu[i] ~ normal(mu_Mu[i], sqrt(sigma2_Mu[i] ./ Lambda[i]));
      Tau[i] = Mu[i] .* Lambda[i];
    }
    
    for(i in 1:N){
      G[i] ~ dirichlet(a * rho);
      U[i] = T[Y[i]] * G[i];
      lambda_X[i] = U[i]' * Lambda;
      tau_X[i] = U[i]' * Tau;
      X[i] ~ normal(tau_X[i] ./ lambda_X[i], sqrt(1 ./ lambda_X[i]));
    }

  }
  "

  # calculate needed parameters (avoid too many inputs)
  N <- length(Y)
  nlabel <- dim(Ts)[1]
  dg <- dim(Ts)[3]
  ntopic <- dim(Ts)[2]
  d <- dim(X)[2]
  K <- ntopic
  
  # we assume all hyperparameters for topics are matrices, if they are not, we will duplicate it to make a matrix
  # mu_Mu = SPM::trans_to_matrix(mu_Mu, ntopic, d)
  # sigma2_Mu = SPM::trans_to_matrix(sigma2_Mu, ntopic, d)
  # alpha_Lambda = SPM::trans_to_matrix(alpha_Lambda, ntopic, d)
  # beta_Lambda = SPM::trans_to_matrix(beta_Lambda, ntopic, d)
  mu_Mu = trans_to_matrix(mu_Mu, ntopic, d)
  sigma2_Mu = trans_to_matrix(sigma2_Mu, ntopic, d)
  alpha_Lambda = trans_to_matrix(alpha_Lambda, ntopic, d)
  beta_Lambda = trans_to_matrix(beta_Lambda, ntopic, d)
  
  # the data to fit
  dat_fit <- list(
    dg = dg,
    N = N,
    nlabel = nlabel,
    ntopic = ntopic,
    d = d,
    
    X = X,
    Y = Y,
    
    T = Ts,
    
    b = b, 
    alpha = alpha,
    mu_Mu = mu_Mu,
    sigma2_Mu = sigma2_Mu,
    alpha_Lambda = alpha_Lambda,
    beta_Lambda = beta_Lambda
  )
  
  # sampling
  if(VI){
    model <- rstan::stan_model(model_code = SPM_Normal_stancode)
    fit_train <-rstan::vb(model, data = dat_fit, seed = seed)
  }else{
    fit_train <-rstan::stan(model_code = SPM_Normal_stancode,
                            data = dat_fit,
                            chains = nchain,
                            iter = ntrace,
                            seed = seed)
  }
  
  trace <- as.matrix(fit_train)
  
  # save results
  nsample <- dim(trace)[1]
  nsave <- nsample %/% nskip
  index_save <- (1:nsave) * nskip
  trace <- trace[index_save,]
  
  
  res <- list()
  
  paras <- c("G", "Lambda", "Mu")
  d1s <- c(N, K, K)
  d2s <- c(dg, d, d)
  for (i in 1:3){
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

#' Given the SPM model, predict labels for a new document.
#' @return A list containing the posterior distribution of labels and the estimated labels.
#' @export
SPM_predicting_Normal <- function(X, Lambda, Mu, a, rho, Ts, nsample, seed, w = NULL){
  nlabel <- dim(Ts)[1]
  N <- dim(X)[1]
  Tau <- Mu * Lambda
  
  probs <- matrix(NA, nrow = N, ncol = nlabel) # stores the posterior probability of labels
  Y <- rep(0, N)
  
  if (is.null(w)) w <- rep(1, nlabel) # if no information about the prior of labels, treat it as uniform
  
  set.seed(seed)
  for(i in 1:N){
    G <- gtools::rdirichlet(n = nsample, alpha = a * rho) # draw G
    
    for (y in 1:nlabel){
      U <- G %*% t(Ts[y,,]) # nsample * K
      
      LambdaX <- U %*% Lambda
      TauX <- U %*% Tau
      MuX <- TauX / LambdaX
      SigmaX <- sqrt(1 / LambdaX)
      
      x <- matrix(rep(X[i,], nsample), nrow = nsample, byrow = TRUE)
      p_matrix <- dnorm(x = x, mean = MuX, sd = sqrt(1 / LambdaX), log = TRUE) # get the nsample * d logP matrix
      logps <- apply(p_matrix, 1, sum) # row sum: logP for each G
      probs[i, y] <- matrixStats::logSumExp(logps) + log(w[y]) # should actually - log(nsample)
      
    }
    
    Y[i] <- which.max(probs[i,])
  }
  
  list(posterior = probs, labels = Y)
}

#' Given the SPM model, estimate the membership for new data points.
#' @return A matrix of the estimated memberships.
#' @export
SPM_membership_Normal <- function(X, Y, Lambda, Mu, a, rho, Ts, VI = FALSE, ntrace = 1000, nchain = 2, nskip = 2, seed = 1){
  SPM_Normal_membership_stancode <-"
  data {
      int<lower=0> dg;  // dim(membership)
      int<lower=0> N;  // size of training set
      int<lower=0> nlabel;
      int<lower=0> ntopic;
      int<lower=0> d; //dim of data
      
      real<lower=0> a;
      vector[dg] rho;
  
      matrix[N, d] X;
      int<lower=1> Y[N];         // label, start from 1 
  
      matrix[ntopic, dg] T[nlabel];
      
  
      matrix[ntopic, d] Lambda;
      matrix[ntopic, d] Tau;
  
  }
  parameters {
      simplex[dg] G[N];
  }
  model{
      for (i in 1:N){
          vector[ntopic] u;
          row_vector[d] lambda_X;
          row_vector[d] tau_X;
          
          G[i] ~ dirichlet(a * rho);
          u = T[Y[i]] * G[i];
  
          lambda_X = u' * Lambda;
          tau_X = u' * Tau;
  
          for (j in 1:d){
              X[i, j] ~ normal(tau_X[j] / lambda_X[j], sqrt(1 / lambda_X[j]));
          }
      }
  }
  "
  
  # calculate needed parameters (avoid too many inputs)
  N <- length(Y)
  nlabel <- dim(Ts)[1]
  dg <- dim(Ts)[3]
  ntopic <- dim(Ts)[2]
  d <- dim(X)[2]
  
  Tau <- Mu * Lambda
  
  if(is.matrix(rho)){
    rho <- as.vector(rho)
  }
  
  dat_fit <- list(
    dg = dg,
    N = N,
    nlabel = nlabel,
    ntopic = ntopic,
    d = d,
    
    a = a,
    rho = rho,
    
    X = X,
    Y = Y,
    
    T = Ts,
    Lambda = Lambda,
    Tau = Tau
  )
  
  # sampling
  if(VI){
    model <- rstan::stan_model(model_code = SPM_Normal_membership_stancode)
    fit_estimate <-rstan::vb(model, data = dat_fit, seed = seed)
  }else{
    fit_estimate <-rstan::stan(model_code = SPM_Normal_membership_stancode,
                            data = dat_fit,
                            chains = nchain,
                            iter = ntrace,
                            seed = seed)
  }

  trace <- as.matrix(fit_estimate)
  
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
  
  G
}

#' Train a BPM model for normal distribution.
#' @return A list containing all estimated parameters.
#' @export
BPM_training_Normal <- function(X, b, alpha, mu_Mu, sigma2_Mu, alpha_Lambda, beta_Lambda, ntopic, VI = FALSE, ntrace = 1000, nchain = 2, nskip = 2, seed = 1){
  BPM_Normal_stancode <-"
  data {
    int<lower=0> N;  // size of training set
    int<lower=0> ntopic;
    int<lower=0> d; //dim of data
    
    matrix[N, d] X;
    
    real<lower=0> b;
    vector<lower=0>[ntopic] alpha;
    matrix[ntopic, d] mu_Mu;
    matrix[ntopic, d] sigma2_Mu;
    matrix[ntopic, d] alpha_Lambda;
    matrix[ntopic, d] beta_Lambda;
  }
  parameters {
    real<lower=0> a;
    simplex[ntopic] rho;
    
    simplex[ntopic] U[N];
    
    matrix[ntopic, d] Mu;
    matrix<lower=0>[ntopic, d] Lambda;
  }
  model{
    matrix[ntopic, d] Tau;
    
    a ~ exponential(b);
    rho ~ dirichlet(alpha);
    
    // topics
    for (i in 1:ntopic){
      for(j in 1:d){
        Lambda[i, j] ~ gamma(alpha_Lambda[i, j], beta_Lambda[i, j]);
        Mu[i, j] ~ normal(mu_Mu[i, j], sqrt(sigma2_Mu[i, j] / Lambda[i, j]));
        Tau[i, j] = Mu[i, j] * Lambda[i, j];
      }
    }
    
    for (i in 1:N){
      row_vector[d] lambda_X;
      row_vector[d] tau_X;
      
      U[i] ~ dirichlet(a * rho);
      
      lambda_X = U[i]' * Lambda;
      tau_X = U[i]' * Tau;
      
      for (j in 1:d){
        X[i, j] ~ normal(tau_X[j] / lambda_X[j], sqrt(1 / lambda_X[j]));
      }
    }
  }
  "
  
  # calculate needed parameters (avoid too many inputs)
  N <- dim(X)[1]
  d <- dim(X)[2]
  K <- ntopic
  
  # we assume all hyperparameters for topics are matrices, if they are not, we will duplicate it to make a matrix
  mu_Mu = trans_to_matrix(mu_Mu, ntopic, d)
  sigma2_Mu = trans_to_matrix(sigma2_Mu, ntopic, d)
  alpha_Lambda = trans_to_matrix(alpha_Lambda, ntopic, d)
  beta_Lambda = trans_to_matrix(beta_Lambda, ntopic, d)
  
  # the data to fit
  dat_fit <- list(
    N = N,
    ntopic = ntopic,
    d = d,
    
    X = X,
    
    b = b, 
    alpha = alpha,
    mu_Mu = mu_Mu,
    sigma2_Mu = sigma2_Mu,
    alpha_Lambda = alpha_Lambda,
    beta_Lambda = beta_Lambda
  )
  
  # sampling
  if(VI){
    model <- rstan::stan_model(model_code = BPM_Normal_stancode)
    fit_train <-rstan::vb(model, data = dat_fit, seed = seed)
  }else{
    fit_train <-rstan::stan(model_code = BPM_Normal_stancode,
                            data = dat_fit,
                            chains = nchain,
                            iter = ntrace,
                            seed = seed)
  }

  trace <- as.matrix(fit_train)
  
  # save results
  nsample <- dim(trace)[1]
  nsave <- nsample %/% nskip
  index_save <- (1:nsave) * nskip
  trace <- trace[index_save,]
  
  res <- list()
  
  paras <- c("U", "Lambda", "Mu")
  d1s <- c(N, K, K)
  d2s <- c(ntopic, d, d)
  for (i in 1:3){
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
  
  tem <- rep(0, ntopic)
  for (i in 1:ntopic){
    varname <- sprintf("rho[%s]", i)
    tem[i] <- mean(trace[, varname])
  }
  res[["rho"]] <- tem
  
  res
}

#' Given the BPM model, estimate the membership for new data points.
#' @return A matrix of the estimated memberships.
#' @export
BPM_membership_Normal <- function(X, Lambda, Mu, a, rho, VI = FALSE, ntrace = 1000, nchain = 2, nskip = 2, seed = 1){
  BPM_Normal_membership_stancode <-"
  data {
      int<lower=0> N;  // size of training set
      int<lower=0> ntopic;
      int<lower=0> d; //dim of data
      
      real<lower=0> a;
      vector[ntopic] rho;
  
      matrix[N, d] X;

      matrix[ntopic, d] Lambda;
      matrix[ntopic, d] Tau;
  
  }
  parameters {
      simplex[ntopic] U[N];
  }
  model{
      for (i in 1:N){
          row_vector[d] lambda_X;
          row_vector[d] tau_X;
          
          U[i] ~ dirichlet(a * rho);
  
          lambda_X = U[i]' * Lambda;
          tau_X = U[i]' * Tau;
  
          for (j in 1:d){
              X[i, j] ~ normal(tau_X[j] / lambda_X[j], sqrt(1 / lambda_X[j]));
          }
      }
  }
  "
  
  # calculate needed parameters (avoid too many inputs)
  N <- dim(X)[1]
  d <- dim(X)[2]
  ntopic <- dim(Lambda)[1]
  Tau <- Mu * Lambda
  K <- ntopic
  
  dat_fit <- list(
    N = N,
    ntopic = ntopic,
    d = d,
    a = a,
    rho = rho,
    X = X,
    Lambda = Lambda,
    Tau = Tau
  )
  
  # sampling
  if(VI){
    model <- rstan::stan_model(model_code = BPM_Normal_membership_stancode)
    fit_estimate <-rstan::vb(model, data = dat_fit, seed = seed)
  }else{
    fit_estimate <-rstan::stan(model_code = BPM_Normal_membership_stancode,
                            data = dat_fit,
                            chains = nchain,
                            iter = ntrace,
                            seed = seed)
  }

  trace <- as.matrix(fit_estimate)
  
  # save results
  nsample <- dim(trace)[1]
  nsave <- nsample %/% nskip
  index_save <- (1:nsave) * nskip
  trace <- trace[index_save,]
  
  U <- matrix(NA, nrow = N, ncol = ntopic)
  for (i1 in 1:N){
    for(i2 in 1:ntopic){
      varname <- sprintf("U[%s,%s]", i1, i2)
      U[i1, i2] <- mean(trace[, varname])
    }
  } 
  
  U
}

#' Calculate the generated distribution given the membership and the topics.
#' @return A list of the parameters of the distributions of data given its memberships.
#' @export
get_parameters_Normal <- function(Lambda, Mu, G = NULL, Ts = NULL, Y = NULL, U = NULL){
  Tau <- Mu * Lambda
  K <- dim(Lambda)[1]
  if (is.null(U)){
    N <- length(Y)
    U <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N){
      U[i,] <- Ts[Y[i], ,] %*% G[i,]
    }
  } 
  
  LambdaX <- U %*% Lambda
  TauX <- U %*% Tau
  SigmaX <- sqrt(1 / LambdaX) 
  MuX <- TauX / LambdaX  
  
  list(MuX = MuX, LambdaX = LambdaX, SigmaX = SigmaX, MuX = MuX)
}
