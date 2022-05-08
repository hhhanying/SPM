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

SPM_training_Normal <- function(X, Y, Ts, b, alpha, mu_Mu, sigma2_Mu, alpha_Lambda, beta_Lambda, ntrace, nchain, nskip){
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
  fit_train <-rstan::stan(model_code = SPM_Normal_stancode,
                          data = dat_fit,
                          chains = nchain,
                          iter = ntrace)
  
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

SPM_membership_Normal <- function(X, Y, Lambda, Mu, a, rho, Ts, ntrace, nchain, nskip, seed, w = NULL){
  SPM_Normal_membership_stancode <-"
  data {
      int<lower=0> dg;  // dim(membership)
      int<lower=0> N;  // size of training set
      int<lower=0> nlabel;
      int<lower=0> ntopic;
      int<lower=0> d; //dim of data
      
      real<lower=0> a;
      simplex[dg] rho;
  
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
  fit_estimate <-rstan::stan(model_code = SPM_Normal_membership_stancode,
                             data = dat_fit,
                             chains = nchain,
                             iter = ntrace)
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

BPM_training_Normal <- function(X, b, alpha, mu_Mu, sigma2_Mu, alpha_Lambda, beta_Lambda, ntopic, ntrace, nchain, nskip){
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
  N <- length(Y)
  d <- dim(X)[2]
  
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
  fit_train <-rstan::stan(model_code = BPM_Normal_stancode,
                          data = dat_fit,
                          chains = nchain,
                          iter = ntrace)
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

get_parameters_Normal <- function(X, Y, G, Ts, Lambda, Mu){
  Tau <- Mu * Lambda
  U <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N){
    U[i,] <- Ts[Y[i], ,] %*% G[i,]
  }
  
  LambdaX <- U %*% Lambda
  TauX <- U %*% Tau
  SigmaX <- sqrt(1 / LambdaX) 
  MuX <- TauX / LambdaX  
  
  list(MuX = MuX, LambdaX = LambdaX, SigmaX = SigmaX, MuX = MuX)
}