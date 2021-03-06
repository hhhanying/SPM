% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normal.R
\name{BPM_training_Normal}
\alias{BPM_training_Normal}
\title{Train a BPM model for normal distribution.}
\usage{
BPM_training_Normal(
  X,
  b,
  alpha,
  mu_Mu,
  sigma2_Mu,
  alpha_Lambda,
  beta_Lambda,
  ntopic,
  VI = FALSE,
  ntrace = 1000,
  nchain = 2,
  nskip = 2,
  seed = 1
)
}
\value{
A list containing all estimated parameters.
}
\description{
Train a BPM model for normal distribution.
}
