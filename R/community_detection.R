library(MASS)
library(stats)

## 得到模拟数据的真实参数real_parameter
simulation_alpha <- function(n, K) {
  real_alpha <- rep(1 / K, K)
  return(real_alpha)
}

simulation_mu <- function(n, K, dm) {
  real_mu <- matrix(0, nrow = K, ncol = K)
  for (i in seq(1, K)) {
    real_mu[i, i] <- runif(1, -dm, dm)
    for (j in seq(i+1, K)) {
      if (j < (K+1)){
        real_mu[i, j] <- runif(1, -dm, dm)
        real_mu[j, i] <- real_mu[i, j]
      }
    }
  }
  return(real_mu)
}

simulation_sigma <- function(n, K, dv) {
  real_sigma <- matrix(dv, nrow = K, ncol = K)
  return(real_sigma)
}

simulation_c <- function(n, K, alpha) {
  c <- sample(1:K, size = n, replace = TRUE, prob = alpha)
  return(c)
}

simulation_A <- function(n, mu, sigma, c) {
  A <- matrix(0, nrow = n, ncol = n)
  for (i in seq(1, n-1)) {
    for (j in seq(i+1, n)) {
      ith <- c[i]
      jth <- c[j]
      loc <- mu[ith, jth]
      scale <- sqrt(sigma[ith, jth])
      A[i, j] <- rnorm(1, mean = loc, sd = scale) + rnorm(1, mean = 0, sd = 0.1)
      A[j, i] <- A[i, j]
    }
  }
  return(A)
}

## 计算真实的成员矩阵
compute_Z <- function(n, K, c) {
  Z <- matrix(0, nrow = n, ncol = K)
  for (i in 1:n) {
    Z[i, c[i]] <- 1
  }
  return(Z)
}

## 得到算法迭代前的初始参数initial_parameter
initial_alpha <- function(K) {
  alpha <- rep(1 / K, K)
  return(alpha)
}

initial_mu <- function(K) {
  mu <- matrix(0, nrow = K, ncol = K)
  return(mu)
}

initial_sigma <- function(K) {
  sigma <- diag(K)
  return(sigma)
}

initial_tau <- function(n, K) {
  tau <- matrix(0, nrow = n, ncol = K)
  for (i in 1:n) {
    index <- sample(1:K, size = 1)
    tau[i, index] <- 1
    tau[i, ] <- tau[i, ] + 0.2
  }
  initial_tau <- t(apply(tau, 1, function(x) x / sum(x)))
  return(initial_tau)
}

## 每次迭代更新参数
update_alpha <- function(K, tau) {
  alpha <- colMeans(tau)
  return(alpha)
}

update_mu <- function(n, K, tau, A) {
  mu <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      sum_mu1 <- tau[, k] %*% t(tau[, l]) * A
      sum_mu2 <- tau[, k] %*% t(tau[, l])
      mu[k, l] <- (sum(sum_mu1) - sum(diag(sum_mu1))) / (sum(sum_mu2) - sum(diag(sum_mu2)))
    }
  }
  return(mu)
}

update_sigma <- function(n, K, tau, A, mu) {
  sigma <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      diff <- A - matrix(mu[k, l], nrow = n, ncol = n)
      sum_sigma1 <- tau[, k] %*% t(tau[, l]) * diff^2
      sum_sigma2 <- tau[, k] %*% t(tau[, l])
      sigma[k, l] <- (sum(sum_sigma1) - sum(diag(sum_sigma1))) / (sum(sum_sigma2) - sum(diag(sum_sigma2)))
    }
  }
  return(sigma)
}

update_tau <- function(n, K, tau, A, alpha, mu, sigma) {
  for (i in 1:n) {
    ratio <- numeric(K)
    mini_ratio <- numeric(K)
    for (k in 1:K) {
      temp <- 2 * pi * sigma[k, ]^(-0.5)
      temp1 <- (A[i, ] - mu[k, ])^2 / (-2 * sigma[k, ])
      temp <- temp * exp(temp1)
      temp <- temp^tau
      mini_temp <- temp^(1 / 100)
      mini_alpha <- alpha^(1 / 100)
      mini_ratio[k] <- mini_alpha[k] * prod(mini_temp)
    }
    mini_sum <- sum(mini_ratio)
    ratio <- (mini_ratio / mini_sum)^100
    tau[i, ] <- ratio / sum(ratio)
  }
  return(tau)
}

## 得到估计社区
get_es_c <- function(n, tau) {
  es_c <- max.col(tau)
  return(es_c)
}

#' @title Normal SBM Matrix Generator using R
#' @description Simulating Normal Weighted Stochastic Block Model Matrices using R
#' @param dm the parameter that determines the normal distribution mean for generating the adjacency matrix
#' @param dv the parameter that determines the normal distribution variance for generating the adjacency matrix
#' @return A list containing:
#' \item{alpha}{The probability of community assignments for each node.}
#' \item{mu}{The true mean of the normal distribution used to generate the matrix.}
#' \item{sigma}{The true variance of the normal distribution used to generate the matrix.}
#' \item{c}{A vector of community assignments for the nodes.}
#' \item{A}{A generated normal weighted stochastic block matrix.}
#' @export
real_ac <- function(dm, dv) {
  set.seed(221)
  n <- 1000
  K <- 4
  real_alpha <- simulation_alpha(n, K)
  real_mu <- simulation_mu(n, K, dm)
  real_sigma <- simulation_sigma(n, K, dv)
  c <- simulation_c(n, K, real_alpha)
  A <- simulation_A(n, real_mu, real_sigma, c)
  list(
    alpha = real_alpha,
    mu = real_mu,
    sigma = real_sigma,
    c = c,
    A = A
  )
}

#' @title Community Detection using the Improved EM Algorithm in R
#' @description This function repeatedly estimates the community assignments by applying the Improved EM algorithm to the matrix parameters
#' @param real the input matrix of data
#' @return a vector of the estimated community assignments
#' @export

main <- function(real) {
  n <- 1000
  K <- 4
  c <- real$c
  A <- real$A
  alpha <- initial_alpha(K)
  mu <- initial_mu(K)
  sigma <- initial_sigma(K)
  tau <- initial_tau(n, K)
  for (w in 1:5) {
    alpha <- update_alpha(K, tau)
    mu <- update_mu(n, K, tau, A)
    sigma <- update_sigma(n, K, tau, A, mu)
    tau <- update_tau(n, K, tau, A, alpha, mu, sigma)
  }
  es_c <- get_es_c(n, tau)
  return(es_c)
}


#' @import MASS
#' @import Rcpp
#' @import stats
#' @importFrom boot boot
#' @import bootstrap
#' @import DAAG
#' @import coda
#' @import ggplot2
#' @import pracma
#' @import lpSolve
#' @import microbenchmark

#' @useDynLib SA24204154
NULL

