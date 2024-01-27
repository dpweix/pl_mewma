library("tidyverse")
library("spcov")
library("CovTools")
library("mgcov")
library("parameters")
library("here")
source(here("scripts", "pl_fd_helper.R"))

# Penalized Likelihood Fault Detection Methods
# Simultaneous Mean and Covariance

# Hawkins, Maboudou-Tchao (2008) ------------------------------------------

hawkins_2008 <- function(dat, mu_0, sigma_0, beta = .2) {
  A <- get_A(sigma_0)
  U <- calc_U(dat, mu_0, A)
  S <- calc_S(U, sigma_0, beta)
  
  pstat_hawkins(S)
}

# Wang 2014, Method 1 -----------------------------------------------------
# Uses an unmodified spcov
wang_2014_1 <- function(dat, mu_0, sigma_0,
                        alpha = .2, beta = .2,
                        lambda_1 = .2, lambda_2 = .2) {
  X <- as.matrix(dat)
  sig_inv <- solve(sigma_0)
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  mu_lt <- mean_SPARSE(mu_t, mu_0, sig_inv, pl_mu_1, lambda_1 = lambda_1)
  sigma_lt     <- cov_MEWMC(X, mu_lt, beta = beta)
  sigma_hat_lt <- cov_SPARSE(sigma_lt, sigma_0, lambda = lambda_2)
  sigma_0t     <- cov_MEWMC_0(X, mu_0, sigma_0, beta = beta)
  
  pmap_dfr(list(sigma_hat_lt, sigma_lt, sigma_0t, list(sig_inv)), pstat_wang_1)
}

# Wang 2014, Method 2 -----------------------------------------------------
# Uses an unmodified spcov
wang_2014_2 <- function(dat, mu_0, sigma_0,
                        alpha = .2, beta = .2,
                        lambda_1 = .2, lambda_2 = .2) {
  X <- as.matrix(dat)
  sig_inv <- solve(sigma_0)
  
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  mu_lt <- mean_SPARSE(mu_t, mu_0, sig_inv, pl_mu_1, lambda_1 = lambda_1)
  sigma_0t     <- cov_MEWMC_0(X, mu_0, sigma_0, beta = beta)
  sigma_hat_0t <- cov_SPARSE(sigma_0t, sigma_0, lambda_2 = lambda_2)
  
  pmap_dfr(list(pmap(mu_lt, c), list(mu_0), sigma_0t, sigma_hat_0t, list(sig_inv)), 
           pstat_wang_2)
}


# Test Method 1 -----------------------------------------------------------

# Hawkins + LASSO based Sparse Covariance Estimation
MC_LASSO <- function(dat, mu_0, sigma_0,
                          beta = .2, lambda_s = .2) {
  A <- get_A(sigma_0)
  U <- calc_U(dat, mu_0, A)
  S <- calc_S(U, sigma_0, beta)
  S <- sparsify_S(S, lambda_s)
  pstat_hawkins(S)
}

# Hawkins + COMET based Sparse Covariance Estimation
MC_COMET <- function(dat, mu_0, sigma_0,
                          beta = .2, cutoff = .1, n_w = 50) {
  A <- get_A(sigma_0)
  U <- calc_U(dat, mu_0, A)
  S <- calc_S(U, sigma_0, beta)
  S <- get_S_hat(threshold_S(S, cutoff = cutoff),
                 get_U_window(U, n_w = n_w),
                 U)
  
  pstat_hawkins(S)
}

# Wang + COMET based sparse covariance estimation
MAC_COMET <- function(dat, mu_0, sigma_0,
                          alpha = .2, beta = .2, cutoff = .1, n_w = 50) {
  A <- get_A(sigma_0)
  X <- as.matrix(dat)
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  U <- calc_U(dat, mu_t, A)
  S <- calc_S(U, sigma_0, beta)
  S_hat <- get_S_hat(threshold_S(S, cutoff = cutoff),
                     get_U_window(U, n_w = n_w),
                     U)
  
  pstat_test_2(mu_t, mu_0, S_hat, solve(sigma_0))
}

# Wang + COMET based sparse covariance estimation
MAC_COMET1 <- function(dat, mu_0, sigma_0,
                      alpha = .2, beta = .2, cutoff = .1, n_w = 50) {
  A <- get_A(sigma_0)
  X <- as.matrix(dat)
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  U <- calc_U(dat, mu_t, A)
  S <- calc_S(U, sigma_0, beta)
  S_hat <- get_S_hat(threshold_S(S, cutoff = cutoff),
                     get_U_window(U, n_w = n_w),
                     U)
  
  pstat <- pstat_test_3(lapply(seq_len(nrow(mu_t)), function(i) mu_t[i, ]),
               mu_0,
               S_hat,
               solve(sigma_0))
}
# test_method_4 <- function(dat, mu_0, sigma_0, window_upper_limit,
#                           alpha = .2, beta = .2, cutoff = .1,
#                           n_min = 50, n_max = 100) {
#   A <- get_A(sigma_0)
#   X <- as.matrix(dat)
#   mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
#   U <- calc_U(dat, mu_t, A)
#   S <- calc_S(U, sigma_0, beta)
#   
#   S_star <- threshold_S(S, cutoff = cutoff)
#   U_window <- create_adaptive_window(dat, n_min, n_max,
#                                      mu_0 = mu_0, sigma_0 = sigma_0,
#                                      upper_limit = window_upper_limit,
#                                      condition = "mean")
#   
#   S_hat <- get_S_hat(S_star, U_window, U)
#   pstat_test_2(mu_t, mu_0, S_hat, solve(sigma_0))
# }
# 
# test_method_5 <- function(dat, mu_0, sigma_0, window_upper_limit,
#                           alpha = .2, beta = .2, cutoff = .1,
#                           n_min = 50, n_max = 100) {
#   p <- ncol(dat)
#   A <- get_A(sigma_0)
#   X <- as.matrix(dat)
#   mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
#   U <- calc_U(dat, mu_t, A)
#   S <- calc_S(U, sigma_0, beta)
#   sigma_S_0 <- cov(matrix(map(S, ks::vech) |> unlist(), ncol = p*(p-1), byrow = TRUE))
#   
#   S_star <- threshold_S(S, cutoff = cutoff)
#   U_window <- create_adaptive_window(dat, n_min, n_max,
#                                      mu_0 = mu_0, sigma_0 = sigma_0,
#                                      S = S, sigma_S_0 = sigma_S_0,
#                                      upper_limit = window_upper_limit,
#                                      condition = "variance")
#   
#   S_hat <- get_S_hat(S_star, U_window, U)
#   pstat_test_2(mu_t, mu_0, S_hat, solve(sigma_0))
# }



# Apply methods in single function ----------------------------------------

# dat: a matrix or data frame with the monitored variables
# n: the size of the training set (the first n observations are assumed IC)

apply_pl_fd <- function(df, n = nrow(df)/2, method = "MAC_COMET") {
  
  if(method == "hawkins") {
    mu_0    <- colMeans(df[1:n, ])
    sigma_0 <- cov(df[1:n, ])
    
    pstat <- hawkins_2008(df, mu_0, sigma_0, selected_method$beta)$pstat
    
  } else if(method == "wang_1") {
    mu_0    <- colMeans(df[1:n, ])
    sigma_0 <- cov(df[1:n, ])
    
    pstat <- wang_2014_1(df, mu_0, sigma_0, 
                         alpha = selected_method$alpha,
                         beta = selected_method$beta,
                         lambda_1 = selected_method$lambda_1,
                         lambda_2 = selected_method$lambda_2)$pstat
    
  } else if(method == "MC_LASSO") {
    mu_0    <- colMeans(df[1:n, ])
    sigma_0 <- spcov(cov(df[1:n, ]), cov(df[1:n, ]), lambda = selected_method$lambda)
    
    pstat <- MC_LASSO(df, mu_0, sigma_0,
                      lambda = selected_method$lambda,
                      lambda_s = selected_method$lambda_s)$pstat
    
  } else if(method == "MC_COMET") {
    mu_0    <- colMeans(df[1:n, ])
    S_0 <- cov(df[1:n, ])
    S_0[(abs(S_0) < selected_method$cutoff)] <- 0
    sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
    
    pstat <- MC_COMET(df, mu_0, sigma_0,
                      lambda = selected_method$lambda,
                      cutoff = selected_method$cutoff,
                      n_w = selected_method$n_w)$pstat
    
  } else if(method == "MAC_COMET") {
    mu_0    <- colMeans(df[1:n, ])
    S_0 <- cov(df[1:n, ])
    S_0[(abs(S_0) < selected_method$cutoff)] <- 0
    sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
    
    pstat <- MAC_COMET(df, mu_0, sigma_0,
                       alpha = selected_method$alpha,
                       beta = selected_method$beta,
                       cutoff = selected_method$cutoff,
                       n_w = selected_method$n_w)
    T_1 <- pstat$T_1
    T_2 <- pstat$T_2
  } 
  
  if(method %in% c("hawkins", "wang_1", "MC_LASSO", "MC_COMET")) {
    h  <- estimate_h(pstat[1:n], arl_ic)
    rl <- get_rl_oc(pstat[(n+1):(2*n)], h = h)
    
    list(method = method,
         pstat = pstat,
         params_data = list(n = n, arl_ic = arl_ic),
         params_method = selected_method,
         h = h,
         rl = rl)
    
  } else if(method %in% c("wang_2", "MAC_COMET")) {
    h1  <- estimate_h(T_1[1:n], arl_ic)
    h2  <- estimate_h(T_2[1:n], arl_ic)
    
    rl1 <- get_rl_oc(T_1[(n+1):(2*n)], h = h1)
    rl2 <- get_rl_oc(T_2[(n+1):(2*n)], h = h2)
    
    rl <- min(rl1, rl2)
    
    list(method = method,
         pstat = tibble(T_1 = T_2, T_2 = T_2),
         params_data = list(n = n, arl_ic = arl_ic),
         params_method = selected_method,
         h1 = h1,
         h2 = h2,
         rl = rl)
  }
}
