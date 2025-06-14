# General FD --------------------------------------------------------------
get_arl_ic <- function(vec_pstat, h) {
  denom <- sum(vec_pstat > h)
  numer <- length(vec_pstat)
  ifelse(denom == 0, numer, numer/denom)
}

get_rl_oc <- function(vec_pstat, h) {
  first_fault <- which(vec_pstat > h)[1]
  ifelse(is.na(first_fault), length(vec_pstat), first_fault)
}

get_rl_oc2 <- function(vec_pstat1, vec_pstat2, h1, h2) {
  first_fault1 <- which(vec_pstat1 > h1)[1]
  first_fault2 <- which(vec_pstat2 > h2)[1]
  first_fault <- min(first_fault1, first_fault2)
  ifelse(is.na(first_fault), length(vec_pstat1), first_fault)
}

get_detection_rate <- function(vec_pstat, h) {
  sum(vec_pstat > h)/length(vec_pstat)
}

get_detection_rate2 <- function(vec_pstat1, vec_pstat2, h1, h2) {
  sum(vec_pstat1 > h1 | vec_pstat2 > h2)/length(vec_pstat1)
}

objective_function <- function(h, arl, vec_pstat) {
  abs(get_arl_ic(vec_pstat, h) - arl)
  #abs(get_rl_oc(vec_pstat, h) - arl)
}

objective_function_fd <- function(h, alpha, vec_pstat) {
  abs(sum(vec_pstat > h)/length(vec_pstat) - alpha)
}

estimate_h <- function(vec_pstat, arl) {
  vec_pstat <- vec_pstat[which(!is.na(vec_pstat))]
  h <- as.numeric(quantile(vec_pstat, 1 - 1/arl)) #MRL
  optim(h, objective_function, arl = arl, vec_pstat = vec_pstat, 
        method = "BFGS")$par
}

estimate_h_fd <- function(vec_pstat, alpha) {
  #quantile(vec_pstat, 1-alpha, na.rm = TRUE) |> as.numeric()
  
  vec_pstat <- vec_pstat[which(!is.na(vec_pstat))]
  h <- as.numeric(quantile(vec_pstat, 1 - alpha)) # h for MRL
  optim(h, objective_function_fd, alpha = alpha, vec_pstat = vec_pstat, 
        method = "BFGS")$par
}

objective_function2 <- function(h, arl, vec_pstat1, vec_pstat2) {
  h1 <- h[1]
  h2 <- h[2]
  
  denom <- sum(vec_pstat1 > h1 | vec_pstat2 > h2)
  numer <- length(vec_pstat1)
  # list(estimate = ifelse(denom == 0, numer, numer/denom),
  #      theoretical = arl)
  abs(ifelse(denom == 0, numer, numer/denom) - arl)
}

objective_function2_fd <- function(h, alpha, vec_pstat1, vec_pstat2) {
  h1 <- h[1]
  h2 <- h[2]
  
  numer <- sum(vec_pstat1 > h1 | vec_pstat2 > h2)
  denom <- length(vec_pstat1)
  abs(numer/denom - alpha)
}

estimate_h2 <- function(vec_pstat1, vec_pstat2, arl) {
  vec_pstat1 <- vec_pstat1[which(!is.na(vec_pstat1))]
  vec_pstat2 <- vec_pstat2[which(!is.na(vec_pstat2))]
  h1 <- as.numeric(quantile(vec_pstat1, 1 - 1/arl)) #MRL
  h2 <- as.numeric(quantile(vec_pstat2, 1 - 1/arl)) #MRL
  h <- c(h1, h2)
  optim(h, objective_function2, arl = arl,
        vec_pstat1 = vec_pstat1, vec_pstat2 = vec_pstat2,
        method = "BFGS")$par
}

estimate_h2_fd <- function(vec_pstat1, vec_pstat2, alpha) {
  # h <- c(quantile(vec_pstat1, 1-alpha/2, na.rm = TRUE) |> as.numeric(),
  #        quantile(vec_pstat2, 1-alpha/2, na.rm = TRUE) |> as.numeric())
  
  vec_pstat1 <- vec_pstat1[which(!is.na(vec_pstat1))]
  vec_pstat2 <- vec_pstat2[which(!is.na(vec_pstat2))]
  h1 <- as.numeric(quantile(vec_pstat1, 1 - alpha/2)) #MRL
  h2 <- as.numeric(quantile(vec_pstat2, 1 - alpha/2)) #MRL
  h <- c(h1, h2)
  optim(h, objective_function2_fd, alpha = alpha,
        vec_pstat1 = vec_pstat1, vec_pstat2 = vec_pstat2,
        method = "BFGS")$par
}

estimate_h_bootstrap <- function(vec_pstat, arl, B = 500) {
  map_dbl(1:B, \(i) {
    estimate_h(vec_pstat = sample(vec_pstat, length(vec_pstat), replace = TRUE),
               arl = arl)
    }) |> mean() # test vs median
}

demean = function(dat){ # Unexported method from mgcov (but necessary)
  meanmat = matrix(rep(colMeans(dat), nrow(dat)), ncol = ncol(dat), byrow = TRUE)
  datdm = dat - meanmat
  return(datdm)
}

split_matrix <- function(X) {
  lapply(seq_len(nrow(X)), function(i) X[i, ])
}

# Hawkins, Maboudou-Tchao Methods -----------------------------------------

get_A <- function(sigma_0) {
  solve(t(chol(sigma_0)))
}

calc_U <- function(dat, mu_0, A) {
  X <- as.matrix(dat)
  if(length(mu_0) == ncol(X))    t(apply(X, 1, \(x) {A %*% (x - mu_0)}))
  else if(nrow(mu_0) == nrow(X)) t(A %*% t(X - mu_0))
  else print("Incorrect dimension size")
}


calc_S <- function(U, sigma_0, beta) {
  S <- vector("list", nrow(U))
  S[[1]] <- beta*(U[1, ] %*% t(U[1, ])) + (1-beta)*sigma_0
  for(t in 2:nrow(U)) {
    S[[t]] <- beta*(U[t, ] %*% t(U[t, ])) + (1-beta)*S[[t-1]]
  }
  S
}

sparsify_S <- function(S, lambda_s = .2) {
  lapply(S, \(S) {spcov::spcov(S, S, lambda = lambda_s, step.size = 100)$Sigma })
}

threshold_S <- function(S, cutoff = .1) {
  purrr::map(S, \(X) {X[(abs(X) < cutoff)] <- 0; X})
}

get_U_window <- function(U, n_w = 50) {
  purrr::map(1:nrow(U), \(i) {
    if(i < n_w) NA
    else (i-n_w+1):i
  })
}

get_S_hat <- function(S_star, U_window, U) {
  purrr::map2(S_star, U_window, \(S, i) {
    if(anyNA(i) || is.null(i)) NA
    else covchaud(S, U[i, ])$mat
  })
}

pstat_hawkins <- function(S) {
  p <- nrow(dplyr::last(S))
  map_dfr(S, \(S) {
    if(anyNA(S)) tibble(pstat = NA,
                        term_1 = NA,
                        term_2 = NA)
    else tibble(pstat = sum(diag(S)) - log(det(S)) - p,
                term_1 = sum(diag(S)),
                term_2 = -log(det(S)))
      
    })
}

pstat_test_1 <- function(S) {
  p <- nrow(dplyr::last(S))
  map_dfr(S, \(S) {
    if(anyNA(S)) tibble(pstat = NA,
                        term_1 = NA,
                        term_2 = NA)
    else tibble(pstat = sum(diag(S)) - log(det(S)) - p,
                term_1 = sum(diag(S)),
                term_2 = -log(det(S)))
    
  })
}

pstat_test_2 <- function(mu_t, mu_0, S_hat, sig_inv) {
  p <- length(mu_0)
  T_1 <- apply(mu_t, 1, \(x) {t(x - mu_0) %*% sig_inv %*% (x - mu_0)})
  
  map_dfr(S_hat, \(S) {
    if(anyNA(S)) tibble(T_2 = NA,
                        term_1 = NA,
                        term_2 = NA)
    else tibble(T_2 = sum(diag(S)) - log(det(S)) - p,
                term_1 = sum(diag(S)),
                term_2 = -log(det(S)))
    
  }) |> mutate(T_1 = T_1, .before = T_2)
}

pstat_test_3 <- function(mu_t, mu_0, S_t, S_hat, sig_inv) {
  # mu_t is a list
  p <- length(mu_0)
  
  pmap_dfr(list(x = mu_t, S_t = S_t, S_hat = S_hat), \(x, S_t, S_hat) {
    if(anyNA(S_hat)) tibble(T_1 = NA,
                        T_2 = NA,
                        term_1 = NA,
                        term_2 = NA)
    else tibble(T_1 = t(x - mu_0) %*% solve(S_hat) %*% (x - mu_0),
                T_2 = - log(sum(abs(S_hat))) - 
                  sum(diag(S_t %*% solve(S_hat))) +
                  sum(diag(sig_inv %*% S_t)),
                term_1 = - log(det(S_hat)),
                term_2 = - sum(diag(S_t %*% solve(S_hat))))
  })
}

# Wang Methods ------------------------------------------------------------
mean_MEWMA <- function(dat, mu_0, alpha = .2) {
  X <- as.matrix(dat)
  mu_t <- matrix(nrow = nrow(dat), ncol = length(mu_0))
  
  mu_t[1, ] <- alpha*X[1, ] + (1-alpha)*mu_0
  for(i in 2:nrow(X)) {
    mu_t[i, ] <- alpha*X[i-1, ] + (1-alpha)*mu_t[i-1, ]
  }
  
  mu_t
}

pl_mu_1 <- function(mu, mu_t, sig_inv, lambda_1) {
  (mu-mu_t) %*% sig_inv %*% (mu-mu_t) + lambda_1*sum(abs(mu))
}

mean_SPARSE <- function(mu_t, mu_0, sig_inv, fun, lambda_1 = 0.2) {
  1:nrow(mu_t) |> 
    purrr::map_dfr(\(i){
      optim(mu_0, fun, mu_t = mu_t[i, ], sig_inv = sig_inv, lambda_1 = lambda_1,
            method = "BFGS",
            control = list(maxit = 1000))$par
    })
}

cov_MEWMC <- function(dat, mu_lt, sigma_0, beta = .2) {
  X <- as.matrix(dat)
  mu_lt <- as.matrix(mu_lt)
  sigma_lt <- vector("list", nrow(dat))
  
  sigma_lt[[1]] <- beta*(X[1, ] - mu_lt[1, ]) %*% t(X[1, ] - mu_lt[1, ]) + (1-beta)*sigma_0
  for(i in 2:nrow(X)) {
    sigma_lt[[i]] <- beta*(X[i-1, ] - mu_lt[i-1, ]) %*% t(X[i-1, ] - mu_lt[i-1, ]) + (1-beta)*sigma_lt[[i-1]]
  }
  sigma_lt
}

cov_SPARSE <- function(sigma_lt, S0, lambda_2 = .2) {
  sigma_lt |> 
    map(\(x) {
      spcov::spcov(S0, x, lambda = lambda_2, step.size = 100)$Sigma # Uses spcov, unmodified
    })
}

cov_MEWMC_0 <- function(dat, mu_0, sigma_0, beta = .2) {
  X <- as.matrix(dat)
  sigma_0t <- vector("list", nrow(dat))
  
  sigma_0t[[1]] <- beta*(X[1, ] - mu_0) %*% t(X[1, ] - mu_0) + (1-beta)*sigma_0
  for(i in 2:nrow(dat)) {
    sigma_0t[[i]] <- beta*(X[i, ] - mu_0) %*% t(X[i, ] - mu_0) + (1-beta)*sigma_0t[[i-1]]
  }
  sigma_0t
}

pstat_wang_1 <- function(sigma_hat_lt, sigma_lt, sigma_0t, sig_inv) {
  list(pstat = 
         -log(sum(abs(sigma_hat_lt))) -
         sum(diag(sigma_lt %*% solve(sigma_hat_lt))) +
         sum(diag(sig_inv %*% sigma_0t)),
       term_1 = -log(sum(abs(sigma_hat_lt))),
       term_2 = -sum(diag(sigma_lt %*% solve(sigma_hat_lt))),
       term_3 = sum(diag(sig_inv %*% sigma_0t)))
}

pstat_wang_2 <- function(mu_lt, mu_0, sigma_0t, sigma_hat_0t, sig_inv) {
  list(T_1 = as.numeric(t(mu_lt - mu_0) %*% sig_inv %*% (mu_lt - mu_0)),
       T_2 = 
         -log(sum(abs(sigma_hat_0t))) -
         sum(diag(sigma_0t %*% solve(sigma_hat_0t))) +
         sum(diag(sig_inv %*% sigma_0t)),
       term_1 = -log(sum(abs(sigma_hat_0t))),
       term_2 = -sum(diag(sigma_0t %*% solve(sigma_hat_0t))),
       term_3 = sum(diag(sig_inv %*% sigma_0t)))
}


# Adaptive Window Size ----------------------------------------------------

create_adaptive_window <- function(dat, n_min = 50, n_max = 100,
                                   mu_0 = NA, sigma_0 = NA, 
                                   S = NA, sigma_S_0 = NA, 
                                   upper_limit = NA,
                                   condition = "mean") {
  X <- as.matrix(dat)
  n_t <- list()
  n_t[[n_min]] <- 1:n_min
  test <- TRUE
  
  for(t in (n_min+1):nrow(X)) {
    if(condition == "mean") {
      test <- upper_limit >= as.numeric(t(X[t-1, ] - mu_0) %*% sigma_0 %*% (X[t-1, ] - mu_0))
    } else if(condition == "variance") {
      test <- upper_limit >= as.numeric(t(ks::vech(sigma_0)) %*% sigma_S_0 %*% ks::vech(S[[t-1]]))
    }
    
    if(test) {
      n_t[[t]] <- (t-min(n_max, length(n_t[[t-1]])+1) + 1):t
    } else {
      #n_t[[t]] <- (t-max(n_min, length(n_t[[t-1]])-1) + 1):t
      n_t[[t]] <- (t-n_min+1):t
    }
  }
  
  n_t
}

# sigma_0 <- cov(x_ic)
# mu_0 <- colMeans(x_ic)
# alpha <- .2
# 
# A <- get_A(sigma_0)
# X <- as.matrix(x_ic)
# mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
# U <- calc_U(x_ic, mu_t, A)
# S <- calc_S(U, sigma_0, beta = .2)
# 
# get_U_window(U, n_w = 100)
# create_adaptive_window(x_ic, n_min = 50, n_max = 100)
# 
# 
# map_dbl(S, \(S) {t(ks::vech(sigma_0)) %*% ks::vech(S)})
# 
# S_hat <- get_S_hat(threshold_S(S, cutoff = .1),
#                    get_U_window(U, n_w = 100),
#                    U)
# pstat_test_2(mu_t, mu_0, S_hat, solve(sigma_0))
# 
# vec <- apply(X-mu_0, MARGIN = 1, FUN = \(x) t(x) %*% sigma_0 %*% x)
# 
# h_arl <- estimate_h(vec, 200)
# h_mrl <- quantile(vec, 1-1/200)
# 
# a <- create_adaptive_window(x_ic, mu_0 = mu_0, sigma_0 = sigma_0, upper_limit = h_arl,
#                             condition = "mean")
# 
# sigma_S_0 <- cov(matrix(map(S, ks::vech) |> unlist(), ncol = p*(p-1), byrow = TRUE))
# 
# vec <- map_dbl(S, \(S) {t(ks::vech(sigma_0)) %*% sigma_S_0 %*% ks::vech(S)})
# 
# h_arl <- estimate_h(vec, 200)
# h_mrl <- quantile(vec, 1-1/200)
# 
# b <- create_adaptive_window(x_ic, mu_0 = mu_0, sigma_0 = sigma_0, 
#                             S = S, sigma_S_0 = sigma_S_0, upper_limit = h_arl,
#                             condition = "variance")
# 
# map_dbl(a, length) |> plot()
# map_dbl(b, length) |> plot()


# create_adaptive_window(x_ic)[[106]]
# 
# n_t <- vector("integer", nrow(x_ic))
# 
# 
# X <- as.matrix(x_ic, sig_inv)
# 
# Z <- vector("numeric", lengch = nrow(X))
# Z[1] <- 0
# 
# for(i in 2:nrow(X)) {
#   Z[i] <- (t(X[i, ]) %*% sig_inv %*% (X[i-1, ]))
# }
# 
# plot(Z)