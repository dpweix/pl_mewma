library("tidyverse")
theme_set(cowplot::theme_cowplot())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.placement = "outside",
             strip.background = element_blank())
library("mvtnorm")
library("spcov")
library("CovTools")
library("mgcov")

# Li 2013, Monitoring Covariance ------------------------------------------
n <- 1000
p <- 3
n_w <- 30

# Generate sample data
set.seed(321)
x_ic <- rmvnorm(n, rep(0, p),   diag(p)) |> as_tibble()
x_f1 <- rmvnorm(n, rep(0, p), 2*diag(p)) |> as_tibble()

# Sparse Estimation
omega_l <- CovTools::PreEst.glasso(as.matrix(x_ic), method = list(type = "fixed", param = .2))$C
omega_l |> round(2)

# Charting statistic
# Works by considering the sample covariance matrix from a moving
# window of size n_w
charting_stat1 <- function(dat, omega_l) {
  S <- cov(as.matrix(dat))
  sum(diag(S)) - sum(diag(omega_l %*% S)) + log(det(omega_l))
}

# Separate data into n_w size chunks
rolling_sep <- function(dat, n_w) {
  n_w:nrow(dat) |> 
    purrr::map(\(i) {
      dat[(i-n_w+1):i, ]
    })
}

x_split <- rolling_sep(bind_rows(x_ic, x_f1), n_w)

# Find the charting stat for each group
c_stat <- x_split |> 
  map_dbl(\(x) {
    charting_stat1(x, omega_l)
  })

tibble(c_stat = c_stat,
       index = 1:length(c_stat)) |> 
  ggplot(aes(index, c_stat)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 970, color = "red", linetype = "dashed")


# Wang, Yeh, Li, 2014, Monitoring Mean and Covariance ---------------------
n <- 1000
p <- 3
n_w <- 30

mat_f4 <- matrix(c(1, 0, .5, 0, 1, 0, .5, 0, 1), byrow = TRUE, nrow = 3)

# Generate sample data (IC, var shift, mean shift, mean + var shift)
set.seed(321)
x_ic <- rmvnorm(n, rep(0, p), diag(p))   |> as_tibble()
x_f1 <- rmvnorm(n, rep(0, p), 2*diag(p)) |> as_tibble()
x_f2 <- rmvnorm(n, rep(1, p), diag(p))   |> as_tibble()
x_f3 <- rmvnorm(n, rep(1, p), 2*diag(p)) |> as_tibble()
x_f4 <- rmvnorm(n, rep(0, p), mat_f4)    |> as_tibble()

### Chart 1 ###
mu_0 <- map_dbl(x_ic, mean)
sigma_0 <- cov(as.matrix(x_ic))
sig_inv <- solve(sigma_0)

# MEWMA mean estimate
mean_MEWMA <- function(dat, mu_0, alpha = .2) {
  X <- as.matrix(dat)
  mu_t <- matrix(nrow = nrow(dat), ncol = length(mu_0))
  
  mu_t[1, ] <- alpha*X[1, ] + (1-alpha)*mu_0
  for(i in 2:nrow(X)) {
    mu_t[i, ] <- alpha*X[i-1, ] + (1-alpha)*mu_t[i-1, ]
  }
  mu_t
}

mu_t <- mean_MEWMA(x_ic, mu_0)

# Sparse mean estimation
# PL mean estimation 1
pl_mu_1 <- function(mu, mu_t, sig_inv, lambda) {
  (mu-mu_t) %*% sig_inv %*% (mu-mu_t) + lambda*sum(abs(mu))
}

mean_SPARSE <- function(mu_t, mu_0, sig_inv, fun, lambda = 0.2) {
  1:nrow(mu_t) |> 
    map_dfr(\(i){
      optim(mu_0, fun, mu_t = mu_t[i, ], sig_inv = sig_inv, lambda = lambda,
            method = "BFGS",
            control = list(maxit = 1000))$par
    })
}

mu_lt <- mean_SPARSE(mu_t, mu_0, sig_inv, pl_mu_1)

# MEWMC covariance estimate
cov_MEWMC <- function(dat, mu_lt, omega = .2) {
  X <- as.matrix(dat)
  mu_lt <- as.matrix(mu_lt)
  sigma_lt <- vector("list", nrow(dat))
  
  sigma_lt[[1]] <- omega*(X[1, ] - mu_lt[1, ]) %*% t(X[1, ] - mu_lt[1, ]) + (1-omega)*sigma_0
  for(i in 2:nrow(X)) {
    sigma_lt[[i]] <- omega*(X[i-1, ] - mu_lt[i-1, ]) %*% t(X[i-1, ] - mu_lt[i-1, ]) + (1-omega)*sigma_lt[[i-1]]
  }
  sigma_lt
}

sigma_lt <- cov_MEWMC(x_ic, mu_lt)

sigma_0

# For sparse covariance estimation use
cov_SPARSE <- function(sigma_lt, S0, lambda = .2) {
  sigma_lt |> 
    map(\(x) {
      spcov::spcov(S0, x, lambda = .2, step.size = 100)$Sigma #Uses spcov, unmodified
    })
}

cov_SPARSE_COMET

sigma_hat_lt <- cov_SPARSE(sigma_lt, sigma_0)

# MEWMC covariance estimate using mu_0 instead of mu_lt
cov_MEWMC_0 <- function(dat, mu_0, sigma_0, omega = .2) {
  X <- as.matrix(dat)
  sigma_0t <- vector("list", nrow(dat))
  
  sigma_0t[[1]] <- omega*(X[1, ] - mu_0) %*% t(X[1, ] - mu_0) + (1-omega)*sigma_0
  for(i in 2:nrow(dat)) {
    sigma_0t[[i]] <- omega*(X[i, ] - mu_0) %*% t(X[i, ] - mu_0) + (1-omega)*sigma_0t[[i-1]]
  }
  sigma_0t
}

sigma_0t <- cov_MEWMC_0(x_ic, mu_0, sigma_0)
sigma_hat_0t <- cov_SPARSE(sigma_0t, sigma_0)

# Plotting statistic
# Given the following matrices we can calc the pstat
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

pstat_wang_1(sigma_hat_lt[[1]], sigma_lt[[1]], sigma_0t[[1]], sig_inv)

pstat_1 <- 
  pmap_dbl(list(sigma_hat_lt, sigma_lt, sigma_0t, list(sig_inv)), pstat_wang_1)

plot(pstat_1)

# Apply method
wang_2014_1 <- function(dat, mu_0, sigma_0,
                        alpha = .2, omega = .2,
                        lambda1 = .2, lambda2 = .2) {
  X <- as.matrix(dat)
  sig_inv <- solve(sigma_0)
  
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  mu_lt <- mean_SPARSE(mu_t, mu_0, sig_inv, pl_mu_1, lambda = lambda1)
  sigma_lt     <- cov_MEWMC(X, mu_lt, omega = omega)
  sigma_hat_lt <- cov_SPARSE(sigma_lt, sigma_0, lambda = lambda2)
  sigma_0t     <- cov_MEWMC_0(X, mu_0, sigma_0, omega = omega)
  
  pmap_dfr(list(sigma_hat_lt, sigma_lt, sigma_0t, list(sig_inv)), pstat_wang_1)
  
}

wang_2014_2 <- function(dat, mu_0, sigma_0,
                        alpha = .2, omega = .2,
                        lambda1 = .2, lambda2 = .2) {
  X <- as.matrix(dat)
  sig_inv <- solve(sigma_0)
  
  mu_t  <- mean_MEWMA(X, mu_0, alpha = alpha)
  mu_lt <- mean_SPARSE(mu_t, mu_0, sig_inv, pl_mu_1, lambda = lambda1)
  sigma_0t     <- cov_MEWMC_0(X, mu_0, sigma_0, omega = omega)
  sigma_hat_0t <- cov_SPARSE(sigma_0t, sigma_0, lambda = lambda2)
  
  pmap_dfr(list(pmap(mu_lt, c), list(mu_0), sigma_0t, sigma_hat_0t, list(sig_inv)), 
           pstat_wang_2)
}

pmap_dfr(list(sigma_hat_lt, sigma_lt, sigma_0t, list(sig_inv)), pstat_wang_1)
pmap_dfr(list(pmap(mu_lt, c), list(mu_0), sigma_0t, sigma_hat_0t, list(sig_inv)), pstat_wang_2)

pstat_wang_1_f1 <- wang_2014_1(bind_rows(x_ic, x_f1), mu_0, sigma_0)
pstat_wang_1_f2 <- wang_2014_1(bind_rows(x_ic, x_f2), mu_0, sigma_0)
pstat_wang_1_f3 <- wang_2014_1(bind_rows(x_ic, x_f3), mu_0, sigma_0)
pstat_wang_1_f4 <- wang_2014_1(bind_rows(x_ic, x_f4), mu_0, sigma_0)

pstat_wang_2_f1 <- wang_2014_2(bind_rows(x_ic, x_f1), mu_0, sigma_0)
pstat_wang_2_f2 <- wang_2014_2(bind_rows(x_ic, x_f2), mu_0, sigma_0)
pstat_wang_2_f3 <- wang_2014_2(bind_rows(x_ic, x_f3), mu_0, sigma_0)
pstat_wang_2_f4 <- wang_2014_2(bind_rows(x_ic, x_f4), mu_0, sigma_0)

df_wang_1 <-
  bind_rows(list(f1 = pstat_wang_1_f1, 
                 f2 = pstat_wang_1_f2, 
                 f3 = pstat_wang_1_f3,
                 f4 = pstat_wang_1_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

df_wang_2 <-
  bind_rows(list(f1 = pstat_wang_2_f1, 
                 f2 = pstat_wang_2_f2, 
                 f3 = pstat_wang_2_f3,
                 f4 = pstat_wang_2_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

# Chart 1 Plots
df_wang_1 |> 
  ggplot(aes(index, pstat)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x")

df_wang_1 |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1) + 
  labs(y = "", color = "")

# Chart 2 Plots
df_wang_2 |> 
  pivot_longer(contains("T_"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x")

df_wang_2 |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1) + 
  labs(y = "", color = "")

