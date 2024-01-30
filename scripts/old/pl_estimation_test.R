library("tidyverse")
library("here")
library("spcov")
library("mgcov")


# Helper functions --------------------------------------------------------
source(here("scripts", "pl_fd_helper.R"))
source(here("scripts", "data_generation.R"))

frobenius_loss <- function(S1, S2) {
  sqrt(sum((S1 - S2)^2))
}

entropy_loss <- function(S1, S2) {
  S1_inv <- solve(S1)
  -log(det(S2 %*% S1_inv)) + sum(diag(S2 %*% S1_inv)) - nrow(S1)
}

run_case <- function(sigma_function, rho, n_w, beta, n_reps, loss = "frobenius") {
  map_dfr(rho, \(rho) {
    map2_dfr(n_w, beta, \(n_w, beta) {
      sigma_0 <- sigma_function(rho)
      map_dfr(1:n_reps, \(i) {
        X <- rmvnorm(n_w, sigma = sigma_0)
        sigma_cov   <- cov(X)
        sigma_mc    <-
          calc_U(X, rep(0, ncol(X)), diag(ncol(X))) |>
          calc_S(sigma_0 = diag(ncol(X)), beta = beta) |> 
          last()
        sigma_lasso <- spcov(sigma_cov, sigma_cov, lambda = .2, step.size = 100)$Sigma
        sigma_comet <- COmet(X, lambda = c(0.05, .1, .2, .3))
        sigma_comet <- sigma_comet$cov_list[[which.min(sigma_comet$bic)]]
        
        tibble(rho = rho, 
               n_w = n_w,
               sample_cov = frobenius_loss(sigma_0, sigma_cov),
               MEWMC      = frobenius_loss(sigma_0, sigma_mc),
               LASSO      = frobenius_loss(sigma_0, sigma_lasso),
               COMET      = frobenius_loss(sigma_0, sigma_comet))
      })
    })
  })
}

summarise_case <- function(df) {
  df |> 
    pivot_longer(-c("rho", "n_w"), names_transform = as_factor) |> 
    summarise(mean_frob_loss = mean(value), .by = c(n_w, name, rho))
}

plot_summary <- function(df) {
  df |> 
    ggplot(aes(name, mean_frob_loss)) +
    geom_bar(stat = "identity") +
    facet_grid(rho ~ n_w)
}

S1 <- sigma_2(1, 0)
X <- rmvnorm(100, sigma = S1)
S2 <- cov(X)
S3 <- spcov(S2, S2, lambda = .2, step.size = 100)$Sigma
S4 <- COmet(X, lambda = .2)$cov_list[[1]]
S5 <-
  calc_U(X, rep(0, ncol(X)), diag(ncol(X))) |>
  calc_S(sigma_0 = diag(ncol(X)), beta = .01) |> 
  last()

frobenius_loss(S1, S2)
frobenius_loss(S1, S3)
frobenius_loss(S1, S4)
frobenius_loss(S1, S5)

entropy_loss(S1, S2)
entropy_loss(S1, S3)
entropy_loss(S1, S4)
entropy_loss(S1, S5)

# Covariances to estimate -------------------------------------------------
source(here("scripts", "data_generation.R"))
source(here("scripts", "sim_study_settings.R"))

# Basically make a test and evaluate at n_w for multiple n_w. 
# See how Frobenius loss changes as n_w increases

rho <- c(0, rho)
n_w <- c(1, 3, 5, 10)*(3*(3+1)/2)
n_w <- c(1, 3, 5, 10)*(10*(10+1)/2)
beta <- c(0.01, .1, .2, .4) |> rev()

# Case: p =  3, sigma_1, rho = 0, .2, .4, .6, .8
df1     <- run_case(sigma_1, rho, n_w, beta, n_reps = 20)
df1_sum <- summarise_case(df1)
df1_sum |> plot_summary()

# Case: p = 10, sigma_2, rho = 0, .2, .4, .6, .8
df2     <- run_case(\(x) {sigma_2(1, x)}, rho, n_w, beta, n_reps = 20)
df2_sum <- summarise_case(df2)
df2_sum |> plot_summary()

# Case: p = 10, sigma_3, rho = 0, .2, .4, .6, .8
df3     <- run_case(sigma_3, rho, n_w, beta, n_reps = 20)
df3_sum <- summarise_case(df3)
df3_sum |> plot_summary()

# p = 3
sigma_0 <- diag(3)
sigma_1(rho)
# p = 10
diag(10)
sigma_2
sigma_3

p <- 3
p <- 10

(p * (p+1))/2
