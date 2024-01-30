library("tidyverse")
library("mvtnorm")
library("sn")
library("here")
source(here("scripts", "pl_fd_methods_v01.R"))
# Setup for Data Generation -----------------------------------------------
n <- 1000

p1 <- 3
p2 <- 10

mu1_0 <- c(0, 0, 0)
mu1_1 <- c(.2, .2, 0)
mu2_0 <- rep(0, p2)
mu2_1 <- c(.2, .2, rep(0, p2 - 2))

sigma1_0 <- 
  matrix(c(1, 0, 0,
           0, 1, 0,
           0, 0, 1), ncol = p1, byrow = TRUE)

sigma1_1 <- 
  matrix(c(1.5, 0, 0,
           0, 1.5, 0,
           0, 0, 1.5), ncol = p1, byrow = TRUE)

sigma1_2 <- 
  matrix(c(1, 0.5, 0,
           0.5, 1, 0,
           0, 0, 1), ncol = p1, byrow = TRUE)



sigma2_0 <- 
  matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ncol = p2, byrow = TRUE)

sigma2_1 <- 
  matrix(c(1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1.5, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 1.5, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ncol = p2, byrow = TRUE)

sigma2_2 <- 
  matrix(c(1, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0,
           0.5, 1, 0.5, 0, 0, 0, 0, 0, 0, 0,
           0.5, 0.5, 1, 0, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ncol = p2, byrow = TRUE)

alpha1 <- rep(10, p1)
alpha2 <- rep(10, p2)
nu1 <- 40


# Data Generation ---------------------------------------------------------
set.seed(321)

x1_ic <- rmvnorm(n, mu1_0, sigma1_0) |> as_tibble()
x1_f1 <- rmvnorm(n, mu1_1, sigma1_0) |> as_tibble()
x1_f2 <- rmvnorm(n, mu1_0, sigma1_1) |> as_tibble()
x1_f3 <- rmvnorm(n, mu1_0, sigma1_2) |> as_tibble()

x2_ic <- rmvnorm(n, mu2_0, sigma2_0) |> as_tibble()
x2_f1 <- rmvnorm(n, mu2_1, sigma2_0) |> as_tibble()
x2_f2 <- rmvnorm(n, mu2_0, sigma2_1) |> as_tibble()
x2_f3 <- rmvnorm(n, mu2_0, sigma2_2) |> as_tibble()

x3_ic <- rmst(n, mu1_0, sigma1_0, alpha1, nu1) |> as_tibble()
x3_f1 <- rmst(n, mu1_1, sigma1_0, alpha1, nu1) |> as_tibble()
x3_f2 <- rmst(n, mu1_0, sigma1_1, alpha1, nu1) |> as_tibble()
x3_f3 <- rmst(n, mu1_0, sigma1_2, alpha1, nu1) |> as_tibble()

x4_ic <- rmst(n, mu2_0, sigma2_0, alpha2, nu1) |> as_tibble()
x4_f1 <- rmst(n, mu2_1, sigma2_0, alpha2, nu1) |> as_tibble()
x4_f2 <- rmst(n, mu2_0, sigma2_1, alpha2, nu1) |> as_tibble()
x4_f3 <- rmst(n, mu2_0, sigma2_2, alpha2, nu1) |> as_tibble()

fault_names <- 
  c("p=3, change in mean",
    "p=3, change in variance",
    "p=3, change in covariance",
    "p=10, change in mean",
    "p=10, change in variance",
    "p=10, change in covariance",
    "skew t, p=3, change in mean",
    "skew t, p=3, change in variance",
    "skew t, p=3, change in covariance",
    "skew t, p=10, change in mean",
    "skew t, p=10, change in variance",
    "skew t, p=10, change in covariance")

df_faults <- 
  list(bind_rows(x1_ic, x1_f1),
       bind_rows(x1_ic, x1_f2),
       bind_rows(x1_ic, x1_f3),
       bind_rows(x2_ic, x2_f1),
       bind_rows(x2_ic, x2_f2),
       bind_rows(x2_ic, x2_f3),
       bind_rows(x3_ic, x3_f1),
       bind_rows(x3_ic, x3_f2),
       bind_rows(x3_ic, x3_f3),
       bind_rows(x4_ic, x4_f1),
       bind_rows(x4_ic, x4_f2),
       bind_rows(x4_ic, x4_f3)) |> set_names(fault_names)

hawkins_results <- 
  map(df_faults, \(x) {
    hawkins_2008(x, colMeans(x[1:n, ]), cov(x[1:n, ]))
  })

hawkins_results_plots <- 
  map2(hawkins_results, fault_names, \(x, plot_name) {
    x |> 
      ggplot(aes(1:length(pstat), pstat)) +
      geom_point(shape = 21) +
      labs(x = "index", title = plot_name)
  })

hawkins_results_plots
