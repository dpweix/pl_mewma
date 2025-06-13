library("tidyverse")
library("here")

source(here("scripts", "data_generation.R"))
source(here("scripts", "pl_fd_methods.R"))

n <- 1000
b <- 1.8

# Generate data
df_3 <- gen_dat("s4", n, b)
df_10<- gen_dat("s6", n, b)

# Estimate params from IC data
mu0_3    <- colMeans(df_3[1:(n/4), ])
sigma0_3 <- cov(df_3[1:(n/4), ])

mu0_10    <- colMeans(df_10[1:(n/4), ])
sigma0_10 <- cov(df_10[1:(n/4), ])


# p = 3: Check for each method
time_1 <- Sys.time()
MC_LASSO(df_3, mu0_3, sigma0_3)
time_2 <- Sys.time()

time_2 - time_1

# p = 10: Check for each method
time_1 <- Sys.time()
MC_LASSO(df_10, mu0_10, sigma0_10)
time_2 <- Sys.time()

time_2 - time_1
