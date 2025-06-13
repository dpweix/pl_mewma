library("tidyverse")
library("here")
source(here("scripts", "data_generation.R"))
source(here("scripts", "pl_fd_methods.R"))

# Given a single scenario/method, saves results for one run of each scenario

# scenario <- "s1"
# i        <- 7
# arl_ic   <- 200
# n        <- 1000
# arg      <- 1
# 
# params_method <- list(method = "hawkins", lambda = .2)

# Generate data -----------------------------------------------------------
set.seed(i)
df <- gen_dat(scenario, n, arg, data_only = TRUE)

# Calculate and pstat -----------------------------------------------------
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
  
} else if(method == "wang_2") {
  mu_0    <- colMeans(df[1:n, ])
  sigma_0 <- cov(df[1:n, ])
  
  pstat <- wang_2014_2(df, mu_0, sigma_0, 
                       alpha = selected_method$alpha,
                       beta = selected_method$beta,
                       lambda_1 = selected_method$lambda_1,
                       lambda_2 = selected_method$lambda_2)
  
  T_1 <- pstat$T_1
  T_2 <- pstat$T_2
  
} else if(method == "MC_LASSO") {
  mu_0    <- colMeans(df[1:n, ])
  sigma_0 <- cov(df[1:n, ])#spcov(cov(df[1:n, ]), cov(df[1:n, ]),
                   #lambda = selected_method$lambda_s, step.size = 100)$Sigma
  
  pstat <- MC_LASSO(df, mu_0, sigma_0,
                    beta = selected_method$beta,
                    lambda_s = selected_method$lambda_s)$pstat
  
} else if(method == "MC_COMET") {
  mu_0    <- colMeans(df[1:n, ])
  #S_0 <- cov(df[1:n, ])
  #S_0[(abs(S_0) < selected_method$cutoff)] <- 0
  #sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
  sigma_0 <- cov(df[1:n, ])
  
  pstat <- MC_COMET(df, mu_0, sigma_0,
                    beta = selected_method$beta,
                    cutoff = selected_method$cutoff,
                    n_w = ifelse(scenario %in% paste0("s", 1:5),
                                 selected_method$n_w[1],
                                 selected_method$n_w[2]))$pstat
  
} else if(method == "MAC_COMET") {
  mu_0    <- colMeans(df[1:n, ])
  #S_0 <- cov(df[1:n, ])
  #S_0[(abs(S_0) < selected_method$cutoff)] <- 0
  #sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
  sigma_0 <- cov(df[1:n, ])
  
  pstat <- MAC_COMET(df, mu_0, sigma_0,
                     alpha = selected_method$alpha,
                     beta = selected_method$beta,
                     cutoff = selected_method$cutoff,
                     n_w = ifelse(scenario %in% paste0("s", 1:5),
                                  selected_method$n_w[1],
                                  selected_method$n_w[2]))
  T_1 <- pstat$T_1
  T_2 <- pstat$T_2
} else if(method == "MAC_COMET1") {
  mu_0    <- colMeans(df[1:n, ])
  S_0 <- cov(df[1:n, ])
  S_0[(abs(S_0) < selected_method$cutoff)] <- 0
  sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
  
  pstat <- MAC_COMET1(df, mu_0, sigma_0,
                     alpha = selected_method$alpha,
                     beta = selected_method$beta,
                     cutoff = selected_method$cutoff,
                     n_w = ifelse(scenario %in% paste0("s", 1:5),
                                  selected_method$n_w[1],
                                  selected_method$n_w[2]))
  
  T_1 <- pstat$T_1
  T_2 <- pstat$T_2
} 

# Save results -------------------------------------------------------------
if(method %in% c("hawkins", "wang_1", "MC_LASSO", "MC_COMET")) {
  h  <- estimate_h(pstat[1:n], arl_ic)
  rl <- get_rl_oc(pstat[(n+1):(2*n)], h = h)
  
  h_B  <- estimate_h_bootstrap(pstat[1:n], arl_ic)
  rl_B <- get_rl_oc(pstat[(n+1):(2*n)], h = h)
  
  list(scenario = scenario,
       method = method,
       pstat = pstat,
       params_data = list(n = n, arl_ic = arl_ic),
       params_method = selected_method,
       h = h,
       rl = rl,
       h_B = h_B,
       rl_B = rl_B) |>
    saveRDS(file = here(data_folder, paste0(scenario, "-", method, "-", arg, "-", i, ".rds")))
  
} else if(method %in% c("wang_2", "MAC_COMET", "MAC_COMET1")) {
  h1  <- estimate_h(T_1[1:n], arl_ic)
  h2  <- estimate_h(T_2[1:n], arl_ic)
  
  rl1 <- get_rl_oc(T_1[(n+1):(2*n)], h = h1)
  rl2 <- get_rl_oc(T_2[(n+1):(2*n)], h = h2)
  
  rl <- min(rl1, rl2)
  
  h1_B  <- estimate_h_bootstrap(T_1[1:n], arl_ic)
  h2_B  <- estimate_h_bootstrap(T_2[1:n], arl_ic)
  
  rl1_B <- get_rl_oc(T_1[(n+1):(2*n)], h = h1)
  rl2_B <- get_rl_oc(T_2[(n+1):(2*n)], h = h2)
  
  rl_B <- min(rl1_B, rl2_B)
  
  list(scenario = scenario,
       method = method,
       pstat = tibble(T_1 = T_1, T_2 = T_2),
       params_data = list(n = n, arl_ic = arl_ic),
       params_method = selected_method,
       h1 = h1,
       h2 = h2,
       rl = rl,
       h1_B = h1_B,
       h2_B = h2_B,
       rl_B = rl_B) |> 
    saveRDS(file = here(data_folder, paste0(scenario, "-", method, "-", arg, "-", i, ".rds")))
}