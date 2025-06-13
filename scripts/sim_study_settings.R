# Parameters selection for each scenario ----------------------------------
a   <- c(0, .3, 1, 2)
b   <- c(1, 1.3, 1.8, 2.5)
rho <- c(0, .2, .4, .6)

scenario_params <- list(s1 = a,
                        s2 = b,
                        s3 = rho,
                        s4 = b,
                        s5 = a,
                        s6 = b,
                        s7 = rho,
                        s8 = rho,
                        s9 = a,
                        s10 = b)

method_params <- list(
  hawkins   = list(method = "hawkins",    beta = .1),
  wang_1    = list(method = "wang_1",    alpha = .1,     beta = .1, lambda_1 = .2, lambda_2 = .2),
  wang_2    = list(method = "wang_2",    alpha = .1,     beta = .1, lambda_1 = .2, lambda_2 = .2),
  MC_LASSO  = list(method = "MC_LASSO",   beta = .1, lambda_s = .2),
  MC_COMET  = list(method = "MC_COMET",   beta = .1,   cutoff = .1, n_w = c(60, 60)),
  MAC_COMET = list(method = "MAC_COMET", alpha = .1,     beta = .1, cutoff = .1, n_w = c(60, 60))
)

data_folder <- "data_new_nu"

methods <- "wang_2" 
#methods <- c("hawkins", "MC_COMET", "MAC_COMET")
#methods <- c("wang_1", "wang_2", "MC_LASSO")

n <- 1000
arl_ic <- 180
i_min <- 1
i_max <- 1000