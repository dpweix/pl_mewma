# Parameters selection for each scenario ----------------------------------
a   <- c(0, .3, 1, 2)
b   <- c(1, 1.3, 1.8, 2.5)
rho <- c(0, .2, .4, .6)

alpha <- c(0.05, 0.1, 0.2)
beta  <- alpha

scenario_params <- list(s1 = a,
                        s2 = b)#,
                        # s3 = rho,
                        # s4 = b,
                        # s5 = b,
                        # s6 = a,
                        # s7 = b,
                        # s8 = rho,
                        # s9 = a,
                        # s10 = b)
method_params <- list(
  hawkins = list(method = "hawkins",    beta = .2),
  wang_1  = list(method = "wang_1",    alpha = .2,     beta = .2, lambda_1 = .2, lambda_2 = .2),
  wang_2  = list(method = "wang_2",    alpha = .2,     beta = .2, lambda_1 = .2, lambda_2 = .2),
  test_1  = list(method = "MC_LASSO",   beta = .2, lambda_s = .2),
  test_2  = list(method = "MC_COMET",   beta = .2,   cutoff = .1, n_w = 60),
  test_3  = list(method = "MAC_COMET", alpha = .2,     beta = .2, cutoff = .1, n_w = 60)
)

# data1 is testing methods with new names
# data2 is testing methods with variable alpha, beta
# data3 is testing methods with variable alpha, beta, and n_w (short long) 

# p =  3: n_w = 30,  60
# p = 10: n_w = 60, 120

data_folder <- "data1"

method <- "wang_2"
selected_method <- method_params[[method]]

n <- 1000
arl_ic <- 200
i_min <- 1
i_max <- 20
