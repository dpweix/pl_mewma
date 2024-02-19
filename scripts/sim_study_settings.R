# Parameters selection for each scenario ----------------------------------
a   <- c(0, .3, 1, 2)
b   <- c(1, 1.3, 1.8, 2.5)
rho <- c(0, .2, .4, .6)

#alpha <- c(0.05, 0.1, 0.2)
#beta  <- alpha

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
  #MAC_COMET1= list(method ="MAC_COMET1", alpha = .1,     beta = .1, cutoff = .1, n_w = c(60, 60))
)

# data1 is testing methods with new names
# data2 is testing methods with alpha, beta = .1 and n_w = 275, arl = 180
# data3 is testing methods with alpha, beta = .01 for Wang and .1 fro COMET, n_w = 165, arl = 180
# data4 is data3 bu n_w = 18 for s1:5 and 165 for s6:10
# data5 is bootstrapped h estimation
# data6 is get_arl_ic for case of IC ARL estimation
# data_test1 is n_w = 60 for both large/small

# p =  3: n_w = 30,  60
# p = 10: n_w = 60, 120

data_folder <- "data"

methods <- "MAC_COMET"
#methods <- c("hawkins", "wang_1", "wang_2")
#methods <- c("wang_1", "wang_2", "MC_LASSO")
#selected_method <- method_params[[methods]]

n <- 1000
arl_ic <- 180
i_min <- 1
i_max <- 200