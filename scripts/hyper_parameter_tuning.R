
# Tuning parameters -------------------------------------------------------

# Hawkins: beta
# wang_2014: alpha, beta, lambda_1, lambda_2
# MC_LASSO: beta, lambda_s
# MC_COMET: beta, cutoff, n_w
# MAC_COMET: alpha, beta, cutoff, n_w

# Of primary interest is a method for tuning MAC_COMET.
# alpha \in (0, 1], determines memory length of mean estimate.
# beta \in (0, 1], determines memory length of variance estimate.
# cutoff \in \R^+, determines which values in Sigma are set to 0.
# n_w \in \Z^+, determines number of observations used in variance estimation.


# Tuning alpha and beta ---------------------------------------------------



# Tuning COMET ------------------------------------------------------------

# Test Frobenius norm for n_w at different levels.
# We want to combine accuracy quick response times.
n_w <- c(20, 40, 60, 80, 100)
cutoff <- c(.01, .05, .1, .2, .4)

params <- expand_grid(n_w, cutoff)

tuning_results <- generate_n_samples(10)

summarize_tuning_results <- tuning_results |> 
  summarize(across(c(s3, s7, s12), mean), .by = c(n_w, cutoff))

tuning_results <- bind_cols(params, tuning_results)

tuning_results |> 
  pivot_longer(s3:s12, names_transform = as.factor) |> 
  mutate(cutoff = as.factor(cutoff)) |> 
  ggplot(aes(n_w, value, color = cutoff)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  labs(y = "Frobenius Norm")


# Helper functions --------------------------------------------------------

# Frobenius norm
f_norm <- function(S1, S2) {
  sum(sqrt((S1 - S2)^2))
}

# Test fit
test_fit <- function(n_w, cutoff, scenario) {
  if(scenario == "s3") {
    df <- gen_dat_s3(n_w, .2)[1:n_w, ]
    S2 <- sigma_1(.2)
  } else if(scenario == "s7") {
    df <- gen_dat_s7(n_w, .2)[1:n_w, ]
    S2 <- sigma_2(1, .2)
  } else if(scenario == "s12") {
    df <- gen_dat_s12(n_w, .2)[1:n_w, ]
    S2 <- sigma_2_20(1, .2)
  }
  
  X <- as.matrix(df)
  S <- var(X)
  S[abs(S) < cutoff] <- 0
  
  S1 <- covchaud(S, X)$mat
  
  f_norm(S1, S2)
}

# Generate tuning data
generate_n_samples <- function(n) {
  1:n |> 
    map_dfr(\(i) {
      tibble(params,
             params |> 
               pmap_dfr(\(n_w, cutoff) {
                 gen_dat_s3(n_w, .2)[1:n_w, ]
                 gen_dat_s7(n_w, .2)[1:n_w, ]
                 gen_dat_s12(n_w, .2)[1:n_w, ]
                 
                 tibble(s3 = test_fit(n_w, cutoff, "s3"),
                        s7 = test_fit(n_w, cutoff, "s7"),
                        s12= test_fit(n_w, cutoff, "s12"))
               }))
    })
}



