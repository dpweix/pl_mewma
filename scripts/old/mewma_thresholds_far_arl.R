mewma <- function(train_data, test_data,
                  lambda = 0.1, ic_arl = 200) {
  
  mean_vec <- colMeans(train_data)
  cov_mat <- cov(train_data)
  
  train_x <- scale(train_data, center = mean_vec, scale = F)
  # center the testing (phase 2) data using the parameter estimates from the training (phase 1) data
  x <- scale(test_data, center = mean_vec, scale = F)
  
  p <- ncol(x)
  q_train <- matrix(nrow = nrow(train_x), ncol = p)
  # monitoring statistics during the training period
  for(i in 1:nrow(train_x)) {
    if(i == 1) {
      q_train[i, ] <- lambda * train_x[i, ] + (1 - lambda)*0 #q_0 is the true IC mean, usually the 0 vector (without loss of generality for simulation study)
    } else {
      q_train[i, ] <- lambda * train_x[i, ] + (1 - lambda)*q_train[i-1, ]
    }
  }
  sigma_q <- (lambda / (2 - lambda)) * cov_mat # asymptotic variance, used in Lowry et al. (1992)
  
  t2_train <- numeric()
  for(i in 1:nrow(train_x)) {
    t2_train[i] <- t(q_train[i, ]) %*% solve(sigma_q) %*% q_train[i, ]
  }
  
  # monitoring statistics during the testing period
  q <- matrix(nrow = nrow(x), ncol = p)
  for(i in 1:nrow(x)) {
    if(i == 1) {
      q[i, ] <- lambda * x[i, ] + (1 - lambda)*0 #q_0 is the true IC mean, usually the 0 vector (without loss of generality for simulation study)
    } else {
      q[i, ] <- lambda * x[i, ] + (1 - lambda)*q[i-1, ]
    }
    # sigma_q <- ((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2*i))) * cov_mat # exact covariance
  }
  
  t2_val <- numeric()
  for(i in 1:nrow(x)) {
    t2_val[i] <- t(q[i, ]) %*% solve(sigma_q) %*% q[i, ]
  }
  
  # use a control limit (h) that ensures IC ARL of ic_arl (such as 200)
  h <- spc::mewma.crit(l = lambda, L0 = ic_arl, p = p) # this assumes multivariate normality
  far_h <- quantile(t2_train, probs = 1 - (1/ic_arl))
  
  list(train_stats = t2_train, mon_stats = t2_val,
       arl_exceed = t2_val > h, arl_threshold = h,
       far_exceed = t2_val > far_h, far_threshold = far_h)
}

set.seed(1234)
sim <- function(n) {
  # Generate IC training data
  train_data <- MASS::mvrnorm(5000, mu = rep(0, 3), Sigma = diag(3))
  # Test on IC data
  test_data <- MASS::mvrnorm(2000, mu = rep(0, 3), Sigma = diag(3))
  
  mewma_fit <- mewma(train_data, test_data, lambda = 0.1, ic_arl = 200)
  
  arl_FAR <- mean(mewma_fit$arl_exceed)
  far_FAR <- mean(mewma_fit$far_exceed)
  tibble(phi = arima(mewma_fit$train_stats, order = c(1, 0, 0))$coef[1],
       arl_FAR = arl_FAR, # FAR if the control limit is fixed to ensure a 200 IC ARL
       far_FAR = far_FAR, # FAR if the control limit is fixed to ensure a 0.005 FAR. This value should be 0.005.
       arl_RL = ifelse(arl_FAR == 0, 2001, which(mewma_fit$arl_exceed)[1]), # RL if the control limit is fixed to ensure a 200 IC ARL. This value should be 200.
       far_RL = ifelse(far_FAR == 0, 2001, which(mewma_fit$far_exceed)[1])) # RL if the control limit is fixed to ensure a 0.005 FAR
}
sim_results <- map_dfr(1:1000, sim)

apply(sim_results, 2, mean) # average AR(1) coefficient, ARL's, and FAR's
apply(sim_results, 2, sd) # uncertainty
hist(sim_results$phi, breaks = 50) # monitoring statistics are strongly autocorrelated
hist(sim_results$arl_FAR, breaks = 50) # FAR is much larger than 0.005
hist(sim_results$far_FAR, breaks = 50) # FAR is close to 0.005
hist(sim_results$arl_RL, breaks = 50) # the RL's taper off as expected
hist(sim_results$far_RL, breaks = 50) # many very large RL values
