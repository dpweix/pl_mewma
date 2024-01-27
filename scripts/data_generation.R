library("tidyverse")
library("magrittr")
library("mvtnorm")
library("sn")
library("here")

# Data generation helper functions ----------------------------------------
format_matrix <- function(X) {
  X |> 
    set_colnames(paste0("x", 1:ncol(X))) |> 
    as_tibble()
}

sigma_1 <- function(rho = .2) {
  matrix(c(1, rho, 0,
           rho, 1, 0,
           0, 0, 1), ncol = 3, byrow = TRUE)
}

mu_1 <- function(a = 1) {
  c(rep(a, 3), rep(0, 7))
}

sigma_2 <- function(b = 1, rho = .2) {
  matrix(c(b, b*rho, b*rho, b*rho, b*rho, 0, 0, 0, 0, 0,
           b*rho, b, b*rho, b*rho, b*rho, 0, 0, 0, 0, 0,
           b*rho, b*rho, b, b*rho, b*rho, 0, 0, 0, 0, 0,
           b*rho, b*rho, b*rho, b, b*rho, 0, 0, 0, 0, 0,
           b*rho, b*rho, b*rho, b*rho, b, 0, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ncol = 10, byrow = TRUE)
}

# sigma_2 <- function(b = 1, rho = .2) {
#   matrix(c(b, b*rho, b*rho, b*rho, 0, 0, 0, 0, 0, 0,
#            b*rho, b, b*rho, b*rho, 0, 0, 0, 0, 0, 0,
#            b*rho, b*rho, b, b*rho, 0, 0, 0, 0, 0, 0,
#            b*rho, b*rho, b*rho, b, 0, 0, 0, 0, 0, 0,
#            0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
#            0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
#            0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
#            0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
#            0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
#            0, 0, 0, 0, 0, 0, 0, 0, 0, 1), ncol = 10, byrow = TRUE)
# }

sigma_3 <- function(rho = .2) {
  matrix(c(1    , rho  , rho^2, rho^3, rho^4, rho^5, rho^6, rho^7, rho^8, rho^9,
           rho  , 1    , rho  , rho^2, rho^3, rho^4, rho^5, rho^6, rho^7, rho^8,
           rho^2, rho  , 1    , rho  , rho^2, rho^3, rho^4, rho^5, rho^6, rho^7,
           rho^3, rho^2, rho  , 1    , rho  , rho^2, rho^3, rho^4, rho^5, rho^6,
           rho^4, rho^3, rho^2, rho  , 1    , rho  , rho^2, rho^3, rho^4, rho^5,
           rho^5, rho^4, rho^3, rho^2, rho  , 1    , rho  , rho^2, rho^3, rho^4,
           rho^6, rho^5, rho^4, rho^3, rho^2, rho  , 1    , rho  , rho^2, rho^3,
           rho^7, rho^6, rho^5, rho^4, rho^3, rho^2, rho  , 1    , rho  , rho^2,
           rho^8, rho^7, rho^6, rho^5, rho^4, rho^3, rho^2, rho  , 1    , rho  ,
           rho^9, rho^8, rho^7, rho^6, rho^5, rho^4, rho^3, rho^2, rho  , 1), ncol = 10, byrow = TRUE)
}

rho_0 <- .2

# Data generation functions -----------------------------------------------

gen_dat_s1 <- function(n, a, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean =   rep(0, 3), sigma = diag(3)) |> format_matrix(),
              OC = rmvnorm(n, mean = a*rep(1, 3), sigma = diag(3)) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:5]
  } else {
    df
  }
}

gen_dat_s2 <- function(n, b, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 3), sigma =   diag(3)) |> format_matrix(),
              OC = rmvnorm(n, mean = rep(0, 3), sigma = b*diag(3)) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:5]
  } else {
    df
  }
}

gen_dat_s3 <- function(n, rho, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 3), sigma = sigma_1(rho_0)) |> format_matrix(),
              OC = rmvnorm(n, mean = rep(0, 3), sigma = sigma_1(rho  )) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:5]
  } else {
    df
  }
}

gen_dat_s4 <- function(n, b, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 3), sigma =   sigma_1(rho_0)) |> format_matrix(),
              OC = rmvnorm(n, mean = rep(0, 3), sigma = b*sigma_1(rho_0)) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:5]
  } else {
    df
  }
}

gen_dat_s5 <- function(n, b, n_drift = 100, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 3), sigma =   sigma_1(rho_0)) |> format_matrix(),
              OC = bind_rows(map_dfr(seq(1, b, length.out = n_drift + 1)[-1], \(b_d) {
                    rmvnorm(1          , mean = rep(0, 3), sigma = b_d*sigma_1(rho_0)) |> format_matrix()}),
                    rmvnorm(n - n_drift, mean = rep(0, 3), sigma =   b*sigma_1(rho_0)) |> format_matrix()
                    ),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:5]
  } else {
    df
  }
}

gen_dat_s6 <- function(n, a, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 10), sigma = diag(10)) |> format_matrix(),
              OC = rmvnorm(n, mean =    mu_1(a), sigma = diag(10)) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:12]
  } else {
    df
  }
}

gen_dat_s7 <- function(n, b, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 10), sigma = sigma_2(b = 1, rho = rho_0)) |> format_matrix(),
              OC = rmvnorm(n, mean = rep(0, 10), sigma = sigma_2(b = b, rho = rho_0)) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:12]
  } else {
    df
  }
}

gen_dat_s8 <- function(n, rho, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmvnorm(n, mean = rep(0, 10), sigma = sigma_3(rho = rho_0)) |> format_matrix(),
              OC = rmvnorm(n, mean = rep(0, 10), sigma = sigma_3(rho = rho  )) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:12]
  } else {
    df
  }
}

gen_dat_s9 <- function(n, a, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmst(n, xi = rep(0, 10), Omega = diag(10), alpha = rep(0, 10), nu = 40) |> format_matrix(),
              OC = rmst(n, xi =    mu_1(a), Omega = diag(10), alpha = rep(0, 10), nu = 40) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:12]
  } else {
    df
  }
}

gen_dat_s10 <- function(n, b, data_only = TRUE) {
  df <- 
    bind_rows(IC = rmst(n, xi = rep(0, 10), Omega = sigma_2(b = 1, rho = rho_0), alpha = rep(0, 10), nu = 40) |> format_matrix(),
              OC = rmst(n, xi = rep(0, 10), Omega = sigma_2(b = b, rho = rho_0), alpha = rep(0, 10), nu = 40) |> format_matrix(),
              .id = "type") |> 
      mutate(index = 1:n(), .before = "x1")
  
  if(data_only) {
    df[, 3:12]
  } else {
    df
  }
}

gen_dat <- function(scenario, n, arg, data_only = TRUE) {
  if(scenario == "s1") gen_dat_s1(n, arg, data_only)
  else if(scenario == "s2") gen_dat_s2(n, arg, data_only)
  else if(scenario == "s3") gen_dat_s3(n, arg, data_only)
  else if(scenario == "s4") gen_dat_s4(n, arg, data_only)
  else if(scenario == "s5") gen_dat_s5(n, arg, n_drift = 100, data_only)
  else if(scenario == "s6") gen_dat_s6(n, arg, data_only)
  else if(scenario == "s7") gen_dat_s7(n, arg, data_only)
  else if(scenario == "s8") gen_dat_s8(n, arg, data_only)
  else if(scenario == "s9") gen_dat_s9(n, arg, data_only)
  else if(scenario == "s10") gen_dat_s10(n, arg, data_only)
}
