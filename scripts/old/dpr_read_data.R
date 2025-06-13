library("tidyverse")
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())
library("lubridate")
library("here")
library("forecast")
library("MTS")
library("magrittr")
source(here("scripts", "pl_fd_helper.R"))
source(here("scripts", "pl_fd_methods.R"))


# Read data ---------------------------------------------------------------
df1 <- readRDS(here("data_dpr", "dpr_fault_exp_1.rds"))


# Model monitored variables -----------------------------------------------
vars_monitored <- c(#"O3-Ozone Concentration (ppm)",
                    "UF-Feed Pressure (psi)",
                    #"UF-Permeate Pressure (psi)",
                    "UF-Feed Flow (GPM)",
                    #"BAF-Turbidity [2] (NTU)",
                    "UF-Turbidity [3] (NTU)")

df1_monitored <- df1[, vars_monitored] |> filter(`UF-Feed Pressure (psi)` > 10)
X1 <- as.matrix(df1_monitored)

X1_fit1 <- VAR(X1)
X1_fit2 <- VARMA(X1, 1, 1)

df1_residuals1 <- X1_fit1$residuals |> 
  set_colnames(vars_monitored) |> 
  as_tibble() 

df1_residuals2 <- X1_fit2$residuals |> 
  set_colnames(vars_monitored) |> 
  as_tibble()


# Comparison --------------------------------------------------------------
# RMSE
X1_fit1$residuals^2 |> colMeans() |> sqrt()
X1_fit2$residuals^2 |> colMeans() |> sqrt()

# TS Plots ----------------------------------------------------------------
df1 |> 
  filter(`UF-Feed Pressure (psi)` > 10) |> 
  pivot_longer(all_of(vars_monitored), names_transform = as_factor) |> 
  ggplot(aes(Date_Time, value)) +
  geom_line() +
  facet_wrap(~ name) +
  labs(title = "Experiment 1: Ozone Shutoff",
       x = "Date", y = "", color = "")

df1 |> 
  mutate(`UF-Feed Pressure (psi)` = runmed(`UF-Feed Pressure (psi)`, 51)) |> 
  pivot_longer(all_of(vars_monitored), names_transform = as_factor) |> 
  ggplot(aes(Date_Time, value)) +
  geom_line() +
  facet_wrap(~ name) +
  labs(title = "Experiment 1: Ozone Shutoff",
       x = "Date", y = "", color = "")

df1 |> 
  pivot_longer(c(`UF-Air Scour Valve State`, `UF-Backwash Pump State`)) |> 
  ggplot(aes(Date_Time, value, color = name)) +
  geom_point(shape = 21) +
  geom_line()

df1 |> 
  filter(`UF-Backwash Pump State` == "ON") |> 
  pivot_longer(all_of(vars_monitored), names_transform = as_factor) |> 
  ggplot(aes(Date_Time, value)) +
  geom_point(shape = 21) +
  facet_wrap(~ name) +
  labs(title = "Experiment 1: Ozone Shutoff",
       x = "Date", y = "", color = "")

df1 |> 
  filter(`UF-Backwash Pump State` == "OFF") |> 
  pivot_longer(all_of(vars_monitored), names_transform = as_factor) |> 
  ggplot(aes(Date_Time, value)) +
  geom_point(shape = 21) +
  geom_smooth() +
  facet_wrap(~ name) +
  labs(title = "Experiment 1: Ozone Shutoff",
       x = "Date", y = "", color = "")

df1_residuals1 |>
  mutate(index = 1:nrow(df1_residuals1)) |> 
  pivot_longer(-index, names_transform = as_factor) |> 
  ggplot(aes(index, value)) +
  geom_point(shape = 21) +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Experiment 1: Ozone Shutoff (VAR(1) Residuals)",
       x = "Date", y = "", color = "")

df1_residuals2|> 
  mutate(Date_Time = df1$Date_Time[-1]) |> 
  pivot_longer(-Date_Time, names_transform = as_factor) |> 
  ggplot(aes(Date_Time, value)) +
  geom_point(shape = 21) +
  facet_wrap(~ name, scales = "free") +
  labs(title = "Experiment 1: Ozone Shutoff (VARMA(1, 1) Residuals)",
       x = "Date", y = "", color = "")

# QQ Plots ----------------------------------------------------------------
df1 |> 
  pivot_longer(vars_monitored, names_transform = as_factor) |> 
  ggplot(aes(sample = value)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~ name, scales = "free")

df1_residuals1 |> 
  pivot_longer(vars_monitored, names_transform = as_factor) |> 
  ggplot(aes(sample = value)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~ name, scales = "free")

mu_0 <- colMeans(df1_residuals1)
sigma_0 <- cov(df1_residuals1[1:5000, ])
A <- get_A(sigma_0)
U <- calc_U(df1_residuals1, mu_0, A)

df1_U1 <- U |> set_colnames(vars_monitored) |> as_tibble()

df1_hawk      <- hawkins_2008(df1_residuals1, mu_0, sigma_0)
df1_MC_COMET  <- MC_COMET(df1_residuals1, mu_0, sigma_0, beta = 0.1, cutoff = 0.1, n_w = 180)
df1_MAC_COMET <- MAC_COMET(df1_residuals1, mu_0, sigma_0, beta = 0.1, cutoff = 0.1, n_w = 180)

df1_hawk |> 
  mutate(index = 1:nrow(df1_hawk)) |> 
  ggplot(aes(index, pstat)) +
  geom_line() +
  geom_point(shape = 21) +
  labs(title = "Experiment 1: Hawkins on VAR(1) Residuals")

df1_MC_COMET |> 
  mutate(index = 1:nrow(df1_MC_COMET)) |> 
  ggplot(aes(index, pstat)) +
  geom_line() +
  geom_point(shape = 21) +
  labs(title = "Experiment 1: MC-COMET on VAR(1) Residuals")

df1_MAC_COMET |> 
  mutate(index = 1:nrow(df1_MAC_COMET)) |> 
  pivot_longer(c(T_1, T_2), names_transform = as_factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  labs(title = "Experiment 1: MAC-COMET on VAR(1) Residuals")


df1_U1[1:17280, ] |> 
  pivot_longer(vars_monitored, names_transform = as_factor) |> 
  ggplot(aes(sample = value)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~ name, scales = "free")

df1_residuals2 |> 
  pivot_longer(vars_monitored, names_transform = as_factor) |> 
  ggplot(aes(sample = value)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~ name, scales = "free")


# ACF plots ---------------------------------------------------------------
ggAcf(df1_monitored)
ggAcf(df1_residuals1 |> select(-Date_Time))
ggAcf(df1_residuals2 |> select(-Date_Time))

