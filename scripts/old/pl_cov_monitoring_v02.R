library("tidyverse")
theme_set(cowplot::theme_cowplot())
theme_update(plot.title = element_text(hjust = 0.5, size = 20),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             axis.title.x = element_text(size = 15),
             axis.title.y = element_text(size = 15),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())
library("mvtnorm")
library("here")
source(here("scripts", "pl_fd_methods_v01.R"))


# Wang, Yeh, Li, 2014, Monitoring Mean and Covariance ---------------------
n      <- 1000
p      <- 3
mat_f4 <- matrix(c(1, 0.9, 0, 0.9, 1, 0, 0, 0, 1), byrow = TRUE, nrow = 3)

# Generate sample data (IC, var shift, mean shift, mean + var shift)
set.seed(321)
x_ic <- rmvnorm(n, rep(0, p), diag(p))   |> as_tibble()
x_f1 <- rmvnorm(n, rep(0, p), 2*diag(p)) |> as_tibble()
x_f2 <- rmvnorm(n, rep(1, p), diag(p))   |> as_tibble()
x_f3 <- rmvnorm(n, rep(1, p), 2*diag(p)) |> as_tibble()
x_f4 <- rmvnorm(n, rep(0, p), mat_f4)    |> as_tibble()

mu_0 <- colMeans(x_ic)
sigma_0   <- cov(x_ic)

# Apply methods
pstat_wang_1_f1 <- wang_2014_1(bind_rows(x_ic, x_f1), mu_0, sigma_0)
pstat_wang_1_f2 <- wang_2014_1(bind_rows(x_ic, x_f2), mu_0, sigma_0)
pstat_wang_1_f3 <- wang_2014_1(bind_rows(x_ic, x_f3), mu_0, sigma_0)
pstat_wang_1_f4 <- wang_2014_1(bind_rows(x_ic, x_f4), mu_0, sigma_0)

pstat_wang_2_f1 <- wang_2014_2(bind_rows(x_ic, x_f1), mu_0, sigma_0)
pstat_wang_2_f2 <- wang_2014_2(bind_rows(x_ic, x_f2), mu_0, sigma_0)
pstat_wang_2_f3 <- wang_2014_2(bind_rows(x_ic, x_f3), mu_0, sigma_0)
pstat_wang_2_f4 <- wang_2014_2(bind_rows(x_ic, x_f4), mu_0, sigma_0)

pstat_hawkins_f1 <- hawkins_2008(bind_rows(x_ic, x_f1), mu_0, sigma_0)
pstat_hawkins_f2 <- hawkins_2008(bind_rows(x_ic, x_f2), mu_0, sigma_0)
pstat_hawkins_f3 <- hawkins_2008(bind_rows(x_ic, x_f3), mu_0, sigma_0)
pstat_hawkins_f4 <- hawkins_2008(bind_rows(x_ic, x_f4), mu_0, sigma_0)

pstat_hawkins_sparse_f1 <- hawkins_2008(bind_rows(x_ic, x_f1), mu_0, sigma_0, lambda_s = .2, method = "spcov")
pstat_hawkins_sparse_f2 <- hawkins_2008(bind_rows(x_ic, x_f2), mu_0, sigma_0, lambda_s = .2, method = "spcov")
pstat_hawkins_sparse_f3 <- hawkins_2008(bind_rows(x_ic, x_f3), mu_0, sigma_0, lambda_s = .2, method = "spcov")
pstat_hawkins_sparse_f4 <- hawkins_2008(bind_rows(x_ic, x_f4), mu_0, sigma_0, lambda_s = .2, method = "spcov")

pstat_hawkins_comet_f1 <- hawkins_2008(bind_rows(x_ic, x_f1), mu_0, sigma_0, lambda_s = .2, method = "comet")
pstat_hawkins_comet_f2 <- hawkins_2008(bind_rows(x_ic, x_f2), mu_0, sigma_0, lambda_s = .2, method = "comet")
pstat_hawkins_comet_f3 <- hawkins_2008(bind_rows(x_ic, x_f3), mu_0, sigma_0, lambda_s = .2, method = "comet")
pstat_hawkins_comet_f4 <- hawkins_2008(bind_rows(x_ic, x_f4), mu_0, sigma_0, lambda_s = .2, method = "comet")

# Get h
h_wang_1  <- estimate_h(pstat_wang_1_f1$pstat[1:1000], 200)
h1_wang_2 <- estimate_h(pstat_wang_2_f1$T_1[1:1000], 200)
h2_wang_2 <- estimate_h(pstat_wang_2_f1$T_2[1:1000], 200)
h_hawkins <- estimate_h(pstat_hawkins_f1$pstat[1:1000], 200)
h_hawkins_sparse <- estimate_h(pstat_hawkins_sparse_f1$pstat[1:1000], 200)
h_hawkins_comet <- estimate_h(pstat_hawkins_comet_f1$pstat[1:1000], 200)

df_wang_1 <-
  bind_rows(list(f1 = pstat_wang_1_f1, 
                 f2 = pstat_wang_1_f2, 
                 f3 = pstat_wang_1_f3,
                 f4 = pstat_wang_1_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

df_wang_2 <-
  bind_rows(list(f1 = pstat_wang_2_f1, 
                 f2 = pstat_wang_2_f2, 
                 f3 = pstat_wang_2_f3,
                 f4 = pstat_wang_2_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

df_hawkins <-
  bind_rows(list(f1 = pstat_hawkins_f1, 
                 f2 = pstat_hawkins_f2, 
                 f3 = pstat_hawkins_f3,
                 f4 = pstat_hawkins_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

df_hawkins_sparse <-
  bind_rows(list(f1 = pstat_hawkins_sparse_f1, 
                 f2 = pstat_hawkins_sparse_f2, 
                 f3 = pstat_hawkins_sparse_f3,
                 f4 = pstat_hawkins_sparse_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

df_hawkins_comet <-
  bind_rows(list(f1 = pstat_hawkins_comet_f1, 
                 f2 = pstat_hawkins_comet_f2, 
                 f3 = pstat_hawkins_comet_f3,
                 f4 = pstat_hawkins_comet_f4), .id = "fault") |> 
  mutate(fault = as_factor(fault),
         index = 1:n(), .by = "fault")

# Plots
width = .5*8000
height = .5*5000

# Plots -------------------------------------------------------------------

fault_labels <- c(f1 = "Fault 1",
                  f2 = "Fault 2",
                  f3 = "Fault 3",
                  f4 = "Fault 4") |> as_labeller()

# Wang Chart 1 Plots
df_wang_1 |> 
  ggplot(aes(index, pstat)) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  geom_hline(yintercept = h_wang_1, color = "darkgreen", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x", labeller = fault_labels) +
  labs(y = "Plotting Statistic",
       title = "Wang Chart 1: Sample Plot of PGLR_t")

ggsave(here("figures", "wang_1_pstat.png"), width = width, height = height, units = 'px')

df_wang_1 |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1, labeller = fault_labels) + 
  labs(y = "Components of the PGLR_t",
       title = "Wang Chart 1: Plotting Statistic Breakdown",
       color = "")

ggsave(here("figures", "wang_1_breakdown.png"), width = width, height = height, units = 'px')

# Chart 2 Plots
df_wang_2 |> 
  pivot_longer(contains("T_"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  geom_hline(yintercept = h1_wang_2, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = h2_wang_2, color = "purple",    linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x", labeller = fault_labels) +
  labs(y = "Plotting Statistic",
       title = "Wang Chart 2: Sample Plot of T_1 and T_2", color = "")

ggsave(here("figures", "wang_2_pstat.png"), width = width, height = height, units = 'px')

df_wang_2 |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1, labeller = fault_labels) + 
  labs(y = "Components of T_2",
       title = "Wang Chart 2: Plotting Statistic Breakdown",
       color = "")

ggsave(here("figures", "wang_2_breakdown.png"), width = width, height = height, units = 'px')

# Hawkins Plot
df_hawkins |> 
  ggplot(aes(index, pstat)) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  geom_hline(yintercept = h_hawkins, color = "darkgreen", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x", labeller = fault_labels) +
  labs(y = "Plotting Statistic",
       title = "Hawkins Chart: Sample Plot of c_t")

ggsave(here("figures", "hawkins_pstat.png"), width = width, height = height, units = 'px')

df_hawkins |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1, labeller = fault_labels) + 
  labs(y = "Components of Hawkins Plotting Statistic",
       title = "Hawkins Chart: Plotting Statistic Breakdown",
       color = "")

ggsave(here("figures", "hawkins_breakdown.png"), width = width, height = height, units = 'px')

# Hawkins Sparse Plot
df_hawkins_sparse |> 
  ggplot(aes(index, pstat)) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  geom_hline(yintercept = h_hawkins_sparse, color = "darkgreen", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x", labeller = fault_labels) +
  labs(y = "Plotting Statistic",
       title = "Hawkins Sparse Chart: Sample Plot of c_t with sparsified S")

ggsave(here("figures", "hawkins_sparse_pstat.png"), width = width, height = height, units = 'px')

df_hawkins_sparse |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1, labeller = fault_labels) + 
  labs(y = "Components of Hawkins Sparse Plotting Statistic",
       title = "Hawkins Sparse Chart: Plotting Statistic Breakdown",
       color = "")

ggsave(here("figures", "hawkins_sparse_breakdown.png"), width = width, height = height, units = 'px')

# Hawkins Comet Plot
df_hawkins_comet |> 
  ggplot(aes(index, pstat)) +
  geom_line(alpha = .5) +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  geom_hline(yintercept = h_hawkins_comet, color = "darkgreen", linetype = "dashed") +
  facet_wrap(~ fault, ncol = 1, scales = "free_x", labeller = fault_labels) +
  labs(y = "Plotting Statistic",
       title = "Hawkins Comet Chart: Sample Plot of c_t with moving window S")

ggsave(here("figures", "hawkins_comet_pstat.png"), width = width, height = height, units = 'px')

df_hawkins_comet |> 
  pivot_longer(contains("term"), names_transform = factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_vline(xintercept = 1000, color = "red", linetype = "dashed") +
  facet_wrap(~ fault, scales = "free_x", ncol = 1, labeller = fault_labels) + 
  labs(y = "Components of Hawkins Comet Plotting Statistic",
       title = "Hawkins Comet Chart: Plotting Statistic Breakdown",
       color = "")

ggsave(here("figures", "hawkins_comet_breakdown.png"), width = width, height = height, units = 'px')
