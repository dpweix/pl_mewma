# File set up -------------------------------------------------------------
library("tidyverse")
theme_set(cowplot::theme_cowplot())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())
library("here")
library("patchwork")
source(here("scripts", "data_generation.R"))
source(here("scripts", "pl_fd_methods.R"))

plot_scenario <- function(df, title = "", type = "color") {
  df <- df |> 
    mutate(index = 1:(2*n),
           type = rep(c("IC", "OC"), each = n) |> as_factor())
  
  if(type == "color") {
    plot <- df |> 
      pivot_longer(matches("x\\d"), names_transform = as_factor) |> 
      ggplot(aes(index, value, color = name)) +
      geom_line() + 
      geom_vline(xintercept = last(which(df$type == "IC")),
                 color = "black", linetype = "dashed")
  } else if(type == "facet") {
    plot <- df |> 
      pivot_longer(matches("x\\d"), names_transform = as_factor) |> 
      ggplot(aes(index, value)) +
      geom_line() + 
      geom_vline(xintercept = last(which(df$type == "IC")),
                 color = "black", linetype = "dashed") +
      facet_wrap(~ name, scales = "free_x")
  }
  
  plot + 
    labs(title = title, color = "", x = "Observation Index", y = "") +
    scale_color_viridis_d(begin = 0, end = .6)
}


# Data/method specifications ----------------------------------------------
source(here("scripts", "sim_study_settings.R"))
param_id <- 3

# Generate data -----------------------------------------------------------
df_s1 <-gen_dat_s1(n,   a[param_id])## |> plot_scenario()
df_s2 <-gen_dat_s2(n,   b[param_id])## |> plot_scenario()
df_s3 <-gen_dat_s3(n, rho[param_id])## |> plot_scenario()
df_s4 <-gen_dat_s4(n,   b[param_id])## |> plot_scenario()
df_s5 <-gen_dat_s5(n,   b[param_id])## |> plot_scenario()
df_s6 <-gen_dat_s6(n,   a[param_id])## |> plot_scenario()
df_s7 <-gen_dat_s7(n,   b[param_id])## |> plot_scenario()
df_s8 <-gen_dat_s8(n, rho[param_id])## |> plot_scenario()
df_s9 <-gen_dat_s9(n,   a[param_id])## |> plot_scenario()
df_s10<-gen_dat_s10(n,  b[param_id])## |> plot_scenario()

# TS Plot of scenarios ----------------------------------------------------
df_s1 |> plot_scenario(title = "Scenario 1")
df_s2 |> plot_scenario(title = "Scenario 2")
df_s3 |> plot_scenario(title = "Scenario 3")
df_s4 |> plot_scenario(title = "Scenario 4")
df_s5 |> plot_scenario(title = "Scenario 5")
df_s6 |> plot_scenario(title = "Scenario 6")
df_s7 |> plot_scenario(title = "Scenario 7")
df_s8 |> plot_scenario(title = "Scenario 8")
df_s9 |> plot_scenario(title = "Scenario 9")
df_s10|> plot_scenario(title = "Scenario 10")

# Scenario 1 --------------------------------------------------------------
hawkins_s1 <- hawkins_2008(df_s1, colMeans(df_s1[1:n, ]), cov(df_s1[1:n, ]))
mc_lasso_s1 <- MC_LASSO(df_s1, colMeans(df_s1[1:n, ]), cov(df_s1[1:n, ]))
mc_comet_s1 <- MC_COMET(df_s1, colMeans(df_s1[1:n, ]), cov(df_s1[1:n, ]))
mac_comet_s1 <- MAC_COMET(df_s1, colMeans(df_s1[1:n, ]), cov(df_s1[1:n, ]))
mac_comet1_s1<- MAC_COMET1(df_s1, colMeans(df_s1[1:n, ]), cov(df_s1[1:n, ]))

h_hawkins <- estimate_h_bootstrap(hawkins_s1$pstat[1:n], arl_ic)
h_mc_lasso <- estimate_h_bootstrap(mc_lasso_s1$pstat[1:n], arl_ic)
h_mc_comet <- estimate_h_bootstrap(mc_comet_s1$pstat[1:n], arl_ic)
h1_mac_comet <- estimate_h_bootstrap(mac_comet_s1$T_1[1:n], arl_ic)
h2_mac_comet <- estimate_h_bootstrap(mac_comet_s1$T_2[1:n], arl_ic)

df_s1_pstat <- df_s1 |> 
  mutate(index = 1:(2*n),
         Hawkins = hawkins_s1$pstat,
         MC_LASSO = mc_lasso_s1$pstat,
         MC_COMET = mc_comet_s1$pstat,
         MAC_COMET_T1 = mac_comet_s1$T_1,
         MAC_COMET_T2 = mac_comet_s1$T_2)


df_s1_pstat |> 
  pivot_longer(-c(x1, x2, x3, index), names_transform = as_factor) |> 
  ggplot(aes(index, value)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ name, scales = "free", nrow = 2)

p_hawkins <- df_s1_pstat |> 
  ggplot(aes(index, Hawkins)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_hawkins, color = "green", linetype = "longdash") +
  labs(title = "Hawkins", x = "", y = "", color = "")

p_mc_lasso <- df_s1_pstat |> 
  ggplot(aes(index, MC_LASSO)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_lasso, color = "green", linetype = "longdash") +
  labs(title = "MC-LASSO", x = "", y = "", color = "")

p_mc_comet <- df_s1_pstat |> 
  ggplot(aes(index, MC_COMET)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_comet, color = "green", linetype = "longdash") +
  labs(title = "MC-COMET", x = "", y = "", color = "")

p_mac_comet <- df_s1_pstat |> 
  pivot_longer(contains("MAC"), names_transform = as_factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h1_mac_comet, color = "firebrick", linetype = "longdash") +
  geom_hline(yintercept = h2_mac_comet, color = "blue", linetype = "longdash") +
  labs(title = "MAC-COMET", x = "", y = "", color = "") + 
  theme(legend.position = c(0.05, .9))

p_hawkins /
  p_mc_lasso /
  p_mc_comet /
  p_mac_comet

ggsave("figures/ts_sample_plots_s1.pdf", width = 3000, height = 3000, unit = "px")

hawkins_s2 <- hawkins_2008(df_s2, colMeans(df_s2[1:n, ]), cov(df_s2[1:n, ]))
mc_lasso_s2 <- MC_LASSO(df_s2, colMeans(df_s2[1:n, ]), cov(df_s2[1:n, ]))
mc_comet_s2 <- MC_COMET(df_s2, colMeans(df_s2[1:n, ]), cov(df_s2[1:n, ]))
mac_comet_s2 <- MAC_COMET(df_s2, colMeans(df_s2[1:n, ]), cov(df_s2[1:n, ]))

h_hawkins <- estimate_h_bootstrap(hawkins_s2$pstat[1:n], arl_ic)
h_mc_lasso <- estimate_h_bootstrap(mc_lasso_s2$pstat[1:n], arl_ic)
h_mc_comet <- estimate_h_bootstrap(mc_comet_s2$pstat[1:n], arl_ic)
h1_mac_comet <- estimate_h_bootstrap(mac_comet_s2$T_1[1:n], arl_ic)
h2_mac_comet <- estimate_h_bootstrap(mac_comet_s2$T_2[1:n], arl_ic)

df_s2_pstat <- df_s2 |> 
  mutate(index = 1:(2*n),
         Hawkins = hawkins_s2$pstat,
         MC_LASSO = mc_lasso_s2$pstat,
         MC_COMET = mc_comet_s2$pstat,
         MAC_COMET_T1 = mac_comet_s2$T_1,
         MAC_COMET_T2 = mac_comet_s2$T_2)


df_s2_pstat |> 
  pivot_longer(-c(x1, x2, x3, index), names_transform = as_factor) |> 
  ggplot(aes(index, value)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ name, scales = "free", nrow = 2)

p_hawkins <- df_s2_pstat |> 
  ggplot(aes(index, Hawkins)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_hawkins, color = "green", linetype = "longdash") +
  labs(title = "Hawkins", x = "", y = "", color = "")

p_mc_lasso <- df_s2_pstat |> 
  ggplot(aes(index, MC_LASSO)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_lasso, color = "green", linetype = "longdash") +
  labs(title = "MC-LASSO", x = "", y = "", color = "")

p_mc_comet <- df_s2_pstat |> 
  ggplot(aes(index, MC_COMET)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_comet, color = "green", linetype = "longdash") +
  labs(title = "MC-COMET", x = "", y = "", color = "")

p_mac_comet <- df_s2_pstat |> 
  pivot_longer(contains("MAC"), names_transform = as_factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h1_mac_comet, color = "firebrick", linetype = "longdash") +
  geom_hline(yintercept = h2_mac_comet, color = "blue", linetype = "longdash") +
  labs(title = "MAC-COMET", x = "", y = "", color = "") + 
  theme(legend.position = c(0.05, .9))

p_hawkins /
  p_mc_lasso /
  p_mc_comet /
  p_mac_comet

ggsave("figures/ts_sample_plots_s2.pdf", width = 3000, height = 3000, unit = "px")


# Scenario 3 --------------------------------------------------------------

hawkins_s3 <- hawkins_2008(df_s3, colMeans(df_s3[1:n, ]), cov(df_s3[1:n, ]))
mc_lasso_s3 <- MC_LASSO(df_s3, colMeans(df_s3[1:n, ]), cov(df_s3[1:n, ]))
mc_comet_s3 <- MC_COMET(df_s3, colMeans(df_s3[1:n, ]), cov(df_s3[1:n, ]))
mac_comet_s3 <- MAC_COMET(df_s3, colMeans(df_s3[1:n, ]), cov(df_s3[1:n, ]))

h_hawkins <- estimate_h_bootstrap(hawkins_s3$pstat[1:n], arl_ic)
h_mc_lasso <- estimate_h_bootstrap(mc_lasso_s3$pstat[1:n], arl_ic)
h_mc_comet <- estimate_h_bootstrap(mc_comet_s3$pstat[1:n], arl_ic)
h1_mac_comet <- estimate_h_bootstrap(mac_comet_s3$T_1[1:n], arl_ic)
h2_mac_comet <- estimate_h_bootstrap(mac_comet_s3$T_2[1:n], arl_ic)

df_s3_pstat <- df_s3 |> 
  mutate(index = 1:(2*n),
         Hawkins = hawkins_s3$pstat,
         MC_LASSO = mc_lasso_s3$pstat,
         MC_COMET = mc_comet_s3$pstat,
         MAC_COMET_T1 = mac_comet_s3$T_1,
         MAC_COMET_T2 = mac_comet_s3$T_2)


df_s3_pstat |> 
  pivot_longer(-c(x1, x2, x3, index), names_transform = as_factor) |> 
  ggplot(aes(index, value)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  facet_wrap(~ name, scales = "free", nrow = 2)

p_hawkins <- df_s3_pstat |> 
  ggplot(aes(index, Hawkins)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_hawkins, color = "green", linetype = "longdash") +
  labs(title = "Hawkins", x = "", y = "", color = "")

p_mc_lasso <- df_s3_pstat |> 
  ggplot(aes(index, MC_LASSO)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_lasso, color = "green", linetype = "longdash") +
  labs(title = "MC-LASSO", x = "", y = "", color = "")

p_mc_comet <- df_s3_pstat |> 
  ggplot(aes(index, MC_COMET)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h_mc_comet, color = "green", linetype = "longdash") +
  labs(title = "MC-COMET", x = "", y = "", color = "")

p_mac_comet <- df_s3_pstat |> 
  pivot_longer(contains("MAC"), names_transform = as_factor) |> 
  ggplot(aes(index, value, color = name)) +
  geom_vline(xintercept = n, color = "red", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = h1_mac_comet, color = "firebrick", linetype = "longdash") +
  geom_hline(yintercept = h2_mac_comet, color = "blue", linetype = "longdash") +
  labs(title = "MAC-COMET", x = "", y = "", color = "") + 
  theme(legend.position = c(0.05, .9))

p_hawkins /
  p_mc_lasso /
  p_mc_comet /
  p_mac_comet

ggsave("figures/ts_sample_plots_s3.pdf", width = 3000, height = 3000, unit = "px")
