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
    labs(title = title, color = "", x = "", y = "") +
    scale_color_viridis_d(begin = 0.05, end = .95, direction = -1)
}


# Data/method specifications ----------------------------------------------
source(here("scripts", "sim_study_settings.R"))
param_id <- 4

# Generate data -----------------------------------------------------------
set.seed(123)
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
df_s1 |> plot_scenario(title = "Scenario 1: Shift in Mean")       -> plot_1
df_s2 |> plot_scenario(title = "Scenario 2: Shift in Variance")   -> plot_2
df_s3 |> plot_scenario(title = "Scenario 4: Shift in Covariance") -> plot_3
df_s4 |> plot_scenario(title = "Scenario 4")
df_s5 |> plot_scenario(title = "Scenario 5")
df_s6 |> plot_scenario(title = "Scenario 6")
df_s7 |> plot_scenario(title = "Scenario 7")
df_s8 |> plot_scenario(title = "Scenario 8")
df_s9 |> plot_scenario(title = "Scenario 9")
df_s10|> plot_scenario(title = "Scenario 10")


# TS plot for paper -------------------------------------------------------
plot_1 /
  plot_2 /
  plot_3 +
  labs(x = "Observation Index") +
  plot_layout(guides = "collect")


ggsave("figures/ts_plot_sim.pdf", width = 4000, height = 2000, units = "px")