library("tidyverse")
library("here")
library("kableExtra")

# Plot theme
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())

scenario_labels <- c("s1" = "Scenario 1, a",
                     "s2" = "Scenario 2, b",
                     "s3" = "Scenario 3, ρ",
                     "s4" = "Scenario 4, b",
                     "s5" = "Scenario 5, a",
                     "s6" = "Scenario 6, b",
                     "s7" = "Scenario 7, \u03C1",
                     "s8" = "Scenario 8, \u03C1",
                     "s9" = "Scenario 9, a",
                     "s10" = "Scenario 10, b") |> as_labeller()

# Load pstats -------------------------------------------------------------
source(here("scripts", "sim_study_settings.R"))
source(here("scripts", "pl_fd_methods.R"))

# pstats <- 
#   map(i_min:i_max, \(i) {
#     map(names(scenario_params), \(s) {
#       map(names(method_params), \(m) {
#         file_names <- here(data_folder, paste0(s, "-", m, "-", scenario_params[[s]], "-", i, ".rds"))
#         result_1 <- read_rds(file_names[[1]])
#         result_2 <- read_rds(file_names[[2]])
#         result_3 <- read_rds(file_names[[3]])
#         result_4 <- read_rds(file_names[[4]])
#         
#         list(param1 = result_1$pstat,
#              param2 = result_2$pstat,
#              param3 = result_3$pstat,
#              param4 = result_4$pstat)
#         
#       }) |> set_names(names(method_params))
#     }) |> set_names(names(scenario_params))
#   })

df_dr <- 
  map_dfr(i_min:i_max, \(i) {
    map_dfr(names(scenario_params), \(s) {
      map_dfr(names(method_params), \(m) {
        file_names <- here(data_folder, paste0(s, "-", m, "-", scenario_params[[s]], "-", i, ".rds"))
        result_1 <- read_rds(file_names[[1]])
        result_2 <- read_rds(file_names[[2]])
        result_3 <- read_rds(file_names[[3]])
        result_4 <- read_rds(file_names[[4]])
        
        if(m %in% c("hawkins", "wang_1", "MC_LASSO", "MC_COMET")) {
          h_param1 <- estimate_h_fd(result_1$pstat[1:1000], .005)
          h_param2 <- estimate_h_fd(result_2$pstat[1:1000], .005)
          h_param3 <- estimate_h_fd(result_3$pstat[1:1000], .005)
          h_param4 <- estimate_h_fd(result_4$pstat[1:1000], .005)
          
          tibble(index = i,
                 method = m,
                 scenario = s,
                 param = scenario_params[[s]],
                 dr = c(get_detection_rate(result_1$pstat[1001:2000], h_param1),
                        get_detection_rate(result_2$pstat[1001:2000], h_param2),
                        get_detection_rate(result_3$pstat[1001:2000], h_param3),
                        get_detection_rate(result_4$pstat[1001:2000], h_param4)))
          
        } else if(m %in% c("wang_2", "MAC_COMET", "MAC_COMET1")) {
          h_param1 <- estimate_h2_fd(result_1$pstat$T_1[1:1000], result_1$pstat$T_2[1:1000], .005)
          h_param2 <- estimate_h2_fd(result_2$pstat$T_1[1:1000], result_2$pstat$T_2[1:1000], .005)
          h_param3 <- estimate_h2_fd(result_3$pstat$T_1[1:1000], result_3$pstat$T_2[1:1000], .005)
          h_param4 <- estimate_h2_fd(result_4$pstat$T_1[1:1000], result_4$pstat$T_2[1:1000], .005)
          
          tibble(index = i,
                 method = m,
                 scenario = s,
                 param = scenario_params[[s]],
                 dr = c(get_detection_rate2(result_1$pstat$T_1[1001:2000],
                                            result_1$pstat$T_2[1001:2000],
                                            h_param1[1], h_param1[2]),
                        get_detection_rate2(result_2$pstat$T_1[1001:2000],
                                            result_2$pstat$T_2[1001:2000],
                                            h_param2[1], h_param2[2]),
                        get_detection_rate2(result_3$pstat$T_1[1001:2000],
                                            result_3$pstat$T_2[1001:2000],
                                            h_param3[1], h_param3[2]),
                        get_detection_rate2(result_4$pstat$T_1[1001:2000],
                                            result_4$pstat$T_2[1001:2000],
                                            h_param4[1], h_param4[2])))
        }
        
      })
    })
  })

df_dr

get_ic_param <- function(scenario = "s1") {
  if(scenario %in% c("s1", "s5", "s9")) {
    return(0)
  } else if(scenario %in% c("s2", "s4", "s6", "s10")) {
    return(1)
  } else if(scenario %in% c("s3", "s7", "s8")) {
    return(0.2)
  }
  NA
}

# Estimate average run lengths
df_adr <- df_dr |> 
  mutate(ic_param = map_dbl(scenario, get_ic_param)) |> 
  mutate(scenario = as_factor(scenario),
         method = case_match(method,
                             "MC_LASSO" ~ "MC-LASSO", 
                             "MC_COMET" ~ "MC-COMET",
                             "MAC_COMET" ~ "MAC-COMET",
                             "hawkins" ~ "MC-Hawkins", 
                             "wang_1" ~ "MC-Wang",
                             "wang_2" ~ "MAC-Wang",
                             .default = method) |> as_factor()) |>
  group_by(scenario, method, param) |> 
  summarize(`Mean Detection Rate` = round(mean(dr), 3),
            `Median Detection Rate` = median(dr),
            `SE of Detection Rate` = sd(dr)/sqrt(i_max),
            ic_param = mean(ic_param)) |>
  ungroup()


# Generate plots ----------------------------------------------------------

# Detection rate plot
df_adr |> 
  group_by(scenario) |> 
  ggplot(aes(param, `Mean Detection Rate`, color = method, linetype = method)) +
  geom_point(aes(x = ic_param, y = 0), shape = 17, size = 3, color = "grey30") +
  geom_segment(aes(x = ic_param, xend = ic_param, y = 0, yend = 1), 
               color = "grey30", linetype = "dashed") +
  geom_point(shape = 21) +
  geom_line() +
  facet_wrap(~ scenario, scales = "free", nrow = 2, ncol = 5, labeller = scenario_labels) +
  scale_x_continuous(n.breaks = 6) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(color = "", linetype = "", x = "", title = "Mean Detection Rate by Scenario, Method, and Parameter")

ggsave(here("figures", "dr-results.pdf"), width = 4000, height = 2000, units = "px", device = cairo_pdf)


# Detection rate table
df_adr |> 
  mutate(value = paste0(`Mean Detection Rate`, " (", round(`SE of Detection Rate`, 1) ,")")) |> 
  select(scenario, method, param, value) |> 
  pivot_wider(names_from = method, values_from = value) |> 
  kbl(format = "latex", booktabs = TRUE, linesep =  c('', '', '', '\\addlinespace'))


# RL Histograms
scenario_n <- "s3"

df_dr |> 
  mutate(method = as_factor(method)) |> 
  filter(scenario == scenario_n) |> 
  ggplot(aes(dr)) +
  geom_histogram() +
  facet_wrap(method ~ param, scale = "free", ncol = 4) +
  labs(title = paste0("RL Histogram for Scenario ", str_remove(scenario_n, "s")))

ggsave(here("figures", paste0(scenario_n, "-histogram.pdf")), width = 4000, height = 4000, units = "px", device = cairo_pdf)

# Histogram plot theme
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside")

method <- "hawkins"

df_dr |> 
  filter(method == method) |> 
  ggplot(aes(dr)) +
  geom_histogram() +
  facet_wrap(scenario ~ param, ncol = 5, scales = "free") +
  labs(x = "Run Length", y = "Count", title = paste0("Detection Rate Histograms: ", method))

ggsave(here("figures", paste0(method, "-dr-histogram.pdf")), width = 4000, height = 5500, units = "px")

# Generate tables ---------------------------------------------------------