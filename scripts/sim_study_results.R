library("tidyverse")
library("here")
library("kableExtra")

# ARL plot theme
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
                     "s5" = "Scenario 5, b",
                     "s6" = "Scenario 6, a",
                     "s7" = "Scenario 7, b",
                     "s8" = "Scenario 8, \u03C1",
                     "s9" = "Scenario 9, a",
                     "s10" = "Scenario 10, b") |> as_labeller()

# Load results ------------------------------------------------------------
source(here("scripts", "sim_study_settings.R"))

# Collect the run lengths
df_rl <- 
  map_dfr(names(scenario_params), \(s) {
    map_dfr(names(method_params), \(m) {
      map_dfr(scenario_params[s], \(param) {
        map_dfr(param, \(value) {
          map_dfr(i_min:i_max, \(i) {
            tibble(scenario = s,
                   method = m,
                   param = value,
                   rl = readRDS(
                     here(data_folder, paste0(s, "-", m, "-", value, "-", i, ".rds")))$rl
                   )
          })
        })
      })
    })
  }) ##|> 
  #mutate(scenario = as_factor(scenario),
  #       method = as_factor(method))

# Trim the outlier run lengths
df_rl_trim <- df_rl |> 
  group_by(scenario, method, param) |> 
  filter(between(rl, quantile(rl, .05), quantile(rl, .95)))

# Estimate average run lengths
df_arl <- df_rl |> 
  mutate(scenario = as_factor(scenario),
         method = case_match(method,
                             "hawkins" ~ "Hawkins", 
                             "wang_1" ~ "Wang_1",
                             "wang_2" ~ "Wang_2",
                             .default = method) |> as_factor()) |>
  group_by(scenario, method, param) |> 
  summarize(ARL = mean(rl),
            MRL = median(rl)) |>
  ungroup() |> 
  mutate(param_min = min(param), 
         arl_max = max(ARL), .by = scenario)
  

# Generate plots ----------------------------------------------------------

# ARL plot
df_arl |> 
  group_by(scenario) |> 
  #filter(method %in% c("MAC_COMET")) |> 
  ggplot(aes(param, ARL, color = method, linetype = method)) +
  geom_point(shape = 21) +
  #geom_point(aes(x = param_min, y = 20), shape = 17, size = 3, color = "black") +
  geom_line() +
  facet_wrap(~ scenario, scales = "free", nrow = 2, ncol = 5, labeller = scenario_labels) +
  scale_x_continuous(n.breaks = 6) +
  #scale_y_continuous(limits = c(0, df_arl$arl_max[1]), breaks = c(0, 60, 120, 180, 240, 300)) +
  labs(color = "", linetype = "", x = "", title = "ARL Values by Scenario, Method, and Parameter")

ggsave(here("figures", "arl-results.pdf"), width = 4000, height = 2000, units = "px", device = cairo_pdf)
ggsave(here("figures", "arl-results100.pdf"), width = 4000, height = 2000, units = "px", device = cairo_pdf)


# ARL Table
df_arl |> 
  select(scenario, method, param, ARL) |> 
  pivot_wider(names_from = method, values_from = ARL) |> 
  kbl(format = "latex", booktabs = TRUE, linesep =  c('', '', '', '\\addlinespace'))


# RL Histograms
scenario_n <- "s10"

df_rl |> 
  mutate(method = as_factor(method)) |> 
  filter(scenario == scenario_n) |> 
  ggplot(aes(rl)) +
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

df_rl |> 
  filter(method == method) |> 
  ggplot(aes(rl)) +
  geom_histogram() +
  facet_wrap(scenario ~ param, ncol = 5, scales = "free") +
  labs(x = "Run Length", y = "Count", title = paste0("ARL Histograms: ", method))

ggsave(here("figures", paste0(method, "-rl-histogram.pdf")), width = 4000, height = 5500, units = "px")

# Generate tables ---------------------------------------------------------