library("tidyverse")
library("here")


# Load results ------------------------------------------------------------
source(here("scripts", "sim_study_settings.R"))
source(here("scripts", "pl_fd_methods.R"))


# View plotting statistics ------------------------------------------------
pstats <- 
  map(i_min:i_max, \(i) {
    map(names(scenario_params), \(s) {
      map(methods, \(m) {
        file_names <- here(data_folder, paste0(s, "-", m, "-", scenario_params[[s]], "-", i, ".rds"))
        result_1 <- read_rds(file_names[[1]])
        result_2 <- read_rds(file_names[[2]])
        result_3 <- read_rds(file_names[[3]])
        result_4 <- read_rds(file_names[[4]])
        
        list(param1 = result_1$pstat,
             param2 = result_2$pstat,
             param3 = result_3$pstat,
             param4 = result_4$pstat)
        
      }) |> set_names(methods)
    }) |> set_names(names(scenario_params))
  })


scenario <- "s1"
method <- "hawkins"
i <- 16
tibble(pstat = pstats[[i]][[scenario]][[method]]$param2,
       h = estimate_h(pstats[[i]][[scenario]][[method]]$param1[1:1001], 180),
       index = 1:2000) |> 
  ggplot(aes(index, pstat)) +
  geom_line() +
  geom_point(shape = 21) +
  geom_hline(aes(yintercept = h), color = 'blue') +
  geom_vline(xintercept = 1000, color = 'red')

acf(pstats[[i]][[scenario]][[method]]$param1)


pstats[[i]][[scenario]][[method]]$param1 |> 
  mutate(index = 1:2000) |> 
  pivot_longer(contains("T_")) |> 
  ggplot(aes(index, value, color = name)) +
  geom_line() +
  geom_point(shape = 21)
    

# Get RL with new control limit -------------------------------------------
df_rl_1 <- 
  map_dfr(i_min:i_max, \(i) {
    map_dfr(names(scenario_params), \(s) {
      map_dfr(names(method_params), \(m) {
        file_names <- here(data_folder, paste0(s, "-", m, "-", scenario_params[[s]], "-", i, ".rds"))
        result_1 <- read_rds(file_names[[1]])
        result_2 <- read_rds(file_names[[2]])
        result_3 <- read_rds(file_names[[3]])
        result_4 <- read_rds(file_names[[4]])
        
        if(m %in% c("hawkins", "wang_1", "MC_LASSO", "MC_COMET")) {
          if(s %in% c("s3", "s8")) {
            h_new <- estimate_h(result_2$pstat[1001:2000], 180)#quantile(result_2$pstat[1001:2000], 1-1/180)
          } else {
            h_new <- estimate_h(result_1$pstat[1001:2000], 180)
          }
          
          tibble(index = i,
                 method = m,
                 scenario = s,
                 param = scenario_params[[s]],
                 rl = c(get_rl_oc(result_1$pstat[1001:2000], h_new),
                        get_rl_oc(result_2$pstat[1001:2000], h_new),
                        get_rl_oc(result_3$pstat[1001:2000], h_new),
                        get_rl_oc(result_4$pstat[1001:2000], h_new)))
          
        } else if(m %in% c("wang_2", "MAC_COMET")) {
          if(s %in% c("s3", "s8")) {
            h1_new <- estimate_h(result_2$pstat$T_1[1001:2000], 180)#quantile(result_2$pstat$T_1[1001:2000], 1-1/180)
            h2_new <- estimate_h(result_2$pstat$T_2[1001:2000], 180)
          } else {
            h1_new <- estimate_h(result_1$pstat$T_1[1001:2000], 180)
            h2_new <- estimate_h(result_1$pstat$T_2[1001:2000], 180)
          }
          
          tibble(index = i,
                 method = m,
                 scenario = s,
                 param = scenario_params[[s]],
                 rl = c(min(get_rl_oc(result_1$pstat$T_1[1001:2000], h1_new),
                            get_rl_oc(result_1$pstat$T_2[1001:2000], h2_new)),
                        min(get_rl_oc(result_2$pstat$T_1[1001:2000], h1_new),
                            get_rl_oc(result_2$pstat$T_2[1001:2000], h2_new)),
                        min(get_rl_oc(result_3$pstat$T_1[1001:2000], h1_new),
                            get_rl_oc(result_3$pstat$T_2[1001:2000], h2_new)),
                        min(get_rl_oc(result_4$pstat$T_1[1001:2000], h1_new),
                            get_rl_oc(result_4$pstat$T_2[1001:2000], h2_new))))
        }
        
      })
    })
  })

df_arl_1 <- df_rl_1 |> 
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

df_arl_1 |> 
  group_by(scenario) |> 
  ggplot(aes(param, MRL, color = method, linetype = method)) +
  geom_point(shape = 21) +
  #geom_point(aes(x = param_min, y = 20), shape = 17, size = 3, color = "black") +
  geom_line() +
  facet_wrap(~ scenario, scales = "free", nrow = 2, ncol = 5, labeller = scenario_labels) +
  scale_x_continuous(n.breaks = 6) +
  #scale_y_continuous(limits = c(0, df_arl_1$arl_max[1]), breaks = c(0, 60, 120, 180, 240, 300)) +
  labs(color = "", linetype = "", x = "", title = "ARL Values by Scenario, Method, and Parameter")

here(data_folder, paste0(s, "-", m, "-", scenario_params[s][1], "-", i, ".rds"))