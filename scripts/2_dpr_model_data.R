library("tidyverse")
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())
library("glue")
library("lubridate")
library("here")
library("magrittr")
library("patchwork")
source(here("scripts", "pl_fd_helper.R"))
source(here("scripts", "pl_fd_methods.R"))


# Read data ---------------------------------------------------------------
df3 <- readRDS(here("data_dpr", "dpr_fault_exp_3.rds"))

# Identify UF States
uf_states <- c("Manual Control", "Startup",
               "Start BW", "Pre-BW Air Scour", "Backwash", "Post-BW Air Scour",
               "Purge", "End BW", "Run", "other")

# Clean UF State Data
uf_state_vec <- map_chr(df3$`UF State`, \(state) {
  return <- 
    state |> 
    str_extract(uf_states) |> 
    na.omit() |> 
    as.character()
  
  if(length(return) == 1) {
    return(return)
  } else {
    "other"
  }
})

uf_state_vec_reduced <-
  str_replace_all(uf_state_vec, c(#`Start BW` = "AS_pre",
    `Pre-BW Air Scour` = "Pre-BW Air Scour",
    `Post-BW Air Scour` = "Post-BW Air Scour",
    `Backwash` = "Backwash",
    `End BW` = "Purge"))

uf_state_iter <- vector("integer", nrow(df3))
temp <- 0
for(i in 2:length(uf_state_iter)) {
  if(uf_state_vec_reduced[i] == uf_state_vec_reduced[i-1]) {
    uf_state_iter[i] <- temp
    temp <- temp + 1
  } else {
    uf_state_iter[i] <- 0
    temp <- 1
  }
}

time_exp_3 <- c(ymd_hms("2023-06-16 15:00:00 UTC"),
                ymd_hms("2023-06-21 11:00:19 UTC"),
                ymd_hms("2023-06-22 10:00:00 UTC"))

time_exp_3_limits <- c(ymd_hms("2023-06-16 19:00:00 UTC"),
                       ymd_hms("2023-06-22 10:00:00 UTC"))

df3 <- df3 |> 
  mutate(`TMP (psi)` = `UF-Feed Pressure (psi)` - `UF-Permeate Pressure (psi)`,
         `UF State` = as_factor(uf_state_vec_reduced),
         uf_state_iter = uf_state_iter) |> 
  filter(between(Date_Time, time_exp_3[1], time_exp_3[3]))

# Pre, during, and post fault data frames (by fault_occuring variable)
df3_pre_raw <- df3[1:(first(which(df3$fault_occurring == 1))-1), ]
df3_drn_raw <- df3[-c(1:(first(which(df3$fault_occurring == 1))-1)), ]

begin_fault <- df3[first(which(df3$fault_occurring == 1)), ]$Date_Time
end_fault   <- df3[last(which(df3$fault_occurring  == 1)), ]$Date_Time

# Select variables --------------------------------------------------------
vars_to_monitor <- c("O3-Gas Feed Concentration (mg/L)",
                     "BAF-Turbidity [2] (NTU)",
                     "UF-Feed Flow (GPM)",
                     "TMP (psi)",
                     "UF-Permeate Flow (GPM)",
                     "UF-Reject Pressure (psi)",
                     "UF-Current Draw (A)",
                     "UF-Turbidity [3] (NTU)")

vars_to_monitor_clean <- c("O3 Gas Feed Concentration",
                           "BAF Turbidity",
                           "UF Feed Flow",
                           "TMP",
                           "UF Permeate Flow",
                           "UF Reject Pressure",
                           "UF Current Draw",
                           "UF Turbidity")

clean_df <- function(df) {
  df[, c(vars_to_monitor, "UF State", "uf_state_iter", "Date_Time")] |> 
    #mutate(TMP = map_dbl(TMP, \(x) ifelse(x < 0, 0, x))) |> 
    filter(!(`UF State` %in% c("Manual Control", "Startup", "Start BW", "other")))
}

df3_clean <- clean_df(df3)

df3_pre <- clean_df(df3_pre_raw)
df3_drn <- clean_df(df3_drn_raw)

# Heat maps ---------------------------------------------------------------

# Make with pre and during data
ggheat <- function(df, method = "Sample") {
  if(method == "Sample") {
    df <- cor(df[, vars_to_monitor])
  } else if(method == "LASSO") {
    df <- spcov(cov(df[, vars_to_monitor]),
                cov(df[, vars_to_monitor]),
                lambda = 0.2,
                step.size = 100)$Sigma |> cov2cor()
  } else if(method == "COMET") {
    df <- COmet(as.matrix(df[, vars_to_monitor]), lambda = 0.05)$cov_list[[1]] |> cov2cor()
  }
  
  colnames(df) <- vars_to_monitor_clean
  
  df |> as_tibble() |> 
    mutate(other_var = as_factor(vars_to_monitor_clean), .before = "TMP") |> 
    pivot_longer(-other_var, names_transform = as_factor) |> 
    ggplot(aes(name, other_var, fill = value)) +
    geom_raster() +
    scale_y_discrete(limits = rev(vars_to_monitor_clean)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) #+
  #labs(title = glue("{method} Covariance Heat Map"), x = "", y = "", fill = "")
}

generate_ggheat_plots <- function(df, method = "Sample", data_lab = "") {
  map(1:length(uf_states_reduced), \(i) {
    df |> 
      filter(`UF State` == uf_states_reduced[[i]]) |> 
      ggheat(method = method) +
      labs(title = glue("{data_lab} {method} Correlation: {uf_states_reduced_nice_names[[i]]}"), x = "", y = "", fill = "")
  })
}

ggheat_comparison <- function(p_list1, p_list2) {
  (p_list1[[1]] + p_list2[[1]]) /
    #(p_list1[[2]] + p_list2[[2]]) /
    (p_list1[[3]] + p_list2[[3]]) /
    (p_list1[[4]] + p_list2[[4]]) /
    (p_list1[[5]] + p_list2[[5]]) +
    plot_layout(guides = "collect")
}

# Sample cov heat map by UF State
heat_pre <- generate_ggheat_plots(df3_pre, method = "Sample", data_lab = "Pre-Fault")
heat_drn <- generate_ggheat_plots(df3_drn, method = "Sample", data_lab = "Faulty")

heat_pre_fit1 <- generate_ggheat_plots(df3_pre_fit1$residuals, method = "Sample", data_lab = "Pre-Fault (Centered)")
heat_drn_fit1 <- generate_ggheat_plots(df3_drn_fit1$residuals, method = "Sample", data_lab = "Faulty (Centered)")

heat_pre_fit2 <- generate_ggheat_plots(df3_pre_fit2$residuals, method = "Sample", data_lab = "Pre-Fault (Detrended)")
heat_drn_fit2 <- generate_ggheat_plots(df3_drn_fit2$residuals, method = "Sample", data_lab = "Faulty (Detrended)")

# COMET cov heat map by UF State
heat_pre_comet <- generate_ggheat_plots(df3_pre, method = "COMET", data_lab = "Pre-Fault")
heat_drn_comet <- generate_ggheat_plots(df3_drn, method = "COMET", data_lab = "Faulty")

heat_pre_fit1_comet <- generate_ggheat_plots(df3_pre_fit1$residuals, method = "COMET", data_lab = "Pre-Fault (Centered)")
heat_drn_fit1_comet <- generate_ggheat_plots(df3_drn_fit1$residuals, method = "COMET", data_lab = "Faulty (Centered)")

heat_pre_fit2_comet <- generate_ggheat_plots(df3_pre_fit2$residuals, method = "COMET", data_lab = "Pre-Fault (Detrended)")
heat_drn_fit2_comet <- generate_ggheat_plots(df3_drn_fit2$residuals, method = "COMET", data_lab = "Faulty (Detrended)")

# Save heatmaps
w_heat <- 4000
h_heat <- 5000

ggheat_comparison(heat_pre, heat_drn); ggsave(here("figures", "heat_raw_sample.png"), width = w_heat, height = h_heat, units = "px")
ggheat_comparison(heat_pre_fit1, heat_drn_fit1); ggsave(here("figures", "heat_fit1_sample.png"), width = w_heat, height = h_heat, units = "px")

ggheat_comparison(heat_pre_comet, heat_drn_comet); ggsave(here("figures", "heat_raw_comet.png"), width = w_heat, height = h_heat, units = "px")
ggheat_comparison(heat_pre_fit1_comet, heat_drn_fit1_comet); ggsave(here("figures", "heat_fit1_comet.png"), width = w_heat, height = h_heat, units = "px")


heat_pre[[1]] +
  heat_drn[[1]] +
  heat_pre_comet[[1]] +
  heat_drn_comet[[1]] +
  plot_layout(nrow = 2, ncol = 2, guides = "collect")

ggsave(here("figures", "heat_run.png"), width = 3500, height = 2500, units = "px")

heat_pre[[5]] +
  heat_drn[[5]] +
  heat_pre_comet[[5]] +
  heat_drn_comet[[5]] +
  plot_layout(nrow = 2, ncol = 2, guides = "collect")

ggsave(here("figures", "heat_purge.png"), width = 3500, height = 2500, units = "px")

# Monitoring helper function ----------------------------------------------
monitor_df <- function(dat, method = "hawkins", n_trn = floor(nrow(df)/2),
                       alpha = .1, beta = .1,
                       lambda_1 = .2, lambda_2 = .2, #cutoff = .1
                       lambda_s = .2, cutoff = .1, n_w = 275) {#n_w = 60
  # Remove excess variables
  df <- dat[, vars_to_monitor]
  
  # Calculate pstat based on 
  if(method == "hawkins") {
    mu_0    <- colMeans(df[1:n_trn, ])
    sigma_0 <- cov(df[1:n_trn, ])
    
    pstat <- hawkins_2008(df, mu_0, sigma_0, beta = beta)
    
  } else if(method == "wang_1") {
    mu_0    <- colMeans(df[1:n_trn, ])
    sigma_0 <- cov(df[1:n_trn, ])
    
    pstat <- wang_2014_1(df, mu_0, sigma_0, 
                         alpha = alpha,
                         beta = beta,
                         lambda_1 = lambda_1,
                         lambda_2 = lambda_2)
    
  } else if(method == "wang_2") {
    mu_0    <- colMeans(df[1:n_trn, ])
    sigma_0 <- cov(df[1:n_trn, ])
    
    pstat <- wang_2014_2(df, mu_0, sigma_0, 
                         alpha = alpha,
                         beta = beta,
                         lambda_1 = lambda_1,
                         lambda_2 = lambda_2)
    
    pstat <- tibble(T_1 = pstat$T_1,
                    T_2 = pstat$T_2)
    
  } else if(method == "MC_LASSO") {
    mu_0    <- colMeans(df[1:n_trn, ])
    sigma_0 <- cov(df[1:n_trn, ])#spcov(cov(df[1:n, ]), cov(df[1:n, ]),
    #lambda = lambda_s, step.size = 100)$Sigma
    
    pstat <- MC_LASSO(df, mu_0, sigma_0,
                      beta = beta,
                      lambda_s = lambda_s)
    
  } else if(method == "MC_COMET") {
    mu_0    <- colMeans(df[1:n_trn, ])
    #S_0 <- cov(df[1:n, ])
    #S_0[(abs(S_0) < cutoff)] <- 0
    #sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
    sigma_0 <- cov(df[1:n_trn, ])
    
    pstat <- MC_COMET(df, mu_0, sigma_0,
                      beta = beta,
                      cutoff = cutoff,
                      n_w = n_w)
    
  } else if(method == "MAC_COMET") {
    mu_0    <- colMeans(df[1:n_trn, ])
    #S_0 <- cov(df[1:n, ])
    #S_0[(abs(S_0) < selected_method$cutoff)] <- 0
    #sigma_0 <- covchaud(S_0, as.matrix(df[1:n, ]))$mat
    sigma_0 <- cov(df[1:n_trn, ])
    
    pstat <- MAC_COMET(df, mu_0, sigma_0,
                       alpha = alpha,
                       beta = beta,
                       cutoff = cutoff,
                       n_w = n_w)
    
    pstat <- tibble(T_1 = pstat$T_1,
                    T_2 = pstat$T_2)
  }
  
  if(method %in% c("hawkins", "wang_1", "MC_LASSO", "MC_COMET")) {
    df_return <- 
      pstat |> 
      mutate(Date_Time = dat$Date_Time,
             `UF State` = dat$`UF State`,
             h = estimate_h_fd(pstat[1:floor(n()/2)], alpha))
    
  } else if(method %in% c("wang_2", "MAC_COMET")) {
    df_return <- 
      pstat |> 
      mutate(Date_Time = dat$Date_Time,
             `UF State` = dat$`UF State`,
             h1 = estimate_h2_fd(T_1[1:floor(n()/2)],
                                 T_2[1:floor(n()/2)], alpha)[1],
             h2 = estimate_h2_fd(T_1[1:floor(n()/2)],
                                 T_2[1:floor(n()/2)], alpha)[2])
  }
  
  df_return
  
}

get_name <- function(method) {
  case_match(method,
             "MC_LASSO" ~ "MC-LASSO", 
             "MC_COMET" ~ "MC-COMET",
             "MAC_COMET" ~ "MAC-COMET",
             "hawkins" ~ "MC-Hawkins", 
             "wang_1" ~ "MC-Wang",
             "wang_2" ~ "MAC-Wang",
             .default = method)
}

get_fault_sections <- function(df) {
  # Identifies pre, during, and post fault sections
  fault_section <- vector("integer", nrow(df))
  temp <- 0
  
  for(i in 2:nrow(df)) {
    if(df$fault_occurring[i-1] != df$fault_occurring[i]) {
      temp <- temp + 1
    }
    fault_section[i] <- temp
  }
  
  fault_section
}

# var should be a single character vector
get_df_rect <- function(df, var, y_min = NA, y_max = NA) {
  # Determine ymin and max
  if(!is.na(y_min)) {
    ymin = y_min
    ymax = y_max
  } else if(is.numeric(unlist(df[, var]))) {
    ymin = min(df[, var], na.rm = TRUE)
    ymax = max(df[, var], na.rm = TRUE)
  } else {
    ymin = first(unique(unlist(df[, var])))
    ymax = last(unique(unlist(df[, var])))
  }
  
  # Create df to define rectagles
  summarise(df,
            fault_occurring = first(fault_occurring),
            xmin = max(c(min(Date_Time), time_exp_3_limits[1])), xmax = min(c(max(Date_Time), time_exp_3_limits[2])),
            ymin = ymin, ymax = ymax,
            .by = fault_section) |> 
    mutate(fault_occurring = as_factor(fault_occurring))
}

# var should be a single character vector
plot_ts <- function(df, var, type = "line", y_min = NA, y_max = NA) {
  # Create plot template
  plot <- 
    df |> 
    ggplot() + 
    geom_rect(data = get_df_rect(df, var, y_min, y_max), 
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fault_occurring |> 
                    case_match("0" ~ "No Fault", "1" ~ "Faulty") |> 
                    as_factor()),
              alpha = .5) +
    scale_x_datetime(limits = c(),### NOT DONE
                     date_breaks = "12 hours",
                     date_labels = "%m/%d%n%H:%M") +
    labs(x = "", 
         y = ifelse(str_detect(var, "\\("),
                    str_extract(var, "\\(.+\\)"),
                    ""),
         fill = "",
         title = str_replace(str_split(var, "\\(|\\[")[[1]], "-", " ")) + 
    scale_fill_brewer(palette = "Paired")
  #scale_fill_viridis_d(begin = 0, end = 1, option = "viridis")
  
  # Add line or points
  if(type == "line") {
    plot + 
      geom_line(aes(Date_Time, !! sym(var))) + 
      ylim(y_min, y_max)
  } else if(type == "point") {
    plot + geom_point(aes(Date_Time, !! sym(var)), shape = 21)
  }
}

gg_pstat <- function(df, method = "hawkins", data = "Raw") {
  df[[method]] <- df[[method]] |> arrange(Date_Time)
  
  df[[method]] <-
    df[[method]]  |> 
    mutate(fault_occurring = between(Date_Time, begin_fault, end_fault),
           fault_section = 
             get_fault_sections(
               df[[method]] |> mutate(fault_occurring = between(Date_Time, begin_fault, end_fault))
             ))
  
  if(method %in% c("wang_2", "MAC_COMET")) {
    p1 <- 
      df[[method]] |> 
      ggplot() +
      geom_rect(data = get_df_rect(df[[method]], "T_1"),
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fault_occurring |>
                      case_match("FALSE" ~ "No Fault", "TRUE" ~ "Faulty") |>
                      as_factor()),
                alpha = .5) +
      scale_fill_brewer(palette = "Paired") +
      geom_line(aes(Date_Time, T_1, color = `UF State`, group = 1)) +
      scale_color_viridis_d() +
      scale_x_datetime(limits = time_exp_3_limits,
                       date_breaks = "12 hours",
                       date_labels = "%m/%d%n%H:%M") +
      geom_line(aes(x = Date_Time, y = h1), color = "black") + 
      labs(x = "", y = "", fill = "", title = glue("{get_name(method)} Charting Statistic: T_1"))
    
    p2 <- 
      df[[method]] |> 
      ggplot() +
      geom_rect(data = get_df_rect(df[[method]], "T_2"),
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fault_occurring |>
                      case_match("FALSE" ~ "No Fault", "TRUE" ~ "Faulty") |>
                      as_factor()),
                alpha = .5) +
      scale_fill_brewer(palette = "Paired") +
      geom_line(aes(Date_Time, T_2, color = `UF State`, group = 1)) +
      scale_color_viridis_d() +
      scale_x_datetime(limits = time_exp_3_limits,
                       date_breaks = "12 hours",
                       date_labels = "%m/%d%n%H:%M") +
      geom_line(aes(x = Date_Time, y = h2), color = "black") +
      labs(x = "", y = "", fill = "", title = glue("{get_name(method)} Charting Statistic: T_2"))
    
    list(p1, p2)
  } else {
    df[[method]] |> 
      ggplot() +
      geom_rect(data = get_df_rect(df[[method]], "pstat"),
                aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax,
                    fill = fault_occurring |>
                      case_match("FALSE" ~ "No Fault", "TRUE" ~ "Faulty") |>
                      as_factor()),
                alpha = .5) +
      scale_fill_brewer(palette = "Paired") +
      geom_line(aes(Date_Time, pstat, color = `UF State`, group = 1)) +
      scale_color_viridis_d() +
      scale_x_datetime(limits = time_exp_3_limits,
                       date_breaks = "12 hours",
                       date_labels = "%m/%d%n%H:%M") +
      geom_line(aes(x = Date_Time, y = h), color = "black") +
      labs(x = "", y = "", fill = "", title = glue("{get_name(method)} Charting Statistic"))
  }
}

# Parameters
#n_trn <- 33353
methods <- c("hawkins", "MAC_COMET", "wang_2")
alpha = 0.005
height_1 <- 3500
width_1 <- 4000

# Monitor raw data irrespective of system state ---------------------------
df_raw_irr <- 
  map(methods, \(m) {
    monitor_df(df3_clean[-c(41116:41197), ], method = m)
  }) |> set_names(methods)

saveRDS(df_raw_irr, here("data_dpr", "df_raw_irr.rds"))
df_raw_irr <- readRDS(here("data_dpr", "df_raw_irr.rds"))

wrap_plots(list(gg_pstat(df_raw_irr, method = "hawkins"),
                gg_pstat(df_raw_irr, method = "MAC_COMET")[[1]],
                gg_pstat(df_raw_irr, method = "MAC_COMET")[[2]],
                gg_pstat(df_raw_irr, method = "wang_2")[[1]],
                gg_pstat(df_raw_irr, method = "wang_2")[[2]]
                ),
           ncol = 1, 
           guides = "collect")

ggsave(here("figures", "pstat_raw_irr.pdf"), width = width_1, height = height_1, units = "px")

# Monitor raw data conditional on the system state ------------------------
df_raw_cond_split <- split(df3_clean, df3_clean$`UF State`, drop = TRUE)
df_raw_cond_split[[1]] <- df_raw_cond_split[[1]][-c(36718:36835), ]

# 36761
df_raw_cond <- 
  map(methods, \(m) {
    df_raw_cond_split |> 
      map_dfr(\(data) {
        monitor_df(data, method = m)
      })
  }) |> set_names(methods)

df_raw_cond$MAC_COMET <- df_raw_cond_split |> 
  map_dfr(\(data) {
    monitor_df(data, method = "MAC_COMET")
  })

saveRDS(df_raw_cond, here("data_dpr", "df_raw_cond.rds"))
df_raw_cond <- readRDS(here("data_dpr", "df_raw_cond.rds"))

wrap_plots(list(gg_pstat(df_raw_cond, method = "hawkins"),
                gg_pstat(df_raw_cond, method = "MAC_COMET")[[1]],
                gg_pstat(df_raw_cond, method = "MAC_COMET")[[2]],
                gg_pstat(df_raw_cond, method = "wang_2")[[1]],
                gg_pstat(df_raw_cond, method = "wang_2")[[2]]
                ),
           ncol = 1, guides = "collect")

ggsave(here("figures", "pstat_raw_cond.pdf"), width = width_1, height = height_1, units = "px")


# Detection rates ---------------------------------------------------------
calc_fault_state <- function(date_time, begin_date, end_date) {
  if(date_time <= begin_date) {
    "Pre-Fault"
  } else if(date_time <= end_date) {
    "Faulty"
  } else {
    "Post-Fault"
  }
}

get_dpr_dr <- function(df_list) {
  df_list |> 
    map(\(df) {
      if("pstat" %in% names(df)) {
        df |> 
          mutate(fault_state = map_chr(Date_Time, \(x) calc_fault_state(x, begin_fault, end_fault)),
                 test_state = pstat > h) |> 
          summarise(mean = mean(test_state), .by = fault_state)
      } else {
        df |> 
          mutate(T_2 = replace_na(T_2, 0),
                 fault_state = map_chr(Date_Time, \(x) calc_fault_state(x, begin_fault, end_fault)),
                 test_state = T_1 > h1 | T_2 > h2) |> 
          summarise(mean = mean(test_state), .by = fault_state)
      }
    }) |> 
    bind_rows(.id = "method")
}

get_dpr_conditional_dr <- function(df_list) {
  df_list |> 
    map(\(df) {
      if("pstat" %in% names(df)) {
        df |> 
          mutate(fault_state = map_chr(Date_Time, \(x) calc_fault_state(x, begin_fault, end_fault)),
                 test_state = pstat > h) |> 
          summarise(mean = mean(test_state), .by = c(fault_state, `UF State`))
      } else {
        df |> 
          mutate(T_2 = replace_na(T_2, 0),
                 fault_state = map_chr(Date_Time, \(x) calc_fault_state(x, begin_fault, end_fault)),
                 test_state = T_1 > h1 | T_2 > h2) |> 
          summarise(mean = mean(test_state), .by = c(fault_state, `UF State`))
      }
    }) |> 
    bind_rows(.id = "method")
}

# Detection Rate Table
get_dpr_dr(df_raw_irr) 
get_dpr_dr(df_raw_cond)

# Formatted Detection Rate Table
dr_table <- 
  bind_rows(Irrespective = get_dpr_dr(df_raw_irr), 
            Conditional = get_dpr_dr(df_raw_cond), .id = "type") |>
  filter(fault_state != "Post-Fault",
         method %in% methods) |> 
  mutate(fault_state = as_factor(fault_state),
         across(where(is.numeric), \(x) round(x, 3))) |> 
  relocate(fault_state, .before = method) |> 
  arrange(fault_state) |> 
  pivot_wider(values_from = mean, names_from = type)

dr_table |> 
  kableExtra::kbl(format = "latex", booktabs = TRUE, linesep =  c('', '', '\\addlinespace'))

# Detection Rate Table (Conditional)
dr_conditional_table <- 
  get_dpr_conditional_dr(df_raw_cond) |> 
  filter(fault_state != "Post-Fault",
         method %in% methods) |> 
  mutate(fault_state = as_factor(fault_state),
         across(where(is.numeric), \(x) round(x, 3))) |> 
  relocate(fault_state, .before = method) |> 
  arrange(fault_state) |> 
  pivot_wider(values_from = mean, names_from = `UF State`)

dr_conditional_table |> view()