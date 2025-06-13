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
library("forecast")
library("MTS")
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
                                  `Pre-BW Air Scour` = "AS_pre",
                                  `Post-BW Air Scour` = "AS_post",
                                  `Backwash` = "BW",
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
  

df3 <- df3 |> 
  mutate(TMP = `UF-Feed Pressure (psi)` - `UF-Permeate Pressure (psi)`,
         `UF State` = as_factor(uf_state_vec_reduced),
         uf_state_iter = uf_state_iter)
  #filter(`UF State` == "Run") # Restricts data to when UF is running automatically


# Pre, during, and post fault data frames (BY DATES)
# df3_pre <- df3 |> filter(Date_Time <= ymd_hms("2023-06-19 04:09:00"))
# df3_drn <- df3 |> filter(between(Date_Time, 
#                                  ymd_hms("2023-06-19 04:09:00"),
#                                  ymd_hms("2023-06-21 09:03:00")))
# df3_pst <- df3 |> filter(Date_Time >= ymd_hms("2023-06-21 09:03:00"))

# Pre, during, and post fault data frames (by fault_occuring variable)
df3_pre <- df3[1:(first(which(df3$fault_occurring == 1))-1), ]
df3_drn <- df3[which(df3$fault_occurring == 1), ]
df3_pst <- df3[(last(which(df3$fault_occurring == 1))+1):nrow(df3), ]



# Select variables --------------------------------------------------------
vars_to_monitor <- c("TMP",
                     "UF-Feed Flow (GPM)",
                     "UF-Permeate Flow (GPM)",
                     "UF-Reject Pressure (psi)",
                     "UF-Current Draw (A)",
                     "BAF-Turbidity [2] (NTU)",
                     "UF-Turbidity [3] (NTU)",
                     "EFF-Turbidity [4] (NTU)",
                     "O3-Gas Feed Concentration (mg/L)")

clean_df <- function(df) {
  df[, c(vars_to_monitor, "UF State", "uf_state_iter", "Date_Time")] |> 
    mutate(TMP = map_dbl(TMP, \(x) ifelse(x < 0, 0, x))) |> 
    filter(!(`UF State` %in% c("Manual Control", "Startup", "Start BW", "other")))
    #filter(`UF State` %in% c("Run", "AS/BW"))## |> 
    # mutate(is_backwash = str_detect(`UF State`, "AS/BW"),
    #        is_purge = str_detect(`UF State`, "Purge"))
}


df3_pre_fit <- clean_df(df3_pre)
df3_drn_fit <- clean_df(df3_drn)


# Modeling monitored variables --------------------------------------------
# Center/scale
df3_centerscale <- 
  df3_pre_fit |> 
  mutate(across(all_of(vars_to_monitor), \(x) {(x - mean(x))/sd(x)}))
  
# Linear model for each UF State, AR(1) for others
uf_vars <- vars_to_monitor[1:5]
other_vars <- vars_to_monitor[6:10]

fit_1 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ `UF State`")),
       df3_pre_fit)
  })

# ar_1 <- other_vars |> 
#   map(\(x) {
#     arima(df3_pre_fit[[x]], order = c(1, 0, 0))
#   })

df3_fit_1 <- c(fit_1) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time,
         `UF State` = df3_pre_fit$`UF State`)

# Linear model for each UF State + uf_state_iter, AR(1) for others
fit_2 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ `UF State` + uf_state_iter")),
       df3_pre_fit)
  })

df3_fit_2 <- c(fit_2) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time,
         `UF State` = df3_pre_fit$`UF State`)

# Linear model for each UF State*uf_iter, AR(1) for others
fit_3 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ `UF State` * uf_state_iter")),
       df3_pre_fit)
  })

df3_fit_3 <- c(fit_3) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time,
         `UF State` = df3_pre_fit$`UF State`)

# Linear model by lag(1)
fit_4 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ lag(`", x, "`)")),
       df3_pre_fit)
  })

df3_fit_4 <- c(fit_4) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time[-1],
         `UF State` = df3_pre_fit$`UF State`[-1])

# Linear model by lag(1) + uf state
fit_5 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ lag(`", x, "`) + `UF State`")),
       df3_pre_fit)
  })

df3_fit_5 <- c(fit_5) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time[-1],
         `UF State` = df3_pre_fit$`UF State`[-1])

# Linear model by lag(1) * uf state
fit_6 <- vars_to_monitor |> 
  map(\(x) {
    lm(as.formula(paste0("`",x, "` ~ lag(`", x, "`) * `UF State`")),
       df3_pre_fit)
  })

df3_fit_6 <- c(fit_6) |> 
  map_dfc(\(x) {
    x$residuals
  }) |> set_names(vars_to_monitor) |> 
  mutate(Date_Time = df3_pre_fit$Date_Time[-1],
         `UF State` = df3_pre_fit$`UF State`[-1])
## VMA
## VAR


# Predicting monitored variables ------------------------------------------

# Center/scale
df3_tst_centerscale <- 
  bind_cols(df3_drn_fit |> select(!all_of(vars_to_monitor)),
            map_dfc(vars_to_monitor, \(name) {
              (df3_drn_fit[[name]] - mean(df3_pre_fit[[name]]))/sd(df3_pre_fit[[name]])
            }) |> set_names(vars_to_monitor))



pred_lm <- function(fit) {
  map2_dfc(fit, vars_to_monitor, \(m, var) {
    df3_drn_fit[[var]] - predict(m, df3_drn_fit)
  }) |> set_names(vars_to_monitor) |> 
    mutate(Date_Time = df3_drn_fit$Date_Time,
           `UF State` = df3_drn_fit$`UF State`)
}

df3_tst_fit_1 <- pred_lm(fit_1)
df3_tst_fit_2 <- pred_lm(fit_2)
df3_tst_fit_3 <- pred_lm(fit_3)
df3_tst_fit_4 <- pred_lm(fit_4)
df3_tst_fit_5 <- pred_lm(fit_5)
df3_tst_fit_6 <- pred_lm(fit_6)
  

# df3_tst_fit_6 <-
#   map2_dfc(fit_6, vars_to_monitor, \(m, var) {
#     df3_drn_fit[[var]] - predict(m, df3_drn_fit)
#   }) |> set_names(vars_to_monitor) |> 
#   mutate(Date_Time = df3_drn_fit$Date_Time,
#          `UF State` = df3_drn_fit$`UF State`)

# Time series plots -------------------------------------------------------
ggtime <- function(df) {
  df |> 
    pivot_longer(all_of(vars_to_monitor), names_transform = as_factor) |> 
    ggplot(aes(Date_Time, value, color = `UF State`)) +
    geom_point(shape = 21) +
    facet_wrap(~ name, ncol = 1, scales = "free")
}

ggtime2 <- function(df1, df2, n_vars = 1:9, n_viz = NULL) {
  if(is.numeric(n_viz)) {
    bind_rows(list("Pre-Fault" = df1, 
                   "Faulty" = df2), .id = "DF_NAME") |> 
      mutate(DF_NAME = as_factor(DF_NAME)) |> 
      slice((max(which(DF_NAME == "Pre-Fault"))-n_viz):(max(which(DF_NAME == "Pre-Fault"))+n_viz + 1)) |> 
      pivot_longer(all_of(vars_to_monitor[n_vars]), names_transform = as_factor) |> 
      ggplot(aes(Date_Time, value, color = `UF State`)) +
      geom_point(shape = 21) +
      facet_grid(name ~ DF_NAME, scales = "free") +
      labs(x = "", y = "")
  } else {
    bind_rows(list("Pre-Fault" = df1, 
                   "Faulty" = df2), .id = "DF_NAME") |> 
      mutate(DF_NAME = as_factor(DF_NAME)) |> 
      pivot_longer(all_of(vars_to_monitor[n_vars]), names_transform = as_factor) |> 
      ggplot(aes(Date_Time, value, color = `UF State`)) +
      geom_point(shape = 21) +
      facet_grid(name ~ DF_NAME, scales = "free") +
      labs(x = "", y = "")
  }
}

ggtime2(df3_pre_fit, df3_drn_fit)

df_test <-
bind_rows(list("Pre-Fault" = df3_pre_fit, 
               "Faulty" = df3_pre_fit), .id = "DF_NAME")

max(which(df_test$DF_NAME == "Pre-Fault"))
which.min(df_test$DF_NAME == "Pre-Fault")
  

n_viz <- 1000
w <- 4000
h <- 5000

# Point plot, full data
ggtime(df3_pre) + labs(title = "Pre-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pre-full.png"), width = w, height= h, units = "px")

ggtime(df3_drn) + labs(title = "Faulty Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-drn-full.png"), width = w, height= h, units = "px")

ggtime(df3_pst) + labs(title = "Post-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pst-full.png"), width = w, height= h, units = "px")

# Point plot, reduced data
ggtime(tail(df3_pre, n_viz)) + labs(title = glue("Last {n_viz} Pre-Fault Observations"), y = "")
ggsave(here("figures", "dpr-ts-pre.png"), width = w, height= h, units = "px")

ggtime(tail(df3_drn, n_viz)) + labs(title = glue("First {n_viz} Faulty Observations"), y = "")
ggsave(here("figures", "dpr-ts-drn.png"), width = w, height= h, units = "px")

# Time series plot of residuals -------------------------------------------

# Side by side (full)
ggtime2(df3_pre_fit, df3_drn_fit); ggsave(here("figures", "dpr_ts_full.png"), width = w, height= h, units = "px")
ggtime2(df3_centerscale, df3_tst_centerscale); ggsave(here("figures", "dpr_centerscale.png"), width = w, height= h, units = "px")

# Side by side (reduced)
ggtime2(df3_pre_fit, df3_drn_fit, n_viz = n_viz); ggsave(here("figures", "dpr_ts_reduced.png"), width = w, height= h, units = "px")


ggtime(tail(df3_centerscale, n_viz)) + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of Center/Scale"), y = ""); ggsave(here("figures", "dpr_centerscale_pre.png"), width = w, height= h, units = "px")
ggtime(tail(df3_fit_1, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_1_pre.png"), width = w, height= h, units = "px") 
ggtime(tail(df3_fit_2, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ `UF State` + Iterations"), y = ""); ggsave(here("figures", "dpr_fit_2_pre.png"), width = w, height= h, units = "px")
ggtime(tail(df3_fit_3, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ `UF State` * Iterations"), y = ""); ggsave(here("figures", "dpr_fit_3_pre.png"), width = w, height= h, units = "px")
ggtime(tail(df3_fit_4, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ lag(variable)"), y = ""); ggsave(here("figures", "dpr_fit_4_pre.png"), width = w, height= h, units = "px")
ggtime(tail(df3_fit_5, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ lag(variable) + `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_5_pre.png"), width = w, height= h, units = "px")
ggtime(tail(df3_fit_6, n_viz))  + labs(title = glue("Pre-Fault (Last {n_viz}): Residuals of ~ lag(variable) * `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_6_pre.png"), width = w, height= h, units = "px")

# Testing Data
ggtime(head(df3_drn_fit, n_viz)) + labs(title = "Faulty (Last {n_viz}): Raw data", y = "")
ggtime(head(df3_tst_centerscale, n_viz)) + labs(title = glue("Faultly (First {n_viz}): Residuals of Center/Scale"), y = ""); ggsave(here("figures", "dpr_centerscale_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_1, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_1_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_2, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ `UF State` + Iterations"), y = ""); ggsave(here("figures", "dpr_fit_2_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_3, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ `UF State` * Iterations"), y = ""); ggsave(here("figures", "dpr_fit_3_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_4, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ lag(variable)"), y = ""); ggsave(here("figures", "dpr_fit_4_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_5, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ lag(variable) + `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_5_drn.png"), width = w, height= h, units = "px")
ggtime(head(df3_tst_fit_6, n_viz))  + labs(title = glue("Faultly (First {n_viz}): Residuals of ~ lag(variable) * `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_6_drn.png"), width = w, height= h, units = "px")


# Training Data (FULL)
ggtime(df3_pre_fit) + labs(title = glue("Pre-Fault (All Observations): Raw data"), y = "")
ggtime(df3_centerscale) + labs(title = glue("Pre-Fault (All Observations): Residuals of Center/Scale"), y = ""); ggsave(here("figures", "dpr_centerscale_pre_full.png"), width = w, height= h, units = "px")
ggtime(df3_fit_1)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_1_pre_full.png"), width = w, height= h, units = "px") 
ggtime(df3_fit_2)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ `UF State` + Iterations"), y = ""); ggsave(here("figures", "dpr_fit_2_pre_full.png"), width = w, height= h, units = "px")
ggtime(df3_fit_3)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ `UF State` * Iterations"), y = ""); ggsave(here("figures", "dpr_fit_3_pre_full.png"), width = w, height= h, units = "px")
ggtime(df3_fit_4)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ lag(variable)"), y = ""); ggsave(here("figures", "dpr_fit_4_pre_full.png"), width = w, height= h, units = "px")
ggtime(df3_fit_5)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ lag(variable) + `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_5_pre_full.png"), width = w, height= h, units = "px")
ggtime(df3_fit_6)  + labs(title = glue("Pre-Fault (All Observations): Residuals of ~ lag(variable) * `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_6_pre_full.png"), width = w, height= h, units = "px")

# Testing Data (FULL)
ggtime(df3_drn_fit) + labs(title = "Faulty (All Observations): Raw data", y = "")
ggtime(df3_tst_centerscale) + labs(title = glue("Faultly (All Observations): Residuals of Center/Scale"), y = ""); ggsave(here("figures", "dpr_centerscale_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_1)  + labs(title = glue("Faultly (All Observations): Residuals of ~ `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_1_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_2)  + labs(title = glue("Faultly (All Observations): Residuals of ~ `UF State` + Iterations"), y = ""); ggsave(here("figures", "dpr_fit_2_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_3)  + labs(title = glue("Faultly (All Observations): Residuals of ~ `UF State` * Iterations"), y = ""); ggsave(here("figures", "dpr_fit_3_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_4)  + labs(title = glue("Faultly (All Observations): Residuals of ~ lag(variable)"), y = ""); ggsave(here("figures", "dpr_fit_4_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_5)  + labs(title = glue("Faultly (All Observations): Residuals of ~ lag(variable) + `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_5_drn_full.png"), width = w, height= h, units = "px")
ggtime(df3_tst_fit_6)  + labs(title = glue("Faultly (All Observations): Residuals of ~ lag(variable) * `UF State`"), y = ""); ggsave(here("figures", "dpr_fit_6_drn_full.png"), width = w, height= h, units = "px")

# Heat map of covariance --------------------------------------------------

# Produce heat map
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
  
  df |> as_tibble() |> 
    mutate(other_var = as_factor(vars_to_monitor), .before = "TMP") |> 
    pivot_longer(-other_var, names_transform = as_factor) |> 
    ggplot(aes(name, other_var, fill = value)) +
    geom_raster() +
    scale_y_discrete(limits = rev(vars_to_monitor)) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
    #scale_fill_viridis_c(option = "plasma", limits = c(-1, 1)) + # if using correlation, limits = c(-1, 1)
    theme(axis.text.x = element_text(angle = 45, hjust=1)) #+
    #labs(title = glue("{method} Covariance Heat Map"), x = "", y = "", fill = "")
}

df_comet_test <- COmet(as.matrix(df3_pre[, vars_to_monitor]), lambda = 0.05)$cov_list[[1]]

df_comet_test |> cov2cor() |> as_tibble() |> 
  mutate(other_var = as_factor(vars_to_monitor), .before = "TMP") |> 
  pivot_longer(-other_var, names_transform = as_factor)|> 
  ggplot(aes(name, other_var, fill = value)) +
  geom_raster() +
  scale_y_discrete(limits = rev(vars_to_monitor)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1)) +
  #scale_fill_viridis_c(option = "plasma", limits = c(-1, 1)) + # if using correlation, limits = c(-1, 1)
  theme(axis.text.x = element_text(angle = 45, hjust=1))

w_heat <- 4000
h_heat <- 5000
# Sample cov heat map
# w_heat <- 2200
# h_heat <- 1800
# ggheat(df3_pre) + labs(title = "Pre-Fault Sample Correlation", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_sample_cov_pre.png"), width = w_heat, height= h_heat, units = "px")
# ggheat(df3_drn) + labs(title = "Faulty Sample Correlation", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_sample_cov_drn.png"), width = w_heat, height= h_heat, units = "px")
# 
# # LASSO cov heat map
# ggheat(df3_pre, method = "LASSO") + labs(title = "Pre-Fault LASSO Covariance", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_lasso_cov_pre.png"), width = w_heat, height= h_heat, units = "px")
# ggheat(df3_drn, method = "LASSO") + labs(title = "Faulty LASSO Covariance", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_lasso_cov_drn.png"), width = w_heat, height= h_heat, units = "px")
# 
# # COMET cov heat map
# ggheat(df3_pre, method = "COMET") + labs(title = "Pre-Fault COMET Correlation", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_comet_cov_pre.png"), width = w_heat, height= h_heat, units = "px")
# ggheat(df3_drn, method = "COMET") + labs(title = "Faulty COMET Correlation", x = "", y = "", fill = ""); ggsave(here("figures", "heatmap_comet_cov_drn.png"), width = w_heat, height= h_heat, units = "px")

# Sample cov heat map by UF State
heat_pre <- c("AS_pre", "BW", "AS_post", "Purge", "Run") |> 
  map(\(name) {
    df3_pre_fit |> 
      filter(`UF State` == name) |> 
      ggheat() +
      labs(title = glue("Pre-Fault Sample Correlation: {name}"), x = "", y = "", fill = "")
  })

# Sample cov heat map by UF State
heat_drn <- c("AS_pre", "BW", "AS_post", "Purge", "Run") |> 
  map(\(name) {
    df3_drn_fit |> 
      filter(`UF State` == name) |> 
      ggheat() +
      labs(title = glue("Faulty Sample Correlation: {name}"), x = "", y = "", fill = "")
  })

wrap_plots(heat_pre, ncol = 2, guides = "collect")
ggsave(here("figures", "heatmap_sample_cov_pre.png"), width = w_heat, height= h_heat, units = "px")

wrap_plots(heat_drn, ncol = 2, guides = "collect")
ggsave(here("figures", "heatmap_sample_cov_drn.png"), width = w_heat, height= h_heat, units = "px")

# COMET heat map by UF State
heat_comet_pre <- c("AS_pre", "BW", "AS_post", "Purge", "Run") |> 
  map(\(name) {
    df3_pre_fit |> 
      filter(`UF State` == name) |> 
      ggheat(method = "COMET") +
      labs(title = glue("Pre-Fault COMET Correlation: {name}"), x = "", y = "", fill = "")
  })

# Sample cov heat map by UF State
heat_comet_drn <- c("AS_pre", "BW", "AS_post", "Purge", "Run") |> 
  map(\(name) {
    df3_drn_fit |> 
      filter(`UF State` == name) |> 
      ggheat(method = "COMET") +
      labs(title = glue("Faulty COMET Correlation: {name}"), x = "", y = "", fill = "")
  })

wrap_plots(heat_comet_pre, ncol = 2, guides = "collect")
ggsave(here("figures", "heatmap_comet_cov_pre.png"), width = w_heat, height= h_heat, units = "px")

wrap_plots(heat_comet_drn, ncol = 2, guides = "collect")
ggsave(here("figures", "heatmap_comet_cov_drn.png"), width = w_heat, height= h_heat, units = "px")