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
source(here("scripts", "pl_fd_helper.R"))
source(here("scripts", "pl_fd_methods.R"))


# Read data ---------------------------------------------------------------
df3 <- readRDS(here("data_dpr", "dpr_fault_exp_3.rds"))

# Identify UF States
uf_states <- c("Manual Control", "Startup",
               "Start BW", "Air Scour", "Backwash", "Purge", "End BW",
               "Run", "other")

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

df3 <- df3 |> 
  mutate(TMP = `UF-Feed Pressure (psi)` - `UF-Permeate Pressure (psi)`,
         `UF State` = as_factor(uf_state_vec))## |> 
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
                     "INF-Turbidity [1] (NTU)",
                     "BAF-Turbidity [2] (NTU)",
                     "UF-Turbidity [3] (NTU)",
                     "EFF-Turbidity [4] (NTU)",
                     "O3-Gas Feed Concentration (mg/L)")




# Modeling monitored variables --------------------------------------------
# TMP
df3_pre$TMP[which(df3_pre$TMP <= 0)] <- 0

df3_pre[1001:2000, ] |> 
  ggplot(aes(Date_Time, TMP)) +
  geom_point(shape = 21)

# UF Current Draw
df3_pre[1001:2000, ] |> 
  ggplot(aes(Date_Time, `UF-Current Draw (A)`)) +
  geom_point(shape = 21)

# Potential methods
## Center/scale
df3_centerscale <- 
  df3 |> 
  filter(`UF State` == "Run") |> 
  mutate(across(all_of(vars_to_monitor), \(x) {(x - mean(x))/sd(x)}))
  
## linear model for each backwash
df3_test <- df3_pre[, c(vars_to_monitor, "UF State")]

map(vars_to_monitor, \(x) {
  lm(x ~ ., df3_test)
})


plot(fit_pre$model$`UF-Current Draw (A)`)
plot(fit_pre$residuals)

qqnorm(fit_pre$model$`UF-Current Draw (A)`)
qqnorm(fit_pre$residuals)


## Hunter's method
## VMA
## VAR
## Nothing

df3_centerscale$`UF-Current Draw (A)` |> mean()


# Monitoring --------------------------------------------------------------




# Time series plots -------------------------------------------------------

ggtime <- function(df) {
  df |> 
    pivot_longer(all_of(vars_to_monitor), names_transform = as_factor) |> 
    ggplot(aes(Date_Time, value)) +
    facet_wrap(~ name, ncol = 2, scales = "free")
}

# Point plot, reduced data
ggtime(df3_pre[1001:2000, ]) + geom_point(shape = 21) + labs(title = "Pre-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pre-full.pdf"), width = 4000, height= 4000, units = "px")

ggtime(df3_drn[1001:2000, ]) + geom_point(shape = 21) + labs(title = "Faulty Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-drn-full.pdf"), width = 4000, height= 4000, units = "px")

ggtime(df3_pst[1001:2000, ]) + geom_point(shape = 21) + labs(title = "Post-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pst-full.pdf"), width = 4000, height= 4000, units = "px")

# Line plot, all data
ggtime(df3_pre) + geom_line() + labs(title = "Pre-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pre.pdf"), width = 4000, height= 4000, units = "px")

ggtime(df3_drn) + geom_line() + labs(title = "Faulty Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-drn.pdf"), width = 4000, height= 4000, units = "px")

ggtime(df3_pst) + geom_line() + labs(title = "Post-Fault Time Series Plot", y = "")
ggsave(here("figures", "dpr-ts-pst.pdf"), width = 4000, height= 4000, units = "px")


# Heat map of covariance --------------------------------------------------

# Produce heat map
ggheat <- function(df, method = "Sample") {
  if(method == "Sample") {
    df <- cor(df[, vars_to_monitor])
  } else if(method == "LASSO") {
    df <- spcov(cov(df[, vars_to_monitor]),
                cov(df[, vars_to_monitor]),
                lambda = 0.2,
                step.size = 100)$Sigma## |> cov2cor()
  } else if(method == "COMET") {
    df <- COmet(as.matrix(df[, vars_to_monitor]), lambda = 0.05)$cov_list[[1]]## |> cov2cor()
  }
  
  df |> as_tibble() |> 
    mutate(other_var = as_factor(vars_to_monitor), .before = "TMP") |> 
    pivot_longer(-other_var, names_transform = as_factor) |> 
    ggplot(aes(name, other_var, fill = value)) +
    geom_raster() +
    scale_y_discrete(limits = rev(vars_to_monitor)) + 
    scale_fill_viridis_c(option = "plasma") + # if using correlation, limits = c(-1, 1)
    theme(axis.text.x = element_text(angle = 45, hjust=1)) #+
    #labs(title = glue("{method} Covariance Heat Map"), x = "", y = "", fill = "")
}

# Sample cov heat map
ggheat(df3_pre) + labs(title = "Pre-Fault Sample Covariance", x = "", y = "", fill = "")
ggheat(df3_drn) + labs(title = "Faulty Sample Covariance", x = "", y = "", fill = "")

# LASSO cov heat map
ggheat(df3_pre, method = "LASSO") + labs(title = "Pre-Fault LASSO Covariance", x = "", y = "", fill = "")
ggheat(df3_drn, method = "LASSO") + labs(title = "Faulty LASSO Covariance", x = "", y = "", fill = "")

# COMET cov heat map
ggheat(df3_pre, method = "COMET") + labs(title = "Pre-Fault COMET Covariance", x = "", y = "", fill = "")
ggheat(df3_drn, method = "COMET") + labs(title = "Faulty COMET Covariance", x = "", y = "", fill = "")
