library("tidyverse")
theme_set(cowplot::theme_cowplot())
theme_update(plot.title = element_text(hjust = 0.5, size = 15),
             plot.subtitle = element_text(hjust = 0.5, size = 10),
             strip.text.x = element_text(size = 15),
             strip.placement = "outside",
             strip.background = element_blank())
library("lubridate")
library("patchwork")


# Helper Variables --------------------------------------------------------

# All potential variables of interest
vars_relevant = c("index",
                  "Date_Time",
                  "Date",
                  "Runtime (hr)",
                  "TMP (psi)",
                  "INF-Filter Pressure (psi)",
                  "O3-Ozone Concentration (ppm)",
                  "O3-Current Draw (A)",
                  "BAF-C1 Flow (GPM)",
                  "BAF-C2 Flow (GPM)",
                  "BAF-C3 Flow (GPM)",
                  "BAF-C4 Flow (GPM)",
                  "BAF_total_flow",
                  "UF-Feed Pressure (psi)",
                  "UF-Permeate Pressure (psi)",
                  "UF-Reject Pressure (psi)",
                  "UF-Feed Volume (gal)",
                  "UF-Feed Flow (GPM)",
                  "UF-Permeate Flow (GPM)",
                  "UF-Recycle Flow (GPM)",
                  "UF-Current Draw (A)",
                  "UF-Feed Temperature (C)",
                  "BW-Backwash Tank Volume (gal)",
                  "INF-Turbidity [1] (NTU)",
                  "BAF-Turbidity [2] (NTU)",
                  "UF-Turbidity [3] (NTU)",
                  "EFF-Turbidity [4] (NTU)",
                  "DFW-pH [1]",
                  "EFF-pH [2]",
                  "DFO-Nitrate (mg/L)",
                  "DFO-NOX (mg/L)",
                  "DFW-Dissolved Oxygen (mg/L)",
                  "O3-Gas Feed Concentration (mg/L)",
                  "O3-Offgas Concentration (mg/L)",
                  "BAF/UF/GAC/UV/EFF-TOC (ppm)",
                  "UF-Permeate pH",
                  "EFF-Pressure (psi)",
                  "BAF-C1 Flow Valve (%)",
                  "BAF-C2 Flow Valve (%)",
                  "BAF-C3 Flow Valve (%)",
                  "BAF-C4 Flow Valve (%)",
                  "UF-P1 (%)",
                  "UF-J2 (%)",
                  "UF-J3 (%)",
                  "INF-TOC (ppm)",
                  "INF-TOC (mg/L)",
                  "INF-TIC (ppm)",
                  "INF-Ammonia (mg/L)",
                  "INF-Nitrate (mg/L)",
                  "INF-NOX (ppmv)",
                  "INF-pH",
                  "INF-Dissolved Oxygen (mg/L)",
                  "BAF-Ammonia (mg/L)",
                  "BAF-Nitrate (mg/L)",
                  "BAF-NOX (ppmv)",
                  "BAF-pH",
                  "BAF-Dissolved Oxygen (mg/L)",
                  "O3-Ozone State",
                  "UF-Backwash Pump State",
                  "UF-Air Scour Valve State",
                  "UF-Purge Valve State",
                  "UF-Permeate Drain/BW Fill Valve State",
                  "RealTech #2 Dual Feed State",
                  "fault_occurring",
                  "fault_section")

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
            xmin = min(Date_Time), xmax = max(Date_Time),
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
    geom_vline(xintercept = end_training, linetype = "dashed", color = "red") +
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

# For use with plot_ts
i_split <- str_which(vars_relevant, "O3-Ozone State") - 1

vars_exp <- tibble(var = vars_relevant,
                   type = c(rep("line", i_split),
                            rep("point", length(vars_relevant) - i_split)))

width_1 <- 4000
height_1 <- 5200
height_2 <- 2500


# Experiment 3: No Air Scour During Backwash ------------------------------
df3 <- readRDS("data_dpr/dpr_fault_exp_3.rds")

df3 <- df3 |> 
  mutate(`TMP (psi)`= `UF-Feed Pressure (psi)` - `UF-Permeate Pressure (psi)`,
         fault_section = get_fault_sections(df3)) |> 
  select(all_of(vars_relevant))

time_exp_3 <- c(ymd_hms("2023-06-16 19:00:00 UTC"),
                ymd_hms("2023-06-22 10:00:00 UTC"))

df3_reduced <- df3 |> filter(between(Date_Time, time_exp_3[1], time_exp_3[2]))

end_training <- df3_reduced$Date_Time[floor(nrow(
    df3_reduced[1:(min(which(df3_reduced$fault_occurring == 1)) - 1), ])/2
  )]

# Indicate variables to visualize
vars_relevant_exp_3 <- c("O3-Gas Feed Concentration (mg/L)",
                         "BAF-Turbidity [2] (NTU)",
                         "UF-Feed Flow (GPM)",
                         "TMP (psi)",
                         "UF-Permeate Flow (GPM)",
                         "UF-Reject Pressure (psi)",
                         "UF-Current Draw (A)",
                         "UF-Turbidity [3] (NTU)")

vars_exp_3 <- filter(vars_exp, var %in% vars_relevant_exp_3) |> 
  mutate(var = factor(var, levels = vars_relevant_exp_3)) |> 
  arrange(var) |> 
  mutate(var = as.character(var), y_min = NA, y_max = NA)

vars_exp_3[2, 3] <- 0
vars_exp_3[2, 4] <- 4
vars_exp_3[4, 3] <- -10
vars_exp_3[4, 4] <- 10
vars_exp_3[6, 3] <- 0
vars_exp_3[6, 4] <- 10
vars_exp_3[7, 3] <- 0
vars_exp_3[7, 4] <- 10
vars_exp_3[8, 3] <- 0
vars_exp_3[8, 4] <- 1

# Plot trimmed data -------------------------------------------------------

# All variables
wrap_plots(pmap(vars_exp_3,\(var, type, y_min, y_max) {
  plot_ts(df3_reduced, var, type, y_min, y_max)
}), guides = 'collect', ncol = 1)

ggsave("figures/ts_fault_trimmed.png", width = width_1, height = height_1, units = "px")

# First 4
wrap_plots(pmap(vars_exp_3[1:4, ],\(var, type, y_min, y_max) {
  plot_ts(df3_reduced, var, type, y_min, y_max)
}), guides = 'collect', ncol = 1)

ggsave("figures/ts_fault_trimmed_1.png", width = width_1, height = height_2, units = "px")

# Second 4
wrap_plots(pmap(vars_exp_3[5:8, ],\(var, type, y_min, y_max) {
  plot_ts(df3_reduced, var, type, y_min, y_max)
}), guides = 'collect', ncol = 1)

ggsave("figures/ts_fault_trimmed_2.png", width = width_1, height = height_2, units = "px")
