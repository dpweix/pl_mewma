library("tidyverse")
library("here")

# Run selected method -----------------------------------------------------
source(here("scripts", "sim_study_settings.R"))


for(scenario in names(scenario_params)) {
  scenario <<- scenario
  
  for(i in i_min:i_max) {
    i <<- i
    for(arg in scenario_params[[scenario]]) {
      arg <<- arg
      rstudioapi::jobRunScript(path = here("scripts", "sim_study_single_run.R"),
                               name = paste0(scenario, "-", method, "-", arg, "-", i),
                               importEnv = TRUE)
      
    }
    
    files_to_save <- paste0(scenario, "-", method, "-", scenario_params[[scenario]], "-", i, ".rds")
    files_saved <- FALSE
    
    while(!files_saved) {
      Sys.sleep(3)
      files_saved <- map_lgl(files_to_save,
                             \(x) list.files(here(data_folder)) |>
                               str_detect(x) |> 
                               any()) |> all()
    }
  }
}


# Run selected scenario/method -----------------------------------------------
scenario <- "s1";
selected_method <- method_params[["hawkins"]]
attach(selected_method)

n <- 1000
arl_ic <- 200
i_min <- 1
i_max <- 3

for(i in i_min:i_max) {
  i <<- i
  for(arg in scenario_params[[scenario]]) {
    arg <<- arg
    rstudioapi::jobRunScript(path = here("scripts", "sim_study_single_run.R"),
                             name = paste0(scenario, "-", method, "-", arg, "-", i, ".rds"),
                             importEnv = TRUE)
    
  }
  
  files_to_save <- paste0(scenario, "-", method, "-", scenario_params[[scenario]], "-", i, ".rds")
  files_saved <- FALSE
  
  while(!files_saved) {
    Sys.sleep(3)
    files_saved <- map_lgl(files_to_save,
                           \(x) list.files(here(data_folder)) |>
                             str_detect(x) |> 
                             any()) |> all()
  }
}
