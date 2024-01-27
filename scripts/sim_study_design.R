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
      
      for(m in methods) {
        method <<- m
        rstudioapi::jobRunScript(path = here("scripts", "sim_study_single_run.R"),
                                 name = paste0(scenario, "-", method, "-", arg, "-", i),
                                 importEnv = TRUE)
      }

    }
    
    files_to_save <- map(methods, \(m) {
      paste0(scenario, "-", m, "-", scenario_params[[scenario]], "-", i, ".rds")}) |> 
      unlist()
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
