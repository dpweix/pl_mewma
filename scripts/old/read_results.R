library("tidyverse")
library("here")
source(here("scripts", "pl_fd_methods.R"))


# Read sim study results --------------------------------------------------
results_hawkins <- 
  map_dfr(1:5, \(i) {
    as_tibble(readRDS(here("data", paste0("hawkins-", i, ".rds")))$rl)
  })


map(results_hawkins, \(x) x$rl$s1)

as_tibble(results_hawkins[[1]]$rl)
