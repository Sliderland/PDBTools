library(cmdstanr)
install_cmdstan(version = "2.32.2", overwrite = TRUE)
set_cmdstan_path("/Users/gerpr308/.cmdstan/cmdstan-2.32.2")
folder_path <- "HeapsStanPrograms/"
source("helper.R")
model_names <- c(
  "semiconj.stan",
  "statinvertprior.stan",
  "statprior.stan",
  "statrmlprior.stan"
)
file_locations <- sapply(model_names, function(x) paste0(folder_path, x))
stan_models <- lapply(file_locations, function(x) cmdstan_model(x, compile = FALSE))
lapply(stan_models,
  function(x) x$format(canonicalize = list("deprecations"),
  overwrite_file = TRUE))
